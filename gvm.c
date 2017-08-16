#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <config.h>

#include "yaml.h"

#include "uthash.h"

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/synced_bcf_reader.h"

#include "utils.h"
#include "ref_seq.h"
#include "base_seq_repr.h"
#include "nmparse.h"
#include "bedparse.h"
#include "report.h"
#include "bam_multi_itr.h"
#include "bam_mate_table.h"

// to aid with debugging in certain cases
#ifdef DEBUG
#include "bam_macro_funcs.h"
#endif
// bam_endpos from htslib is expensive and yet does not cache the immutable result.
// I've gone in and implemented this myself and called it bam_endpos2. However, I
// cannot expect everyone to use the patch so this feature test allows the program
// to fall back to the original implementation.
#ifndef HTS_HAS_ENHANCED_ENDPOS
#define bam_endpos2 bam_endpos
#endif

// The synced_bcf_reader api allows reading multiple vcfs at the same time. They
// seem to be accessed individually through indices, so I'm defining those here.
#define SNP_VCF_INDEX    0
#define COSMIC_VCF_INDEX 1

// Settings {{{

/* This program's settings */
struct settings {
	char ref_file[256];
	char bam_list[256];
	char bed_file[256];
	char snp_vcf_path[256];
	char snp_vcf_name[256];
	char cosm_vcf_path[256];
	char normal_base_path[256];
	char out_name[256];

	char chromosome[100];

	uint32_t min_mq;
	uint32_t min_bq;

	uint32_t default_mq;
	uint32_t default_bq;

	uint32_t min_variants;

	double pv_freq;
};

/* global settings object */
struct settings settings;

// }}}

// Context struct definition {{{
struct context {
	struct ref_seq ref_seq_info;
	struct bam_multi_itr *bmi;

	bcf_srs_t *bcf_reader;

	// be careful with these
	bam1_t *bam;
	bam1_t *mbam;

	uint32_t reg_start; // BED region start
	uint32_t reg_end;

	int sample_index;
	int tid;
	char *target_name;

	struct nm_itr *nmi;
	struct nm_tbl *nmt;

	FILE *pos_file;
	FILE *exon_file;
};

// }}}

struct extra_data {
	uint32_t mm_count;
	uint32_t total_softclip;
};


// Recording data from BAM {{{

/* utility */
static inline
int is_mismatch(struct alignment_report rep, struct ref_seq ref_seq_info)
{
	return !( rep.size > 0 &&
			rep.align_type == at_single &&
			rep.data == char_to_b5base(ref_seq_get(ref_seq_info, rep.pos)) );
}

static
int modify_overlap_quals(bam1_t *bam, bam1_t *mbam, 
		struct alignment_report rep,
		uint32_t *mq, uint32_t *bq,
		int *is_overlap)
{

	uint32_t p1, p2, q1, q2, fspos, rspos, fb, rb, fbq, rbq;

	*is_overlap = 0;

	bam1_t *before, *after;

	if (bam->core.pos > mbam->core.pos) {
		after = bam;
		before = mbam;
	} else {
		after = mbam;
		before = bam;
	}

	p1 = before->core.pos + 1;
	p2 = after->core.pos + 1;
	q1 = bam_endpos2(before) + 1;
	q2 = bam_endpos2(after) + 1;

	assert(p2 >= p1);

	*is_overlap = rep.pos >= p2 && rep.pos < q1 && rep.pos >= p1 && rep.pos < q2;

	if (*is_overlap) {
		assert(rep.pos >= p1 && rep.pos < q2);
		if (bam_is_rev(mbam)) {
			if (rep.align_type == at_single) {
				fspos = rep.spos;
				rspos = rep.pos - mbam->core.pos - 1;
				fb = bam_seqi(bam_get_seq(bam), fspos);
				rb = bam_seqi(bam_get_seq(mbam), rspos);
				fbq = bam_get_qual(bam)[fspos];
				rbq = bam_get_qual(mbam)[rspos];
				if (fb == rb) {
					*bq += rbq;
					if (rbq > fbq) rep.data = int_to_b5base(rb);
				} else {
					*bq = 0;
				}

				*mq += mbam->core.qual;
				*mq /= 2;
			} else {
				// no one knows what to do here
			}

		} else {
			return 1; // avoid double counting
		}
	}

	return 0;
}

static
void record_match(	struct context *context,
			struct alignment_report rep,
			void *extra_data_p)
{
	// ignore soft clips
	if (rep.align_type == at_sclip) return;

	struct extra_data *extra_data = extra_data_p;

	uint32_t mq, bq;
	int is_overlap = 0;

	uint32_t sample_index = context->sample_index;

	bam1_t *bam = context->bam;
	bam1_t *mbam = context->mbam;

	uint32_t mm_count = extra_data->mm_count;

	struct variant_table *vtable, *vtentry;
	struct variant_counts *vcounts;

	vtable = context->bmi->itr_list[sample_index].vtable;

	mq = bam->core.qual;
	bq = bam_get_qual(bam)[rep.spos];

	if (mbam != NULL) {
		// We are dealing with a read pair
		int double_count = modify_overlap_quals(bam, mbam, rep, &mq, &bq, &is_overlap);

		if (double_count) return;
	}

	/* Check if the offset already has an entry */
	HASH_FIND_INT(vtable, &rep.pos, vtentry);

	if (vtentry == NULL) {

		/* Prepare the hash entry */
		vtentry = calloc(1, sizeof(struct variant_table));
		if (vtentry == NULL) {
			err_printf("could not allocate memory\n");
			return;
		}
		vtentry->sample_index = sample_index;
		vtentry->tid = context->tid;
		vtentry->offset = rep.pos;
		/* Add it to the hash table */

		HASH_ADD(hh, vtable, offset, sizeof(rep.pos), vtentry);
	}
	
	vtentry->read_count++;

	if (mq >= settings.min_mq && bq >= settings.min_bq) {

		vtentry->read_count_pass++;

		/* check if this is a mismatch */
		if (is_mismatch(rep, context->ref_seq_info)) {
			vtentry->mismatch_count++;
		}// else printf("\n");

		/* Update the entry */
		struct variant_counts *vcounts_table = vtentry->counts;
		
		/* We use a zero'd out key holder to avoid any struct padding nonsense */
		struct variant_counts key_holder;
		memset(&key_holder, 0, sizeof(key_holder));
		/* These are the only three values we're interested in */
		key_holder.report.align_type = rep.align_type;
		key_holder.report.pos = rep.pos;
		key_holder.report.data = rep.data;
		key_holder.report.size = rep.size;
		
		HASH_FIND(hh, vcounts_table, &key_holder.report, sizeof(rep), vcounts);
		if (vcounts == NULL) {
			vcounts = calloc(1, sizeof(struct variant_counts));
			if (vtentry == NULL) {
				err_printf("could not allocate memory\n");
				return;
			}

			vcounts->pos = rep.pos;
			vcounts->report = key_holder.report;

			HASH_ADD(hh, vcounts_table, report, sizeof(rep), vcounts);
		}// else { err_printf("HASH HIT!\n"); }
	
		/* Update the count */
		if (is_overlap) {
			vcounts->count_r++;
			vcounts->count_f++;
		} else if (bam_is_rev(bam)) {
			vcounts->count_r+=2;
		} else {
			vcounts->count_f+=2;
		}

		/* Update the total MQs and BQs */
		vcounts->total_mq += mq;
		vcounts->total_bq += bq;
		vcounts->total_read_pos += rep.spos;
		vcounts->total_pmm += (double) mm_count / bam->core.l_qseq;

		vtentry->counts = vcounts_table;

	}


	/* Update the vtable */
	context->bmi->itr_list[sample_index].vtable = vtable;

}

static
void report_aggregate(	struct context *context,
			struct alignment_report rep, 
			void *extra_data_p)
{
	struct extra_data *extra_data = extra_data_p;
	struct ref_seq ref_seq_info = context->ref_seq_info;

	/* check if this is a mismatch */
	if (is_mismatch(rep, ref_seq_info)) {
		extra_data->mm_count += rep.size;
	}

	if (rep.align_type == at_sclip) {
		extra_data->total_softclip += rep.size;
	}
}


// }}}

// CIGAR string processor {{{
// returns negative on error, otherwise if not
static
int run_cigar(	struct context *context,
		uint32_t cigar_op, uint32_t cigar_nextop,
		uint32_t *base_offset,
		uint32_t *seq_index,
		reporter_func rfunc,
		void *extra_data)
{

	bam1_t *bam = context->bam;
	uint32_t op_type, len, next_op_type;
	uint8_t *seq_data = bam_get_seq(bam);

	struct alignment_report report;
	memset(&report, 0, sizeof(report));


	op_type = bam_cigar_op(cigar_op);
	len = bam_cigar_oplen(cigar_op);
	next_op_type = bam_cigar_op(cigar_nextop);

	switch(op_type) {
	case BAM_CDEL:
		report.align_type = at_del;
		report.pos = *base_offset - 1;
		report.spos = *seq_index - 1;
		report.data = int_to_b5base(bam_seqi(seq_data, *seq_index));
		report.size = len;
		rfunc(context, report, extra_data);
		*base_offset += len;
		break;
	case BAM_CINS:
		report.align_type = at_ins;
		report.pos = *base_offset - 1;
		report.size = len;

		if (*seq_index == 0) {
			report.spos = 0;
			report.data = char_to_b5base(ref_seq_get(context->ref_seq_info, report.pos));
		} else {
			report.spos = *seq_index - 1;
			report.data = int_to_b5base(bam_seqi(seq_data, report.spos));
		}

		if (len + 1 > INDEL_REPR_LEN_CUTOFF) {
			report.data = INDEL_REPR_CUTOFF + len + 1;
		} else {
			report.data = collect_seqi(seq_data, report.spos + 1, len, report.data);
		}

		rfunc(context, report, extra_data);
		*seq_index += len;
		break;
	case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
		report.align_type = at_single;
		while (len-- > 0) {
			report.pos = *base_offset;
			report.spos = *seq_index;
			report.data = int_to_b5base(bam_seqi(seq_data, *seq_index));
			report.size = 1;
			if (len == 0 && (next_op_type == BAM_CINS || next_op_type == BAM_CDEL)) {
				// skip
			} else {
				rfunc(context, report, extra_data);
			}
			(*base_offset)++;
			(*seq_index)++;
		}
		break;
	case BAM_CSOFT_CLIP:
		*seq_index += len;
		//*base_offset += len+2;
		break;
	case BAM_CHARD_CLIP:
		*base_offset += len;
		break;
	case BAM_CBACK:
		// I can't find any documentation on this cigar op so I'm
		// just going to report it as an error.
		return -1;
	}

	return 0;

}

// }}}

// Alignment calculation {{{

/* This function is passed ONE alignment and reports each allele with its
 * position to the reporter function.
 */
static
int calc_alignments(struct context *context, reporter_func rfunc)
{
	bam1_t *bam = context->bam;

	uint32_t offset = bam->core.pos + 1;

	uint32_t cigar_len = bam->core.n_cigar;
	if (cigar_len == 0) return 0;

	uint32_t *cigar_data = (uint32_t *) bam_get_cigar(bam);
	uint32_t *cigar_ptr;
	uint32_t cigar_nextop;

	uint32_t seq_index = 0;

	int result = 0;

	uint32_t i;

	struct extra_data extra_data;

	uint32_t offset_backup, seq_index_backup;

	for (cigar_ptr = cigar_data, i = 0; i < cigar_len; i++, cigar_ptr++) {
		/* I need to run the run_cigar function twice- first time to calculate the
		 * pmm, next time to actually report all the matches. Therefore, these two
		 * values need to be backed up */

		if (i < cigar_len - 1) {
			cigar_nextop = *(cigar_ptr + 1);
		} else {
			cigar_nextop = 0;
		}

		extra_data.mm_count = 0;
		extra_data.total_softclip = 0;

		offset_backup = offset;
		seq_index_backup = seq_index;

		result = run_cigar(
			context,
			*cigar_ptr, cigar_nextop,
			&offset,          /* modified, needs backup */
			&seq_index,       /* modified, needs backup */
			report_aggregate, /* Alternate report function */
			&extra_data       /* modified, no backup */
		);

		/* Restore the backups */
		seq_index = seq_index_backup;
		offset = offset_backup;

		result = run_cigar(
			context,
			*cigar_ptr, cigar_nextop,
			&offset, /* base */
			&seq_index,
			rfunc,
			&extra_data
		);



		if (result == -1) {
			err_printf("An error was encountered while running the cigar string\n");
			break;
		}
	}

	return result;
}
// }}}

// Result reporting {{{

// TODO: do this properly
static
int chr2idx(const char *chrname)
{
	if (strlen(chrname) == 0) return -1;		
	int sum = chrname[0] + 2*chrname[1];
	switch (sum) {
	case '1': return 0;
	case '2': return 1;
	case '3': return 2;
	case '4': return 3;
	case '5': return 4;
	case '6': return 5;
	case '7': return 6;
	case '8': return 7;
	case '9': return 8;
	case '1' + 2*'0': return 9;
	case '1' + 2*'1': return 10;
	case '1' + 2*'2': return 11;
	case '1' + 2*'3': return 12;
	case '1' + 2*'4': return 13;
	case '1' + 2*'5': return 14;
	case '1' + 2*'6': return 15;
	case '1' + 2*'7': return 16;
	case '1' + 2*'8': return 17;
	case '1' + 2*'9': return 18;
	case '2' + 2*'0': return 19;
	case '2' + 2*'1': return 20;
	case '2' + 2*'2': return 21;
	case 'X': return 22;
	case 'Y': return 23;
	}

	return -1;
}

static
int cmp_vcounts(struct variant_counts *a, struct variant_counts *b)
{
	return (a->count_f + a->count_r) - (b->count_f + b->count_r);
}

static
void sort_vcounts(struct variant_table *vtable)
{
	HASH_SORT(vtable->counts, cmp_vcounts);
}

static
void find_ab(struct variant_table *vtable, struct variant_counts **a, struct variant_counts **b)
{
	struct variant_counts *vcounts = vtable->counts, *candidate, *tmp;

	uint32_t biggest = 0, second_biggest = 0;

	if (vtable->a_cache != NULL && vtable->b_cache != NULL) {
		*a = vtable->a_cache;
		*b = vtable->b_cache;
		return;
	}

	*a = NULL;
	*b = NULL;

	HASH_ITER(hh, vcounts, candidate, tmp) {
		uint32_t total = candidate->count_f + candidate->count_r;
		if (total > biggest) {
			second_biggest = biggest;
			biggest = total;
			*b = *a;
			*a = candidate;
		} else if (total > second_biggest) {
			second_biggest = total;
			*b = candidate;
		}
	}

	vtable->a_cache = *a;
	vtable->b_cache = *b;
}

static
int get_af(bcf_sr_t *reader, bcf1_t *line, float **af)
{
	// The call to bcf_get_info_float clobbers the stack for whatever reason!
	// To avoid this, I use the heap as the destination of the received value.
	float *af2 = malloc(sizeof(float));

	int result, af_size = sizeof(**af);
	result = bcf_get_info_float(
			reader->header,
			line,
			"AF",
			&af2,
			&af_size);


	if (result < 1) {
		*af = NULL;
	} else {
		**af = *af2;
	}

	free(af2);

	return result;

}

static
void get_afs(	struct context *context,
		struct alignment_report report,
		float **pop_af, float **cosm_af)
{

	bcf1_t *line;
	bcf_sr_t *reader;

	bcf_sr_seek(context->bcf_reader, context->target_name, report.pos);
	bcf_sr_next_line(context->bcf_reader);

	if (bcf_sr_has_line(context->bcf_reader, SNP_VCF_INDEX)) {
		line = bcf_sr_get_line(context->bcf_reader, SNP_VCF_INDEX);
		if ((uint32_t) line->pos == report.pos) {
			reader = &context->bcf_reader->readers[SNP_VCF_INDEX];
			get_af(reader, line, pop_af);
		} else {
			*pop_af = NULL;
		}
	} else {
		*pop_af = NULL;
	}

	if (bcf_sr_has_line(context->bcf_reader, COSMIC_VCF_INDEX)) {
		line = bcf_sr_get_line(context->bcf_reader, COSMIC_VCF_INDEX);
		if ((uint32_t) line->pos == report.pos) {
			reader = &context->bcf_reader->readers[COSMIC_VCF_INDEX];
			get_af(reader, line, cosm_af);
		} else {
			*cosm_af = NULL;
		}
	} else {
		*cosm_af = NULL;
	}

}

static
float get_af_true(float *af_p, int is_mm)
{
	if (af_p == NULL) {
		return is_mm ? settings.pv_freq : 1 - 3*settings.pv_freq;
	} else {
		return is_mm ? *af_p : 1 - *af_p;
	}

}

static
void dump_variant_info(	struct context *context,	
			FILE *f,
			struct variant_counts *vc,
			uint32_t ref_allele_partial)
{

	struct alignment_report report = vc->report;

	double total = (vc->count_f + vc->count_r) / 2;
	if (total == 0) total = 1;

	uint32_t read = 0;
	switch(report.align_type) {
	case at_single: case at_ins:
		read = seq_append2(report.data, ref_allele_partial);
		break;
	case at_del:
		read = char_to_b5base(ref_seq_get(context->ref_seq_info, report.pos));
		break;
	case at_sclip:
		assert(0 /* soft clip should never be reported! */);
		break;
	}

#ifdef DEBUG
#define VARIANT_TYPE_FORMAT  "%c="
#define VARIANT_TYPE	     report.align_type["SID"]
#else
#define VARIANT_TYPE_FORMAT  "%s"
#define VARIANT_TYPE	     ""
#endif
	fprintf(f, VARIANT_TYPE_FORMAT "%d\t%g\t%g\t%g\t%g\t%g\t%g",
			VARIANT_TYPE,
			read,
			(double) vc->count_f / 2,
			(double) vc->count_r / 2,
			(double) vc->total_bq / total,
			(double) vc->total_mq / total,
			vc->total_pmm / total,
			(double) vc->total_read_pos / total);
#undef VARIANT_TYPE_FORMAT
#undef VARIANT_TYPE
}

static
void dump_nm_data(	struct context *context,
			uint32_t offset )
{
	FILE *f = context->pos_file;

	struct nm_tbl *nmt = context->nmt;
	struct nm_entry *ent = nm_tbl_get(nmt, offset);

	if (ent == NULL) {
		err_printf("no normal metrics data for chromosome %s, offset %d\n",
				settings.chromosome, offset);
		fprintf(f, "NaN\tNaN\tNaN\tNaN");
	} else {

		fprintf(f, "%g\t%g\t%g\t%g",
				ent->norm_read_depth,
				ent->prob_map_err,
				ent->read_pass,
				ent->a_or_b);
	}
}

static
void dump_vcounts(	struct context *context,
			struct variant_table *v,
			int max_delete_size )
{
	FILE *f = context->pos_file;
	struct variant_counts *a = NULL, *b = NULL, dummy;
	struct ref_seq ref_seq_info = context->ref_seq_info;

	int is_mm;

	uint32_t ref_allele, ref_allele_partial;

	float pop_af,
	      cosm_af,
	      a_pop_af,
	      b_pop_af;

	float *pop_af_p = &pop_af,
	      *cosm_af_p  = &cosm_af;

	find_ab(v, &a, &b);


	if (!a || !b) {
		memset(&dummy, 0, sizeof(dummy));
		dummy.total_mq = settings.default_mq;
		dummy.total_bq = settings.default_bq;
	}

	if (!a) a = &dummy;
	if (!b) b = &dummy;

	ref_allele_partial =
		collect_seqc(ref_seq_info, v->offset + 1, max_delete_size, 0);

	ref_allele = seq_append2( char_to_b5base(ref_seq_get(ref_seq_info, v->offset)),
				  ref_allele_partial );

	fprintf(f, "%d\t%d\t%d\t%d\t%d\t%d\t",
		v->sample_index + 1,
		v->tid + 1, /* chromosomes are 0-indexed apparently */
		v->offset,
		v->read_count,
		v->read_count_pass,
		ref_allele);

	dump_variant_info(context, f, a, ref_allele_partial);
	fprintf(f, "\t");
	dump_variant_info(context, f, b, ref_allele_partial);
	fprintf(f, "\t");

	get_afs(context, a->report, &pop_af_p, &cosm_af_p);

	is_mm = is_mismatch(a->report, context->ref_seq_info);
	a_pop_af = get_af_true(pop_af_p, is_mm);

	is_mm = is_mismatch(b->report, context->ref_seq_info);
	b_pop_af = get_af_true(pop_af_p, is_mm);

	if (cosm_af_p == NULL) {
		cosm_af = 0;
	}

	fprintf(f, "%g\t%g\t%g\t", a_pop_af, b_pop_af, cosm_af);

	dump_nm_data(context, v->offset);
	fprintf(f, "\n");
}

static
void dump_blank_vcounts(struct context *context, uint32_t sample_idx, uint32_t offset, int max_delete_size)
{
	int tid = chr2idx(settings.chromosome);
	struct variant_table vt = {0};

	vt.sample_index = sample_idx;
	vt.offset = offset;
	vt.tid = tid;

	dump_vcounts(context, &vt, max_delete_size);
	
}

static inline
int get_delete_size(struct variant_counts *vc)
{
	if (vc == NULL) return -1;

	return vc->report.align_type == at_del ? vc->report.size : 0;
}

static
int get_max_delete_size(struct context *context, uint32_t offset)
{
	uint32_t sample_idx;
	struct variant_table *vtable;
	struct variant_counts *counts, *tmp, *a, *b;

	int max_delete_size = 0, max_of_ab;
	int has_enough = 0;

	for (sample_idx = 0; sample_idx < context->bmi->num_iters; sample_idx++) {
		HASH_FIND_INT(context->bmi->itr_list[sample_idx].vtable, &offset, vtable);
		if (vtable == NULL) continue;


		find_ab(vtable, &a, &b);
		max_of_ab = MAX(get_delete_size(a), get_delete_size(b));
		if (max_of_ab > max_delete_size) max_delete_size = max_of_ab;

		HASH_ITER(hh, vtable->counts, counts, tmp) {
			if (is_mismatch(counts->report, context->ref_seq_info) &&
					counts->count_f + counts->count_r >=
						2 * settings.min_variants) {
				has_enough = 1;
			}
		}

	//	HASH_ITER(hh, vtable->counts, counts, tmp) {
	//		// Remember, count_f + count_r is TWICE the actual count
	//		is_enough =
	//			counts->count_f + counts->count_r >= 2 * settings.min_variants;

	//		if (!is_enough) continue;

	//		cur_delete_size = get_delete_size(counts);

	//		if (cur_delete_size > max_delete_size) {
	//			max_delete_size = cur_delete_size;
	//		}
	//		
	//		if (is_mismatch(counts->report, context->ref_seq_info)) {
	//			has_enough = 1;
	//		}
	//	}
	}

	return has_enough ? max_delete_size : -1;
}

static
void flush_results(struct context *context, uint32_t begin, uint32_t end)
{
	int max_delete_size;
	struct variant_table *vt, *v;

	uint32_t offset, sample_idx;

	// clamp in region
	if (begin < context->reg_start) begin = context->reg_start;
	if (end > context->reg_end + 1) end = context->reg_end + 1;

	for (offset = begin; offset < end; offset++) {
		max_delete_size = get_max_delete_size(context, offset);
		if (max_delete_size < 0) continue;

		for (sample_idx = 0; sample_idx < context->bmi->num_iters; sample_idx++) {
			vt = context->bmi->itr_list[sample_idx].vtable;

			HASH_FIND_INT(vt, &offset, v);
			if (v == NULL) {
				dump_blank_vcounts(context, sample_idx, offset, max_delete_size);
			} else {
				dump_vcounts(context, v, max_delete_size);
			}
		}
	}
}

// }}}

// Checking various properties of alignments {{{
#define BAM_FLAG_GOODALIGN 0x2
#define BAM_FLAG_REVCOMP 0x10
#define BAM_SECONDARY 0x100
#define BAM_FLAG_DUPLICATE 0x400
#define BAM_SUPPLEMENTARY 0x800

static
int check_flags(bam1_t *bam)
{
	uint16_t flagval = bam->core.flag;

	return (flagval & BAM_FLAG_GOODALIGN) &&
		!(flagval & BAM_FLAG_DUPLICATE) &&
		!(flagval & BAM_SECONDARY) &&
		!(flagval & BAM_SUPPLEMENTARY);
}

static
int is_split_read(bam1_t *bam)
{
	return bam->core.tid != bam->core.mtid;
}

// }}}

// Configuration reading {{{
// Many thanks to https://www.wpsoftware.net/andrew/pages/libyaml.html

int config_try_set_option(char *option, char *value, struct settings *s)
{
	uint32_t min_mq, min_bq, default_mq, default_bq, min_variants;
	double pv_freq;
	char *tmp;

	if (strcmp(option, "bamList") == 0) {
		strncpy(s->bam_list, value, sizeof(s->bam_list));
		return 1;
	}

	if (strcmp(option, "refGenome") == 0) {
		strncpy(s->ref_file, value, sizeof(s->ref_file));
		return 1;
	}

	if (strcmp(option, "regionsFile") == 0) {
		strncpy(s->bed_file, value, sizeof(s->bed_file));
		return 1;
	}

	if (strcmp(option, "snpVCFpath") == 0) {
		strncpy(s->snp_vcf_path, value, sizeof(s->snp_vcf_path));
		return 1;
	}

	if (strcmp(option, "snpVCFname") == 0) {
		strncpy(s->snp_vcf_name, value, sizeof(s->snp_vcf_name));
		return 1;
	}

	if (strcmp(option, "cosmicVCF") == 0) {
		strncpy(s->cosm_vcf_path, value, sizeof(s->cosm_vcf_path));
		return 1;
	}

	if (strcmp(option, "NormalBase") == 0) {
		strncpy(s->normal_base_path, value, sizeof(s->normal_base_path));
		return 1;
	}

	if (strcmp(option, "outName") == 0) {
		strncpy(s->out_name, value, sizeof(s->out_name));
		return 1;
	}

	if (strcmp(option, "minBQ") == 0) {
		min_bq = strtol(value, &tmp, 10);
		if (*tmp != '\0') {
			err_printf("Invalid value for minBQ: %s\n", value);
			return 0;
		}

		s->min_bq = min_bq;
		return 1;
	}

	if (strcmp(option, "minMQ") == 0) {
		min_mq = strtol(value, &tmp, 10);
		if (*tmp != '\0') {
			err_printf("Invalid value for minMQ: %s\n", value);
			return 0;
		}

		s->min_mq = min_mq;
		return 1;
	}

	if (strcmp(option, "defaultBQ") == 0) {
		default_bq = strtol(value, &tmp, 10);
		if (*tmp != '\0') {
			err_printf("Invalid value for defaultBQ: %s\n", value);
			return 0;
		}

		s->default_bq = default_bq;
		return 1;
	}

	if (strcmp(option, "defaultMQ") == 0) {
		default_mq = strtol(value, &tmp, 10);
		if (*tmp != '\0') {
			err_printf("Invalid value for defaultMQ: %s\n", value);
			return 0;
		}

		s->default_mq = default_mq;
		return 1;
	}

	if (strcmp(option, "minBCount") == 0) {
		min_variants = strtol(value, &tmp, 10);
		if (*tmp != '\0') {
			err_printf("Invalid value for minBCount: %s\n", value);
			return 0;
		}

		s->min_variants = min_variants;
		return 1;
	}

	if (strcmp(option, "pvFreq") == 0) {
		pv_freq = strtod(value, &tmp);
		if (*tmp != '\0') {
			err_printf("Invalid value for pvFreq: %s\n", value);
			return 0;
		}

		s->pv_freq = pv_freq;
		return 1;
	}

	return 1;

}

int config_read(const char *cfg_fn, struct settings *s)
{
	FILE *f = fopen(cfg_fn, "r");
	int result = 1;

	/* FSM */
	int state = 0; /* 0 = need key, 1 = need value */
	char value_buf[1024];

	if (f == NULL) {
		err_printf("could not open config %s: ", cfg_fn);
		perror("");
		return 0;
	}


	yaml_parser_t parser;
	yaml_event_t event;

	if (!yaml_parser_initialize(&parser)) {
		err_printf("failed to initialize YAML parser\n");
		result = 0;
		goto done;
	}
	
	yaml_parser_set_input_file(&parser, f);

	/* START new code */
	do {
		if (!yaml_parser_parse(&parser, &event)) {
			err_printf("Parser error %d\n", parser.error);
			result = 0;
			goto done1;
		}

		switch(event.type)
		{ 
		case YAML_NO_EVENT:
		case YAML_STREAM_START_EVENT:
		case YAML_STREAM_END_EVENT:
		case YAML_DOCUMENT_START_EVENT: 
		case YAML_DOCUMENT_END_EVENT:	
		case YAML_SEQUENCE_START_EVENT: 
		case YAML_SEQUENCE_END_EVENT:	
		case YAML_MAPPING_START_EVENT:	
		case YAML_MAPPING_END_EVENT:	
			break;

		case YAML_ALIAS_EVENT:
			fprintf(stderr, "Warning: no alias support (anchor %s)\n", event.data.alias.anchor);
			break;
		case YAML_SCALAR_EVENT: 
			if (state == 0) { /* This means we need a key */
				strncpy(value_buf, (char *)event.data.scalar.value, sizeof(value_buf));
				state = 1;
			} else { /* We have a key, now we need a value */
				if (!config_try_set_option(value_buf, (char *) event.data.scalar.value, s)) {
					result = 0;
					goto done1;
				}
				state = 0;
			}

		default: break;

		}


		if (event.type != YAML_STREAM_END_EVENT) {
			yaml_event_delete(&event);
		}
	} while (event.type != YAML_STREAM_END_EVENT);

done1:
	yaml_event_delete(&event);
	yaml_parser_delete(&parser);

done:
	fclose(f);
	return result;
} // }}}

static
void write_exon_line(struct context *context, uint32_t start, uint32_t end)
{
	int i;
	uint32_t offset;
	double mean_read_depth, count;
	FILE *exon_file;
	struct variant_table *v;
	struct bam_multi_itr *bmi;

	struct nm_entry avgs;

	exon_file = context->exon_file;
	bmi = context->bmi;
	avgs = context->nmt->avgs;


	for (i = 0; i < bmi->num_iters; i++) {
		mean_read_depth = 0;
		count = 0.0;
		for (offset = start; offset <= end; offset++) {
			HASH_FIND_INT(bmi->itr_list[i].vtable, &offset, v);
			if (v == NULL) continue;

			count++;

			mean_read_depth =
				mean_read_depth * (count - 1) / count +
				v->read_count_pass / count;
			
		}

		if (i > 0) fprintf(exon_file, "\t");
		fprintf(exon_file, "%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g",
				chr2idx(settings.chromosome) + 1,
				start,
				end,
				mean_read_depth,
				avgs.norm_read_depth,
				avgs.prob_map_err,
				avgs.read_pass,
				avgs.a_or_b);
	}

	fprintf(exon_file, "\n");
}

static
int do_region(struct context *context, uint32_t start, uint32_t end) // {{{
{
	verbose_fprintf(stderr, "--- region:  %d  to  %d  \n", start, end);

	struct bam_multi_itr *bmi = context->bmi;
	if (!bmi_query(bmi, settings.chromosome, start, end)) {
		return 0;
	}

	int result;

	struct bam_mate_table *bmt = NULL, *left_behind, *tmp;
	int32_t sample_index;

	bam1_t *bam, *mbam;

	context->nmt = nm_tbl_create();

	context->reg_start = start;
	context->reg_end = end;

	nm_query(context->nmi, start - 1, end);
	nm_tbl_slurp(context->nmt, context->nmi);

	// Main processing loop {{{
	while ( (sample_index = bmi_next(bmi, &bam)) >= 0 ) {
		if (!check_flags(bam)) continue;

		mbam = NULL;

		if (!is_split_read(bam)) {
			bmt = bmt_register(bmt, bam, &mbam, &result);
			if (result == 1) continue;
		}

		context->sample_index = (uint32_t) sample_index;
		context->bam = mbam;
		context->mbam = bam;

#define GVM_CHECK_RESULT(R) \
			if (R < 0) { \
				err_printf("failed to calculate alignments (err code %d)\n", R); \
				return EXIT_FAILURE; \
			}

		if (mbam != NULL) {
			result = calc_alignments(context, record_match);
			GVM_CHECK_RESULT(result);
		}

		context->bam = bam;
		context->mbam = mbam;

		result = calc_alignments(context, record_match);
		GVM_CHECK_RESULT(result);

		if (mbam != NULL) {
			bam_destroy1(mbam);
		}

	}
	// }}}

	context->mbam = NULL;
	HASH_ITER(hh, bmt, left_behind, tmp) {
		context->bam = left_behind->bam;
		result = calc_alignments(context, record_match);
		GVM_CHECK_RESULT(result);
	}

	flush_results(context, start, end + 1);
	// It should be empty now!
	
	write_exon_line(context, start, end);

#undef GVM_CHECK_RESULT

	bmt_destroy(bmt);
	nm_tbl_destroy(context->nmt);

	nm_cleanup(context->nmi);	
	bmi_cleanup(bmi);

	//fflush(context->pos_file);

	return 1;
} // }}}

int main(int argc, char *argv[]) // {{{
{
	faidx_t *ref_idx;
	char *ref_seq_data, *snp_vcf_fname;
	int len, result, tid;
	FILE *pos_file, *exon_file;
	struct ref_seq ref_seq_info;
	struct bam_multi_itr *bmi;
	struct nm_itr *nmi;
	struct context context;

	// Command line settings {{{
	if (argc < 3) {
		fprintf(stderr, PACKAGE_STRING "\n");
		fprintf(stderr, "usage: %s <config.yaml> <region>\n", argv[0]);
		return EXIT_FAILURE;
	}

	if (!config_read(argv[1], &settings)) {
		return EXIT_FAILURE;
	}
	strncpy(settings.chromosome, argv[2], sizeof(settings.chromosome));

	// }}}

	// BAM file loading {{{
	bmi = bmi_create(settings.bam_list);
	if (bmi == NULL) {
		return EXIT_FAILURE;
	}

	assert(bmi->num_iters > 0);
	tid = bmi_get_tid(bmi, settings.chromosome);
	// }}}

	// Reference genome loading {{{
	/* open reference sequence index */
	verbose_fprintf(stderr, "Reading reference genome index...\n");

	ref_idx = fai_load(settings.ref_file);

	if (ref_idx == NULL) {
		err_printf("Failed to load reference genome index.\n");
		return EXIT_FAILURE;
	}
	
	verbose_fprintf(stderr, "Reading reference genome...\n");

	ref_seq_data = fai_fetch(ref_idx, settings.chromosome, &len);
	ref_seq_info.data = ref_seq_data;
	ref_seq_info.len = len;
	ref_seq_info.offset = 1;

	if (len < 0) {
		err_printf("Failed to load reference genome from index. (err code %d)\n", len);
		return EXIT_FAILURE;
	}
	// }}}

	// VCF loading {{{
	context.bcf_reader = bcf_sr_init();
	snp_vcf_fname = malloc(strlen(settings.snp_vcf_path) +
				strlen(settings.snp_vcf_name) +
				strlen(settings.chromosome) + 1);
	
	if (snp_vcf_fname == NULL) {
		err_printf("failed to allocate memory\n");
		return EXIT_FAILURE;
	}

	strcpy(snp_vcf_fname, settings.snp_vcf_path);
	strcat(snp_vcf_fname, settings.chromosome);
	strcat(snp_vcf_fname, settings.snp_vcf_name);

	bcf_sr_set_regions(context.bcf_reader, settings.chromosome, 0);

	result = bcf_sr_add_reader(context.bcf_reader, snp_vcf_fname);
	if (result == 0) {
		err_printf("failed to load SNP VCF: %s\n", snp_vcf_fname);
		err_printf("errnum = %d\n", context.bcf_reader->errnum);
		return EXIT_FAILURE;
	}

	result = bcf_sr_add_reader(context.bcf_reader, settings.cosm_vcf_path);
	if (result == 0) {
		err_printf("failed to load cosmic VCF: %s\n", settings.cosm_vcf_path);
		err_printf("errnum = %d\n", context.bcf_reader->errnum);
		return EXIT_FAILURE;
	}

	free(snp_vcf_fname);
	// }}}

	// Normal metrics loading {{{
	char *normal_path = malloc(strlen(settings.normal_base_path) +
				   strlen(settings.chromosome) +
				   strlen(".txt.gz") + 1);
				 
	strcpy(normal_path, settings.normal_base_path);
	strcat(normal_path, settings.chromosome);
	strcat(normal_path, ".txt.gz");

	nmi = nm_open(normal_path, 0);

	free(normal_path);


	// }}}
	
	// Opening output files {{{
	char *pos_fn = malloc(strlen(settings.out_name) +
			      strlen(settings.chromosome) +
			      strlen("_pos.txt") + 2);
	
	if (pos_fn == NULL) {
		err_printf("unable to allocate memory.");
		return EXIT_FAILURE;
	}

	char *exon_fn = malloc(strlen(settings.out_name) +
			       strlen(settings.chromosome) +
			       strlen("_exon.txt") + 2);
	
	if (exon_fn == NULL) {
		err_printf("unable to allocate memory.");
		return EXIT_FAILURE;
	}

	sprintf(pos_fn, "%s_%s_pos.txt", settings.out_name, settings.chromosome);
	sprintf(exon_fn, "%s_%s_exon.txt", settings.out_name, settings.chromosome);

#ifdef DEBUG
	pos_file = stdout;
#else
	pos_file = fopen(pos_fn, "w");
#endif
	if (pos_file == NULL) {
		err_printf("unable to open %s", pos_fn);
		perror("");
	}

	exon_file = fopen(exon_fn, "w");
	if (exon_file == NULL) {
		err_printf("unable to open %s", exon_fn);
		perror("");
	}

	free(pos_fn);
	free(exon_fn);

	// }}}

	// Context initialization {{{
	context.ref_seq_info = ref_seq_info;
	context.bmi = bmi;
	context.tid = tid;
	context.target_name = settings.chromosome;
	context.nmi = nmi;

	context.pos_file = pos_file;
	context.exon_file = exon_file;

	if (context.tid == -1) {
		err_printf("warning: no bam files loaded\n");
	}
	// }}}
	
	bedf_forall_region_chr(
			&context,
			settings.bed_file,
		       	settings.chromosome, 
			(region_handle_func) do_region);

	// Cleanup {{{


#ifdef CLEANUP /* There's really no reason to free this right before termination */
	bmi_destroy(bmi);
	free(ref_seq_data);
	fai_destroy(ref_idx);
	bcf_sr_destroy(context.bcf_reader);
	nm_destroy(context.nmi);

	fclose(pos_file);
	fclose(exon_file);
#endif

	// }}}

	return EXIT_SUCCESS;
} // }}}
