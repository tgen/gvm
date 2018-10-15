
/* gvm.c - gvm source code
 *
 * This program is essentially a pileup engine that calculates metrics
 * for each interesting position in the genome. It is used exclusively
 * in the context of the LumosVar variant caller which uses these
 * metrics to make variant calls.
 *
 * gvm works by reading the BAM file read-by-read and applying the
 * CIGAR string to the read to determine the exact offsets of the base
 * pairs in the read.
 *
 * Information about each read is stored in a bunch of hash tables
 * (variant_table.h) which are flushed on completion of each region in
 * the genome (specified by settings.bed_file)
 *
 * To read this code, the best place to start would be the do_region
 * function. However, if you only wanted to read one function then
 * `record_match` the place to be.
 *
 * To learn how to add a new column:
 *   - understand htslib data structures:
 *     https://github.com/samtools/htslib/blob/develop/htslib/sam.h
 *   - understand gvm data structures: ./variant_table.h
 *   - read the comment above `record_match`
 *
 * To learn how to add a new YAML configuration option, read the
 * comment right above (regexp search) "^struct settings"
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"

#include <yaml.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/synced_bcf_reader.h>

#include <uthash.h>

#include "config_keys.h"
#include "utils.h"
#include "ref_seq.h"
#include "base_seq_repr.h"
#include "nmparse.h"
#include "nmcalc.h"
#include "bedparse.h"
#include "report.h"
#include "cigar.h"
#include "bam_multi_itr.h"
#include "bam_mate_table.h"
#include "cosmvcf.h"

#include "gengetopt/cmdline.h"

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

// Nice names for alignment flags
#define BAM_FLAG_GOODALIGN 0x2
#define BAM_FLAG_REVCOMP 0x10
#define BAM_SECONDARY 0x100
#define BAM_FLAG_DUPLICATE 0x400
#define BAM_SUPPLEMENTARY 0x800

// Settings {{{

/* This program's settings
 *
 * Some of these are determined by command line options, and some are
 * determined by the YAML file.
 *
 ***** How to add a YAML configuration option
 *
 * The first step to adding a new YAML configuration option is to add
 * a new member to this struct. The name of this member is the
 * "internal name" of the configuration option.
 *
 * Then, open up config_keys.h and make a #DEFINition
 *
 *    #DEFINE GVM_CONFIG_<internal name> "the option's REAL name (what's in the YAML file)"
 *
 * Then, search (in this file) for the function
 * `config_try_set_option` and add the appropriate line. If it's a
 * string option, then add a line:
 *
 * GVM_CONFIG_CHECK_STROPT(internal_name); // or NUMOPT for numerical option
 *
 * This will ensure that the value from the configuration file is put
 * into the corresponding struct value.
 */

struct settings {
	char conf_path[256];
	char ref_file[256];
	char bam_list[256];
	char bed_file[256];
	char snp_vcf_path[256];
	char snp_vcf_name[256];
	char cosm_vcf_path[256];
	char normal_base_path[256];
	char normal_base_out[256];
	char sex_chrs[256];
	char out_name[256];
	char ploidy_str[256];

	char chromosome[100];

	uint32_t output_pos:1,
		 output_exon:1,
		 output_nmetrics:1,
		 verbose:1,
		 use_bed:1,
		 no_extra_columns:1,
		 dummy:27;

	uint32_t min_mq;
	uint32_t min_bq;

	uint32_t default_mq;
	uint32_t default_bq;

	uint32_t min_b_count;

	double pv_freq;
	double prior_map_error;
};

/* global settings object */
struct settings settings;

// }}}

// Context struct definition {{{

/**
 * The context struct holds the state of the program and is passed
 * around basically everywhere. This is a feeble attempt at avoiding
 * global variables and ridiculous function signatures.
 */
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

	FILE *pos_file;
	FILE *exon_file;
	FILE *nmetrics_file;
};

// }}}

struct extra_data {
	uint32_t mm_count;
	uint32_t total_softclip;
};

// Global blank variant counts dummy value
// (used in dump_vcounts)
static struct variant_counts dummy = { 0 };

// Only call this after all config has been read!
static
void init_dummy(struct variant_table *v)
{
	dummy.total_mq = settings.default_mq;
	dummy.total_bq = settings.default_bq;
	dummy.pos = v->offset;
	dummy.report.pos = v->offset;
	dummy.pop_af = 0.0f;
}

// Recording data from BAM {{{

/* utility */
static inline
int is_mismatch(struct alignment_report rep, struct ref_seq ref_seq_info)
{
	return !( rep.size > 0 &&
			rep.align_type == at_single &&
			rep.data == char_to_b5base(ref_seq_get(ref_seq_info, rep.pos)) );
}

/**
 * Horrible mess of a function that handles overlapping reads
 *
 * The `mq` and `bq` pointers will have their contents MODIFIED to
 * reflect the adjustment in overlap qualities.
 *
 * Note that this function isn't complete- nobody knows what to do to
 * overlap qualities if there is an indel within the overlapping region.
 * See comment here:
 *    https://github.com/samtools/htslib/blob/c4a7c7284d8c144dbb0101d04fd97be008c3b427/sam.c#L2362
 */
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
void update_ab(struct variant_table *v, struct variant_counts *new_entry)
{
	uint32_t total_count, a_total, b_total;
	total_count = new_entry->count_r + new_entry->count_f;

	if (!v->a) {
		v->a = new_entry;
		return;
	}
	if (!v->b) {
		if (new_entry == v->a) return;

		v->b = new_entry;
		return;
	}

	a_total = v->a->count_r + v->a->count_f;
	b_total = v->b->count_r + v->b->count_f;

	if (total_count > a_total) {
		v->b = v->a;
		v->a = new_entry;
	} else if (total_count > b_total && new_entry != v->a) {
		v->b = new_entry;
	}
}

static
uint32_t get_insert_size(bam1_t *bam, bam1_t *mbam)
{
	assert(bam != NULL);

	uint32_t s1, s2, e1, e2;
	s1 = bam->core.pos;
	e1 = bam_endpos2(bam);

	if (mbam == NULL) {
		return e1 - s1;
	}

	s2 = mbam->core.pos;
	e2 = bam_endpos2(mbam);

	return MAX(e1, e2) - MIN(s1, s2);
}

static
uint32_t get_read_orientation(bam1_t *bam)
{
	return ! !(bam->core.flag & 0x10);
}

/**
 * This is the meat of the program. The rest is just fiddling around
 * with data to prepare it or fiddling with what this produces to
 * output it. Pretty much all the relevant calculations occur here.
 *
 * @param context: the current state of the pileup engine-- the current
 * read, offset in the read, offset in the genome, what not.
 *
 * @param rep: the current alignment report. See report.h for more
 * details.
 *
 * @param extra_data_p: parameter holds aggregate information about
 * the current read. It is generated by a function that is called
 * similarly to this one, just before. See `calc_alignments` in
 * cigar.c for more details.
 *
 * The `variant_table` struct holds information about the current
 * position. There is generally one member for each kind of column in
 * the POS file output. Note that there isn't really a 1-1 mapping
 * here but if you examine the definition in `variant_table.h`, it
 * should be clear what everything does. An important note is that
 * averages are stored as totals and are divided by the read depth
 * later.
 *
 * The `variant_counts` struct holds information about how many times
 * each alignment report is seen at each position. This is how the A
 * and B counts are calculated.
 *
 * If you want to add a column, there are some considerations:
 *
 *     Is it PER-POSITION? (example: REF, ALT, etc) Then add a field
 *     to `struct variant_table` and update it here based on `rep`.
 *
 *     Or is it PER-A/B? (example: AmeanSC, BmeanRO, etc) Then add a
 *     field to `struct variant_counts` and update here based on
 *     `rep`.
 *
 * I find that the best way to learn this is to see an example of
 * it. The mean read orientation field is a good example of a field I
 * added after-the-fact. Search for `total_read_orientation` here and
 * in variant_table.h and duplicate what you see for any new PER-A/B
 * fields you might want to add. You would need to modify the
 * `variant_counts` struct as well as the `record_match` and
 * `dump_variant_info` functions. That should be enough, aside from
 * adding any helper functions (such as `get_read_orientation`)
 *
 * For a more advanced PER-A/B field, search `total_pmm` and
 * `mm_count`.
 *
 * For a PER-POSITION field, look in `dump_vcounts` and search for
 * things you see printed out there, for example "read_count")
 *
 * Removing columns is easy, just modify the `fprintf` calls in
 * `dump_vcounts` or `dump_variant_info`.
 */
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
	uint32_t total_softclip = extra_data->total_softclip;

	// REF1: This is the comment to which record_aggregate()
	// refers. Note that the mm_count is decremented by the size
	// of the current report. This is the current mismatch is not
	// counted in PMM.
	if (is_mismatch(rep, context->ref_seq_info)) {
		mm_count -= rep.size;
	}

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

	/* Now `vtentry` represents the record for the position
	 * `rep.pos` in the POS file. It is manipulated to record the new
	 * information.
	 *
	 * If the information is global to the POS file line (for
	 * example, `read_count`), then it is manipulated directly.
	 *
	 * If the information is local to a possible variant (such as
	 * mean mapping quality), it must go in the record for the
	 * variant. The piece of data that distinguishes variants is
	 * the report, so that is the key of the `variant_counts` hash
	 * table. The search for the key occurs close to this comment
	 * below, with the `HASH_FIND` macro call. If the call was
	 * successful, then the variable `vcounts` is a pointer to the
	 * table that maps the variant to its details.
	 */

	vtentry->read_count++;

	if (mq >= settings.min_mq && bq >= settings.min_bq) {

		vtentry->read_count_pass++;

		/* check if this is a mismatch */
		if (is_mismatch(rep, context->ref_seq_info)) {
			vtentry->mismatch_count++;
		}

		/* Update the entry */
		struct variant_counts *vcounts_table = vtentry->counts;

		/* We use a zero'd out key holder to avoid any struct padding nonsense */
		struct variant_counts key_holder = {0};
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

		vcounts->total_softclip += total_softclip;
		vcounts->total_insert_size += get_insert_size(bam, mbam);
		vcounts->total_read_orientation += get_read_orientation(bam);

		assert(mm_count <= (uint32_t) bam->core.l_qseq);
		vcounts->total_pmm += (double) mm_count / bam->core.l_qseq;

		update_ab(vtentry, vcounts);

		/* The pointer may have changed */
		vtentry->counts = vcounts_table;
	}

	/* Update the vtable */
	context->bmi->itr_list[sample_index].vtable = vtable;

}

static
void record_aggregate(	struct context *context,
			struct alignment_report rep,
			void *extra_data_p)
{

	struct extra_data *extra_data = extra_data_p;
	struct ref_seq ref_seq_info = context->ref_seq_info;

	// Count total mismatch
	//
	// NOTE: the mm_count calculated here is not the same that is
	// used to calculate PMM. See record_match() for more
	// details. (find REF1)
	int is_mm = is_mismatch(rep, ref_seq_info);
	if (is_mm) {
		extra_data->mm_count += rep.size;
	}

	if (rep.align_type == at_sclip) {
		extra_data->total_softclip += rep.size;
	}
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

static __attribute__((unused))
void sort_vcounts(struct variant_table *vtable)
{
	HASH_SORT(vtable->counts, cmp_vcounts);
}

static
void get_bcf_entries(	struct context *context,
			uint32_t offset,
			bcf1_t **pop_entry)
{
	bcf1_t *line;
	int pop_done;

	uint32_t pos;

	pop_done = 0;

	// TODO: figure out return values for this function and handle errors
	// properly
	bcf_sr_seek(context->bcf_reader, context->target_name, offset - 1);

	while (!pop_done) {
		if (!bcf_sr_next_line(context->bcf_reader)) {
			if (context->bcf_reader->errnum) {
				err_printf("vcf error: %s\n", bcf_sr_strerror(context->bcf_reader->errnum));
			}
			break;
		}

		if (!pop_done && bcf_sr_has_line(context->bcf_reader, SNP_VCF_INDEX)) {
			line = bcf_sr_get_line(context->bcf_reader, SNP_VCF_INDEX);
			pos = (uint32_t) line->pos;
			if (pos+1 == offset) {
				*pop_entry = line;
			}
			if (pos+1 >= offset) {
				pop_done = 1;
			}
		}
	}

}

static
float get_af_true(	struct context *context,
			int reader_index,
			bcf1_t *entry,
			struct alignment_report rep)
{
	int i, num_alleles, pv_freq_coeff;
	uint32_t ref, alt;

	float pv_freq, result;
	bcf_info_t *info;
	bcf_hdr_t *header;

	pv_freq = settings.pv_freq;

	if (entry == NULL) goto no_af;

	header = context->bcf_reader->readers[reader_index].header;

	num_alleles = entry->n_allele;

	info = bcf_get_info(header, entry, "AF");

	// this macro takes an allele and returns something
	// that can be compared with the report data
#define NORMALIZE_ALLELE(rep, al) \
		( (rep).align_type == at_single ? \
			char_to_b5base( (al)[0] ) : \
			str_to_b5seq( (al) ) )

	ref = NORMALIZE_ALLELE(rep, entry->d.allele[0]);

	if (info == NULL) goto no_af;

	if (rep.data == ref) {
		// if #alts <= 3, then calculate 1 - sum(AF) - (3 - #alts) * pv_freq
		// otherwise,  calculate 1 - sum(AF)
		result = 1;
		pv_freq_coeff = 3;

		for (i = 1; i < num_alleles; i++) {
			result -= ((float *) info->vptr)[i-1];
			pv_freq_coeff -= 1;
		}

		if (pv_freq_coeff < 0) pv_freq_coeff = 0;

		result -= pv_freq_coeff * pv_freq;
		return result;

	} else {
		// Loop through the alleles in entry and report the correct AF
		// If there is no AF, use pv_freq

		for (i = 1; i < num_alleles; i++) {
			alt = NORMALIZE_ALLELE(rep, entry->d.allele[i]);
			if (rep.data == alt) {
				return ((float *) info->vptr)[i-1];
			}
		}

		return pv_freq;
	}

	// Shouldn't be reached
	assert(0);
	return 0;

no_af:
	return is_mismatch(rep, context->ref_seq_info) ? pv_freq : 1 - 3*pv_freq;
}
#undef NORMALIZE_ALLELE


void find_ab(struct variant_table **vtables, uint32_t n_samples,
	     /* out */ struct variant_counts **a, struct variant_counts **b)
{
	uint32_t sample_idx;
	uint32_t total, max, second_max;

	struct variant_table *vt;
	struct variant_counts *counts, *tmp;

	max = 0;
	second_max = 0;

	*a = NULL;
	*b = NULL;

	for (sample_idx = 0; sample_idx < n_samples; sample_idx++) {
		vt = vtables[sample_idx];
		if (vt == NULL) continue;

		HASH_ITER(hh, vt->counts, counts, tmp) {
			total = counts->count_f + counts->count_r;
			if (total > max) {
				second_max = max;
				*b = *a;
				max = total;
				*a = counts;
			} else if (total > second_max) {
				// Check if it's the same allele, just in a different sample.
				// This case is ignored.
				if (a != NULL && counts->report.data == (*a)->report.data) {
					continue;
				}
				second_max = total;
				*b = counts;
			}
		}
	}

}

void get_ab(struct variant_table *v, struct variant_counts **a, struct variant_counts **b)
{
	*a = v && v->a ? v->a : &dummy;
	*b = v && v->b ? v->b : &dummy;
}

__attribute__((unused))
int check_ab(struct variant_counts *a, struct variant_counts *b)
{
	(void) a; (void) b; // to avoid gcc warning

	/* If this fails, then a and b are both being set
	 * but to the same thing which is bad */
	assert((a == &dummy && b == &dummy) || a != b);

	/* If this fails then a and b are not ordered correctly
	 * and are most likely completely wrong. */
	assert(a->count_f + a->count_r >= b->count_f + b->count_r);

	return 1;
}

static
void populate_afs(struct context *context, bcf1_t *pop_entry, struct variant_table *v)
{
	double a_pop_af, b_pop_af;
	struct variant_counts *a, *b;

	get_ab(v, &a, &b);
	assert(check_ab(a, b));

	a_pop_af = get_af_true(context, SNP_VCF_INDEX, pop_entry, a->report);
	b_pop_af = get_af_true(context, SNP_VCF_INDEX, pop_entry, b->report);

	a->pop_af = a_pop_af;
	b->pop_af = b_pop_af;
}

static
uint32_t dump_variant_info(	struct context *context,
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
	const char *format_string;
	if (settings.no_extra_columns) {
		format_string = VARIANT_TYPE_FORMAT "%d\t%g\t%g\t%g\t%g\t%g\t%g";
	} else {
		format_string = VARIANT_TYPE_FORMAT "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g";
	}

	fprintf(f, format_string,
		VARIANT_TYPE,
		read,
		(double) vc->count_f / 2,
		(double) vc->count_r / 2,
		(double) vc->total_bq / total,
		(double) vc->total_mq / total,
		vc->total_pmm / total,
		(double) vc->total_read_pos / total,
		(double) vc->total_softclip / total,
		(double) vc->total_insert_size / total,
	        (double) vc->total_read_orientation / total);
#undef VARIANT_TYPE_FORMAT
#undef VARIANT_TYPE

	return read;
}

static
void dump_nm_data(	struct context *context,
			uint32_t offset )
{
	FILE *f = context->pos_file;

	struct bam_multi_itr *bmi = context->bmi;

	struct nm_tbl *nmt = bmi->itr_list[context->sample_index].nmt;
	struct nm_entry *ent = nm_tbl_get(nmt, offset);

	if (ent == NULL) {
		err_printf("no normal metrics data for chromosome %s, offset %d\n",
				settings.chromosome, offset);
		fprintf(f, "nan\tnan\tnan\tnan");
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
			/* Global A/B alleles for this offset */
			struct alignment_report global_a_report,
			struct alignment_report global_b_report,
			int max_delete_size	)
{
	FILE *f = context->pos_file;

	// Used for zero-d out A/B in case of mismatch with global A/B
	struct variant_counts dummy_a = {0}, dummy_b = {0};
	dummy_a.report = global_a_report;
	dummy_b.report = global_b_report;

	struct variant_counts *local_a = NULL, *local_b = NULL;
	struct ref_seq ref_seq_info = context->ref_seq_info;

	float a_pop_af, b_pop_af;

	uint32_t ref_allele, ref_allele_partial;
	uint32_t a_alt, b_alt;

	int32_t cosm_count_a, cosm_count_b, cosm_count_real;

	// Same zero'd out key holder technique here
	struct variant_counts key_holder_a = {0};
	struct variant_counts key_holder_b = {0};

	// Search for a and b in `v`

	// Prepare the key_holder for searching for a
	key_holder_a.report.align_type = global_a_report.align_type;
	key_holder_a.report.pos = global_a_report.pos;
	key_holder_a.report.data = global_a_report.data;
	key_holder_a.report.size = global_a_report.size;

	// Perform search
	HASH_FIND(hh, v->counts, &key_holder_a.report, sizeof(global_a_report), local_a);

	// Handle not-found case
	if (!local_a) {
		local_a = &dummy_a;
	}


	// Repeat for b
	key_holder_b.report.align_type = global_b_report.align_type;
	key_holder_b.report.pos = global_b_report.pos;
	key_holder_b.report.data = global_b_report.data;
	key_holder_b.report.size = global_b_report.size;

	// Perform search
	HASH_FIND(hh, v->counts, &key_holder_b.report, sizeof(global_b_report), local_b);

	// Handle not-found case
	if (!local_b) {
		local_b = &dummy_b;
	}

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

	a_alt = dump_variant_info(context, f, local_a, ref_allele_partial);
	fprintf(f, "\t");
	b_alt = dump_variant_info(context, f, local_b, ref_allele_partial);
	fprintf(f, "\t");

	a_pop_af = local_a->pop_af;
	b_pop_af = local_b->pop_af;

	cosm_count_a = get_cosmic_count(v->offset, a_alt);
	cosm_count_b = get_cosmic_count(v->offset, b_alt);
	cosm_count_real = MAX(0, MAX(cosm_count_a, cosm_count_b));

	fprintf(f, "%g\t%g\t%d\t", a_pop_af, b_pop_af, cosm_count_real);

	dump_nm_data(context, v->offset);
	fprintf(f, "\n");
}

static
void dump_blank_vcounts(struct context *context,
			uint32_t offset,
			struct alignment_report global_a_report,
			struct alignment_report global_b_report,
			int max_delete_size)
{

	struct variant_table vt = {0};

	vt.sample_index = context->sample_index;
	vt.offset = offset;
	vt.tid = context->tid;

	dump_vcounts(context, &vt,
		     global_a_report, global_b_report,
		     max_delete_size);
}

static
void dump_nmetrics(struct context *context, uint32_t offset, struct nm_entry avgs, int count)
{
	(void)count;

	fprintf(context->nmetrics_file, "%d\t%d\t%g\t%g\t%g\t%g\n",
			context->tid + 1,
			offset,
			avgs.norm_read_depth,
			avgs.prob_map_err,
			avgs.read_pass,
			avgs.a_or_b);
}

static
void dump_blank_nmetrics(struct context *context, uint32_t offset)
{
	fprintf(context->nmetrics_file, "%d\t%d\t0\tNaN\tNaN\tNaN\n",
			context->tid + 1,
			offset);
}

static inline
int get_delete_size(struct variant_counts *vc)
{
	if (vc == NULL) return -1;

	return vc->report.align_type == at_del ? vc->report.size : 0;
}

static
int get_max_delete_size(struct variant_counts *a, struct variant_counts *b)
{
	// The purpose of this function is not immediately obvious.
	// It calculates, for a given offset, the size of the largest
	// event.

	// How have I defined the size of an event? If it's a single
	// base change or an insertion, it's 1. If it's a deletion,
	// it's the number of bases deleted.

	// Why is this necessary? Because when calculating the REF
	// field for pos file output for a position where a deletion
	// is the most common event, the size of the deletion is important.

	// Let's say we have a position where allele A is a deletion
	// of size 3 and allele B is an single point change. Then, the
	// REF and ALT fields might look like this:

	// REF: ACCC    ALT_A: A    ALT_B: TCCC

	// It's clear that to calculate the REF field, the size of the
	// A deletion must be known. That's what this function
	// calculates.


	int max_delete_size = 0, max_of_ab;

	max_of_ab = MAX(get_delete_size(a), get_delete_size(b));
	if (max_of_ab > max_delete_size) max_delete_size = max_of_ab;

	return max_delete_size;
}

static
int has_enough(struct context *context, uint32_t offset)
{
	// The purpose of this function is to determine if, at
	// `offset`, there are enough B's to include this offset in
	// the output. Note: we don't use the global B's.

	uint32_t sample_idx;
	int enough = 0;
	struct variant_counts *counts, *tmp;
	struct variant_table *vtable;

	for (sample_idx = 0; sample_idx < context->bmi->num_iters; sample_idx++) {
		HASH_FIND_INT(context->bmi->itr_list[sample_idx].vtable, &offset, vtable);
		if (vtable == NULL) continue;

		HASH_ITER(hh, vtable->counts, counts, tmp) {
			if (is_mismatch(counts->report, context->ref_seq_info) &&
					counts->count_f + counts->count_r >=
						2 * settings.min_b_count) {
				enough = 1;
			}
		}
	}

	return enough;
}

static
void flush_results(struct context *context, uint32_t begin, uint32_t end)
{
	int max_delete_size, enough;
	struct variant_table *vt, *v;
	struct variant_counts *a, *b;

	struct alignment_report global_a_report = {0};
	struct alignment_report global_b_report = {0};

	uint32_t offset, sample_idx;
	uint32_t n_samples = context->bmi->num_iters;

	struct nm_entry ent = {0};
	struct nm_tbl *tbl;

	int ploidy; // for normal metrics calculation

	// Array to hold variant_tables for each sample_index at a
	// certain position
	struct variant_table **vtables;
	vtables = malloc(n_samples * sizeof(*vtables));

	bcf1_t *pop_entry;

	// clamp in region
	if (begin < context->reg_start) begin = context->reg_start;
	if (end > context->reg_end + 1) end = context->reg_end + 1;

	for (offset = begin; offset < end; offset++) {
		enough = has_enough(context, offset);

		// If normal metrics output isn't a concern and there
		// is nothing of interest at this position, SKIP IT
		//
		// For future reference, this check is crucial to the
		// performance of this program, especially when only
		// doing pos/exon output.
		if (!enough && !settings.output_nmetrics) continue;

		pop_entry = NULL;
		get_bcf_entries(context, offset, &pop_entry);

		// This table is used to average normal metrics over
		// the samples
		tbl = nm_tbl_create();

		// First, populate vtables (array of variant_table for each
		// sample at `offset`)
		for (sample_idx = 0; sample_idx < n_samples; sample_idx++) {
			vt = context->bmi->itr_list[sample_idx].vtable;
			HASH_FIND_INT(vt, &offset, vtables[sample_idx]);
		}

		// Now, determine A/B
		find_ab(vtables, n_samples, &a, &b);
		if (a != NULL) {
			global_a_report = a->report;
		}

		if (b != NULL) {
			global_b_report = b->report;
		}

		max_delete_size = get_max_delete_size(a, b);

		for (sample_idx = 0; sample_idx < n_samples; sample_idx++) {
			context->sample_index = sample_idx;
			v = vtables[sample_idx];

			if (v != NULL) {
				init_dummy(v);
				populate_afs(context, pop_entry, v);
			}

			if (max_delete_size >= 0) {
				if (settings.output_pos) {
					if (v == NULL) {
						dump_blank_vcounts(context, offset,
								   global_a_report,
								   global_b_report,
								   max_delete_size);
					} else {
						dump_vcounts(context, v,
							     global_a_report,
							     global_b_report,
							     max_delete_size);
					}
				}
			}

			if (settings.output_nmetrics) {
				struct variant_table dummy = {0}, *rv;
				rv = v == NULL ? &dummy : v;

				// This is guaranteed to not segfault at least
				ploidy = settings.ploidy_str[sample_idx] - '0';

				nmcalc(context, rv, settings.prior_map_error, ploidy, &ent);
				nm_tbl_add(tbl, ent, 1 /* only compute averages */);

				if (v == &dummy) {
					v = NULL;
				}
			}
		}

		if (settings.output_nmetrics) {
			if (tbl->count == 0) {
				dump_blank_nmetrics(context, offset);
			} else {
				dump_nmetrics(context, offset, tbl->avgs, tbl->count);
			}
		}

		nm_tbl_destroy(tbl);
	}

	free(vtables);

}

// }}}

// Checking various properties of alignments {{{

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
	uint32_t min_mq, min_bq, default_mq, default_bq, min_b_count;
	double pv_freq, prior_map_error;
	char *tmp;

	// See config_keys.h for the real values
	GVM_CONFIG_CHECK_STROPT(bam_list);
	GVM_CONFIG_CHECK_STROPT(ref_file);
	GVM_CONFIG_CHECK_STROPT(bed_file);
	GVM_CONFIG_CHECK_STROPT(snp_vcf_path);
	GVM_CONFIG_CHECK_STROPT(snp_vcf_name);
	GVM_CONFIG_CHECK_STROPT(cosm_vcf_path);
	GVM_CONFIG_CHECK_STROPT(normal_base_path);
	GVM_CONFIG_CHECK_STROPT(normal_base_out);
	GVM_CONFIG_CHECK_STROPT(out_name);

	GVM_CONFIG_CHECK_NUMOPT(min_bq);
	GVM_CONFIG_CHECK_NUMOPT(min_mq);
	GVM_CONFIG_CHECK_NUMOPT(default_bq);
	GVM_CONFIG_CHECK_NUMOPT(default_mq);
	GVM_CONFIG_CHECK_NUMOPT(min_b_count);
	GVM_CONFIG_CHECK_NUMOPT(pv_freq);

	if (settings.output_nmetrics) {
		GVM_CONFIG_CHECK_NUMOPT(prior_map_error);
	}

	return 2;

}

int config_read(const char *cfg_fn, struct settings *s)
{
	FILE *f = fopen(cfg_fn, "r");
	int result = 1;
	int set_result;

	/* FSM */
	int state = 0; /* 0 = need key, 1 = need value */

	char value_buf[1024];
	char *option_value;

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
			err_printf("error in %s at line %lu, col %lu: %s\n",
			           cfg_fn,
			           parser.problem_mark.line + 1,
			           parser.problem_mark.column + 1,
			           parser.problem);

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
			err_printf("Warning: no alias support (anchor %s)\n", event.data.alias.anchor);
			break;
		case YAML_SCALAR_EVENT:
			option_value = (char *) event.data.scalar.value;
			if (state == 0) { /* This means we need a key */
				strncpy(value_buf, option_value, sizeof(value_buf));
				state = 1;
			} else { /* We have a key, now we need a value */

				// 0 indicates failure
				// 1 indicates success and option set
				// 2 indicates success and option ignored
				set_result = config_try_set_option(value_buf, option_value, s);

				if (set_result == 0) {
					result = 0;
					goto done1;
				}

				if (set_result == 1) {
					verbose_fprintf(stderr, "  config: %s = %s\n", value_buf, option_value);
				}

				state = 0;
			}

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
int write_exon_line(struct context *context, uint32_t start, uint32_t end)
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


	for (i = 0; i < bmi->num_iters; i++) {
		avgs = bmi->itr_list[i].nmt->avgs;
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

	return 1;
}

static
int open_out_files(struct context *context)
{
	char *pos_fn, *exon_fn, *nmetrics_fn;
	FILE *pos_file, *exon_file, *nmetrics_file;
	int cannot_continue = 0;

	if (settings.output_pos) {
		pos_fn = malloc(strlen(settings.out_name) +
				strlen(settings.chromosome) +
				strlen("_pos.txt") + 2);

		if (pos_fn == NULL) {
			err_printf("unable to allocate memory.");
			return 1;
		}
		sprintf(pos_fn, "%s_%s_pos.txt", settings.out_name, settings.chromosome);

#ifdef DEBUG
		pos_file = stdout;
#else
		pos_file = fopen(pos_fn, "w");
#endif
		if (pos_file == NULL) {
			err_printf("unable to open %s:", pos_fn);
			perror("");
			cannot_continue = 1;
		}

		free(pos_fn);
	} else {
		pos_file = NULL;
	}


	if (settings.output_exon) {
		exon_fn = malloc(strlen(settings.out_name) +
				 strlen(settings.chromosome) +
				 strlen("_exon.txt") + 2);

		if (exon_fn == NULL) {
			err_printf("unable to allocate memory.");
			return 1;
		}

		sprintf(exon_fn, "%s_%s_exon.txt", settings.out_name, settings.chromosome);

		exon_file = fopen(exon_fn, "w");
		if (exon_file == NULL) {
			err_printf("unable to open %s:", exon_fn);
			perror("");
			cannot_continue = 1;
		}

		free(exon_fn);
	} else {
		exon_file = NULL;
	}

	if (settings.output_nmetrics) {
		if (strlen(settings.ploidy_str) == 0) {
			err_printf("normal metrics requested yet --ploidystr not provided\n");
			return 1;
		}

		if (strlen(settings.normal_base_out) == 0) {
			err_printf("normal metrics requested yet " GVM_CONFIG_normal_base_out " not set in config\n");
			return 1;
		}

		nmetrics_fn = malloc(strlen(settings.normal_base_out) +
				     strlen(settings.chromosome) +
				     strlen(".txt") + 2);

		if (nmetrics_fn == NULL) {
			err_printf("unable to allocate memory.");
			return 1;
		}

		sprintf(nmetrics_fn, "%s%s.txt", settings.normal_base_out, settings.chromosome);
		nmetrics_file = fopen(nmetrics_fn, "w");
		if (nmetrics_file == NULL) {
			err_printf("unable to open %s:", nmetrics_fn);
			perror("");
			cannot_continue = 1;
		}

		free(nmetrics_fn);
	} else {
		nmetrics_file = NULL;
	}

	if (cannot_continue) {
		return 1;
	}

	context->pos_file = pos_file;
	context->exon_file = exon_file;
	context->nmetrics_file = nmetrics_file;

	return 0;
}


#define GVM_CALL_CALC_ALIGN(R) \
	do { \
		memset(&ed, 0, sizeof(ed)); \
		R = calc_alignments( \
				context->bam, \
				context->ref_seq_info, \
				record_match, \
				record_aggregate, \
				context, \
				&ed); \
	} while (0)

#define GVM_CHECK_RESULT(R) \
	if (R < 0) { \
		err_printf("failed to calculate alignments (err code %d)\n", R); \
		return EXIT_FAILURE; \
	}

/**
 * Calculate metrics for positions between `start` and `end`
 *
 * Please mind the macro calls to GVM_CALL_CALC_ALIGN-- this calls the
 * function in cigar.c to process the read by using the cigar string.
 */
static
int do_region(struct context *context, uint32_t start, uint32_t end) // {{{
{
	verbose_fprintf(stderr, "--- region:  %d  to  %d  \n", start, end);

	struct bam_multi_itr *bmi = context->bmi;
	if (!bmi_query(bmi, settings.chromosome, start, end)) {
		return 0;
	}

	int result;

	struct bam_mate_table *left_behind, *tmp;
	int32_t sample_index;

	struct extra_data ed;

	bam1_t *bam, *mbam;

	context->reg_start = start;
	context->reg_end = end;

	//nm_query(context->nmi, start - 1, end);
	//nm_tbl_slurp(context->nmt, context->nmi);

	// Seek the VCF reader to the beginning of the region
	bcf_sr_seek(context->bcf_reader, context->target_name, start);

	if (settings.output_pos || settings.output_nmetrics) {
		// Main processing loop {{{
		while ( (sample_index = bmi_next(bmi, &bam)) >= 0 ) {
			if (!check_flags(bam)) continue;

			mbam = NULL;

			if (!is_split_read(bam)) {
				bmi->itr_list[sample_index].bmt =
					bmt_register(bmi->itr_list[sample_index].bmt, /* The pair table */ bam, &mbam, /* The pair */ &result /* Out: result value */);
				
				if (result == 1) continue;
			}

			context->sample_index = (uint32_t) sample_index;
			context->bam = mbam;
			context->mbam = bam;


			if (mbam != NULL) {
				GVM_CALL_CALC_ALIGN(result);
				GVM_CHECK_RESULT(result);
			}

			context->bam = bam;
			context->mbam = mbam;

			GVM_CALL_CALC_ALIGN(result);
			GVM_CHECK_RESULT(result);

			if (mbam != NULL) {
				bam_destroy1(mbam);
			}

		}
		// }}}

		context->mbam = NULL;
		for (sample_index = 0; sample_index < bmi->num_iters; sample_index++) {
			HASH_ITER(hh, bmi->itr_list[sample_index].bmt, left_behind, tmp) {
				context->bam = left_behind->bam;
				context->sample_index = sample_index;
				GVM_CALL_CALC_ALIGN(result);
				GVM_CHECK_RESULT(result);
			}
		}

		flush_results(context, start, end + 1);
	}

	// It should be empty now!

	if (settings.output_exon) {
		write_exon_line(context, start, end);
	}

	bmi_cleanup(bmi);

	//fflush(context->pos_file);

	return 1;
} // }}}

#undef GVM_CALL_CALC_ALIGN
#undef GVM_CHECK_RESULT

static
int try_parse_region(char *regstr, uint32_t *begin, uint32_t *end)
{
	char *endptr;
	*begin = strtoul(regstr, &endptr, 0);
	if (endptr == regstr) {
		// Could not parse
		return 0;
	}

	if (*endptr != '-') {
		return 0;
	}

	regstr = endptr+1;
	*end = strtoul(regstr, &endptr, 0);

	if (endptr == regstr) {
		return 0;
	}

	if (*endptr != '\0') {
		return 0;
	}

	return 1;
}

int main(int argc, char **argv) // {{{
{
	struct gengetopt_args_info args_info;
	faidx_t *ref_idx;
	char *ref_seq_data, *snp_vcf_fname;
	int len, result, tid;
	struct ref_seq ref_seq_info;
	struct bam_multi_itr *bmi;
	struct context context;

	// The region specified on the commandline
	uint32_t cmdline_region_start, cmdline_region_end;

	region_handle_func regfn;

	// Command line settings {{{
	if (cmdline_parser(argc, argv, &args_info) != 0) {
		return EXIT_FAILURE;
	}

	settings.use_bed = 1;
	cmdline_region_start = 0;
	cmdline_region_end = 0;
	if (args_info.region_given) {
		if (!try_parse_region(args_info.region_arg, &cmdline_region_start, &cmdline_region_end)) {
			err_printf("Unable to parse region string %s\n", args_info.region_arg);
			return EXIT_FAILURE;
		}

		settings.use_bed = 0;
	}

	settings.no_extra_columns = 1;
	if (args_info.no_extra_columns_given) {
		settings.no_extra_columns = 0;
	}

	// Don't have to check if conf_arg or chr_arg are given because they're required
	strncpy(settings.conf_path, args_info.conf_arg, sizeof(settings.conf_path));
	strncpy(settings.chromosome, args_info.chr_arg, sizeof(settings.chromosome));
	// This option is not required, hence the check
	if (args_info.ploidystr_given) {
		strncpy(settings.ploidy_str, args_info.ploidystr_arg, sizeof(settings.ploidy_str));
	}

	/* output settings */
	settings.output_pos = args_info.output_pos_flag;
	settings.output_exon = args_info.output_exon_flag;
	settings.output_nmetrics = args_info.output_normal_flag;

	settings.verbose = args_info.verbose_given;

	verbose_fprintf(stderr, "%s (bugs to %s)\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
	verbose_fprintf(stderr, "Configuration: %s\n", settings.conf_path);
	verbose_fprintf(stderr, "Chromosome: %s\n", settings.chromosome);
	verbose_fprintf(stderr, "Outputs: ");
	if (settings.output_pos) {
		verbose_fprintf(stderr, "pos ");
	}
	if (settings.output_exon) {
		verbose_fprintf(stderr, "exon ");
	}
	if (settings.output_nmetrics) {
		verbose_fprintf(stderr, "normalmetrics ");
	}
	verbose_fprintf(stderr, "\n");

	if (settings.output_pos) {
		verbose_fprintf(stderr, "extra pos columns: %s\n",
				settings.no_extra_columns ? "no" : "yes");
	}

	if (!config_read(settings.conf_path, &settings)) {
		return EXIT_FAILURE;
	}

	cmdline_parser_free(&args_info);

	// }}}

	// BAM file loading {{{
	bmi = bmi_create(settings.bam_list, settings.normal_base_path, settings.chromosome);
	if (bmi == NULL) {
		return EXIT_FAILURE;
	}

	assert(bmi->num_iters > 0);

	// If normal metrics calculation is requested, check that len(ploidystr) == len(bamlist)
	if (settings.output_nmetrics && strlen(settings.ploidy_str) != bmi->num_iters) {
		err_printf("bamlist has length %lu, yet ploidystr has length %d\n",
		           strlen(settings.ploidy_str),
		           bmi->num_iters);
		return EXIT_FAILURE;
	}

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

	bcf_sr_set_regions(context.bcf_reader, settings.bed_file, 1);

	result = bcf_sr_add_reader(context.bcf_reader, snp_vcf_fname);
	if (result == 0) {
		err_printf("failed to load SNP VCF: %s\n", snp_vcf_fname);
		err_printf("error: %s\n", bcf_sr_strerror(context.bcf_reader->errnum));
		return EXIT_FAILURE;
	}

	if (strlen(settings.cosm_vcf_path) > 0) {
		verbose_fprintf(stderr, "Loading cosmic count data...\n");
		load_cosmic_table(settings.cosm_vcf_path, settings.chromosome);
	} else if (settings.output_pos || settings.output_exon) {
		err_printf("pos or exon file requested, but no cosmic vcf provided\n");
		return EXIT_FAILURE;
	}

	free(snp_vcf_fname);
	// }}}

	if (open_out_files(&context) != 0) {
		return EXIT_FAILURE;
	}

	// Context initialization {{{
	context.ref_seq_info = ref_seq_info;
	context.bmi = bmi;
	context.tid = tid;
	context.target_name = settings.chromosome;

	if (context.tid == -1) {
		err_printf("warning: no bam files loaded\n");
	}
	// }}}

	regfn = (region_handle_func) do_region;

	if (settings.use_bed) {
		bedf_forall_region_chr(
				&context,
				settings.bed_file,
				settings.chromosome,
				regfn);
	} else {
		regfn(&context, cmdline_region_start, cmdline_region_end);
	}

	// Cleanup {{{
#ifdef CLEANUP /* There's really no reason to free this right before termination */
	bmi_destroy(bmi);
	free(ref_seq_data);
	fai_destroy(ref_idx);
	bcf_sr_destroy(context.bcf_reader);

	if (context->pos_file != NULL) {
		fclose(context->pos_file);
	}

	if (context->exon_file != NULL) {
		fclose(context->exon_file);
	}
#endif

	// }}}

	return EXIT_SUCCESS;
} // }}}
