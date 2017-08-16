#include <string.h>
#include <stdio.h>

#include "cigar.h"
#include "base_seq_repr.h"
#include "ref_seq.h"
#include "utils.h"

#include "htslib/sam.h"


// returns negative on error, otherwise if not
int run_cigar(  bam1_t *bam,
		struct ref_seq ref_seq_info,

		uint32_t cigar_op, uint32_t cigar_nextop,
		uint32_t *base_offset,
		uint32_t *seq_index,
		reporter_func rfunc,
		void *context,
		void *extra_data)
{

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
			report.data = char_to_b5base(ref_seq_get(ref_seq_info, report.pos));
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

int calc_alignments(	bam1_t *bam,
			struct ref_seq ref_seq_info,
			reporter_func rfunc,
			reporter_func rfunc_alt,
			void *context,
			void *extra_data)
{
	uint32_t offset = bam->core.pos + 1;

	uint32_t cigar_len = bam->core.n_cigar;
	if (cigar_len == 0) return 0;

	uint32_t *cigar_data = (uint32_t *) bam_get_cigar(bam);
	uint32_t *cigar_ptr;
	uint32_t cigar_nextop;

	uint32_t seq_index = 0;

	int result = 0;

	uint32_t i;

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

		offset_backup = offset;
		seq_index_backup = seq_index;

		result = run_cigar(
			bam,
			ref_seq_info,
			*cigar_ptr, cigar_nextop,
			&offset,          /* modified, needs backup */
			&seq_index,       /* modified, needs backup */
			rfunc_alt,        /* Alternate report function */
			context,
			extra_data       /* modified, no backup */
		);

		/* Restore the backups */
		seq_index = seq_index_backup;
		offset = offset_backup;

		result = run_cigar(
			bam,
			ref_seq_info,
			*cigar_ptr, cigar_nextop,
			&offset, /* base */
			&seq_index,
			rfunc,
			context,
			&extra_data
		);



		if (result == -1) {
			err_printf("An error was encountered while running the cigar string\n");
			break;
		}
	}

	return result;
}
