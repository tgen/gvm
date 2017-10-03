#ifndef GVM_CIGAR_H
#define GVM_CIGAR_H

#include <stdint.h>

#include "htslib/sam.h"
#include "ref_seq.h"
#include "report.h"

int run_cigar(  bam1_t *bam,
		struct ref_seq ref_seq_info,

		uint32_t cigar_op, uint32_t cigar_nextop,
		uint32_t *base_offset,
		uint32_t *seq_index,
		reporter_func rfunc,
		void *context,
		void *extra_data);

int calc_alignments(	bam1_t *bam,
			struct ref_seq ref_seq_info,
			reporter_func rfunc,
			reporter_func rfunc_alt,
			void *context,
			void *extra_data);
#endif
