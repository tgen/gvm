#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "uthash.h"
#include "htslib/sam.h"

#include "bam_mate_table.h"

struct bam_mate_table *bmt_register(	struct bam_mate_table *bmt,
					bam1_t *bam,
					bam1_t **out,
					int *cont )
{
	char *qname = bam_get_qname(bam);

	struct bam_mate_table *existing;
	HASH_FIND_STR(bmt, qname, existing);

	*cont = 1;
	if (existing) {
		HASH_DEL(bmt, existing);

		*out = existing->bam;

		free(existing->qname);
		free(existing);

		*cont = 0;
	} else {
		if (bam->core.mpos + MAX_READ_SIZE < bam->core.pos) {
			*cont = 0;
			return bmt;
		}

		if (bam->core.pos + MAX_READ_SIZE < bam->core.mpos) {
			*cont = 0;
			return bmt;
		}

		if (bam_is_rev(bam)) {
			// Either:
			//  * we've seen the forward strand and didn't want to pair
			//    because it didn't overlap
			//  * this came first
			if (bam->core.mpos >= bam_endpos(bam)) {
				return bmt;
			}
		}
		existing = malloc(sizeof(struct bam_mate_table));
		existing->qname = malloc(bam->core.l_qname + 1);
		existing->bam = bam_dup1(bam);
		strncpy(existing->qname, qname, bam->core.l_qname);
		HASH_ADD_STR(bmt, qname, existing);
	}

	return bmt;

}

void bmt_destroy(struct bam_mate_table *bmt)
{
	struct bam_mate_table *el, *tmp;

	HASH_ITER(hh, bmt, el, tmp) {
		HASH_DEL(bmt, el);
		free(el->qname);
		bam_destroy1(el->bam);
		free(el);
	}
}

bam1_t *bmt_get_min(struct bam_mate_table *bmt)
{
	uint32_t min_offset = UINT_MAX;
	bam1_t *min_bam = NULL;
	struct bam_mate_table *el, *tmp;

	HASH_ITER(hh, bmt, el, tmp) {
		if ((uint32_t) el->bam->core.pos < min_offset) {
			min_offset = el->bam->core.pos;
			min_bam = el->bam;
		}
	}

	return min_bam;
}

uint32_t bmt_get_min_offset(struct bam_mate_table *bmt)
{
	bam1_t *min_bam = bmt_get_min(bmt);
	if (min_bam == NULL) return -1;
	return min_bam->core.pos;
}

void bmt_print(struct bam_mate_table *bmt)
{
	struct bam_mate_table *el, *tmp;

	HASH_ITER(hh, bmt, el, tmp) {
		printf("(%c,%d,%d) ", bam_is_rev(el->bam)?'R':'F', el->bam->core.pos, el->bam->core.mpos);
	}

	printf("\n");
}
