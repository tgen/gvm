#ifndef GVM_BAM_MATE_TABLE_H
#define GVM_BAM_MATE_TABLE_H

// just a guess, pray that the actual value is lower than this
#define MAX_READ_SIZE 200

#include "uthash.h"

#include "htslib/sam.h"

struct bam_mate_table {
	char *qname; /* key */
	bam1_t *bam; /* make sure to use bam_dup for this */

	UT_hash_handle hh;
};


/*
 * Dumps a bam into the table. If its mate is already there, remove that and
 * return it.  
 *
 * @param bmt the table to use
 * @param bam the alignment to put in the table
 * @param out pointer to `bam`'s mate, or NULL if no mate
 *
 * @return the modified bam_mate_table
 *
 * If *out is set to NULL, memory is allocated. Otherwise, it is set to a value
 * that was removed from the hashtable and MUST be freed with bam1_destroy
 */

struct bam_mate_table *bmt_register(	struct bam_mate_table *bmt,
					bam1_t *bam,
					bam1_t **out,
					int *cont );

void bmt_destroy(struct bam_mate_table *bmt);

/**
 * @param bmt a mate table
 * @return the lowest offset in the table
 */
bam1_t *bmt_get_min(struct bam_mate_table *bmt);

uint32_t bmt_get_min_offset(struct bam_mate_table *bmt);
void bmt_print(struct bam_mate_table *bmt);

#endif
