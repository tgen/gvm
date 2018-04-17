// variant table data structures
#ifndef GVM_VARIANT_TABLE_H
#define GVM_VARIANT_TABLE_H
#include <stdint.h>

#include "uthash.h"

#include "report.h"

/*
 * 
 * we have a hash table
 * that maps ints (genome offset)
 * to ANOTHER HASH of [struct] variant_counts
 * 
 * for every report, we find out the corresponding reference genome value
 * then we increment its stuff
 * 
 */

struct variant_counts {
	struct alignment_report report; /* hash key */

	uint32_t pos; /* The value of report.pos will alawys be zero so we need this */

	/* These are stored as integers and each read is DOUBLE
	 * COUNTED! This is because there can be half reads (in the
	 * case of overlap, both the f and the r strand counts get 0.5
	 *
	 * Therefore, in order to calculate the actual read depth, you
	 * would do (count_f + count_r) / 2.0 */

	uint32_t count_f;
	uint32_t count_r;

	/* I keep track of the totals, not the means, because it's
	 * easier to keep the total and calculate the mean at the end
	 * instead of constantly updating the mean */

	uint32_t total_bq;
	uint32_t total_mq;
	double total_pmm;

	uint32_t total_softclip;
	uint32_t total_insert_size;
	uint32_t total_read_orientation;

	uint32_t total_read_pos;

	double pop_af;
	uint32_t cosmic_count; // FIXME: unused but maybe useful later?

	UT_hash_handle hh;
};

struct variant_table {
	uint32_t offset; /* hash key */
	uint32_t sample_index;
	int32_t tid;
	uint32_t read_count;
	uint32_t read_count_pass;

	uint32_t mismatch_count;

	struct variant_counts *counts;

	// keep track of top two alleles
	struct variant_counts *a;
	struct variant_counts *b;

	UT_hash_handle hh;
};


void destroy_variant_counts(struct variant_counts *vc);
void destroy_variant_table(struct variant_table *vt);

#endif
