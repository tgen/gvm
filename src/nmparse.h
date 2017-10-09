#ifndef GVM_NMPARSE_H
#define GVM_NMPARSE_H

#include "htslib/hts.h"
#include "htslib/tbx.h"

#include "uthash.h"

#include <stdint.h>

struct nm_itr {
	htsFile *htsfp;	
	tbx_t *t;
	hts_itr_t *itr;

	int tid;
};

struct nm_entry {	
	// char *chr;
	uint32_t pos; /* key */
	double norm_read_depth;
	double prob_map_err;
	double read_pass;
	double a_or_b;

	UT_hash_handle hh;
};

struct nm_tbl {

	struct nm_entry *tbl;

	// to hold the averges
	struct nm_entry avgs;
	// to calculate average
	int count;
	
};

struct nm_itr *nm_open(const char *, int);
void nm_destroy(struct nm_itr *);
void nm_cleanup(struct nm_itr *);
int nm_query(struct nm_itr *, uint32_t, uint32_t);
int nm_parse_next(struct nm_itr *, struct nm_entry *);

struct nm_tbl *nm_tbl_create();
void nm_tbl_destroy(struct nm_tbl *);
int nm_tbl_add(struct nm_tbl *, struct nm_entry, int only_average);
int nm_tbl_slurp(struct nm_tbl *, struct nm_itr *);
struct nm_entry *nm_tbl_get(struct nm_tbl *, uint32_t);

#endif
