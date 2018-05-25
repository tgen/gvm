#ifndef GVM_BAM_MULTI_ITR_H
#define GVM_BAM_MULTI_ITR_H

#include <stdint.h>

#include <htslib/sam.h>
#include <uthash.h>

#include "bam_mate_table.h"
#include "variant_table.h"
#include "nmparse.h"

struct bam_single_itr {
	htsFile *f;
	hts_idx_t *idx;
	hts_itr_t *itr;
	bam_hdr_t *hdr;
	struct nm_itr *nmi;
	struct nm_tbl *nmt;
	struct bam_mate_table *bmt;
	bam1_t *buf;
	struct variant_table *vtable;

	int used;
};

// This is to ensure that each normal metrics file is used only once
struct nm_itr_tbl {
    const char *file_name;
    struct nm_itr *itr;
    UT_hash_handle hh;
};

struct bam_multi_itr {
	uint32_t num_iters:16, num_done:16;

	struct bam_single_itr *itr_list;
        struct nm_itr **nm_itr_list;
        struct nm_tbl **nm_tbl_list;

        struct nm_itr_tbl *nms;
};

void bmi_cleanup(struct bam_multi_itr *bmi);
void bmi_destroy(struct bam_multi_itr *bmi);
struct bam_multi_itr *bmi_create(const char *bam_files_fn, const char *default_normal_base_path, const char *chromosome);
int bmi_get_tid(struct bam_multi_itr *bmi, const char *target_name);
int bmi_query(struct bam_multi_itr *bmi, const char *target_name, uint32_t start, uint32_t end);
int bmi_get_begin(struct bam_multi_itr *bmi);
int bmi_get_end(struct bam_multi_itr *bmi);
int bmi_finished(struct bam_multi_itr *bmi);
int bmi_next(struct bam_multi_itr *bmi, bam1_t **out);


#endif
