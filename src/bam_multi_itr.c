#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>


#include <htslib/sam.h>
#include <uthash.h>

#include "utils.h"
#include "variant_table.h"
#include "bam_multi_itr.h"

void bmi_cleanup(struct bam_multi_itr *bmi)
{
	int i;
	bmi->num_done = 0;

	struct bam_single_itr *s;

	for (i = 0; i < bmi->num_iters; i++) {
		s = &bmi->itr_list[i];
		if (s->itr) hts_itr_destroy(s->itr);
		s->itr = NULL;
		destroy_variant_table(s->vtable);
		s->vtable = NULL;
		s->used = 1;

		if (bmi->nm_itr_list[i]) {
			nm_cleanup(bmi->nm_itr_list[i]);
		}
		if (bmi->nm_tbl_list[i]) {
			nm_tbl_destroy(bmi->nm_tbl_list[i]);
			bmi->nm_tbl_list[i] = NULL;
		}
	}

}

static
void bmi_nm_slurp(struct bam_multi_itr *bmi, uint32_t start, uint32_t end)
{
	int i, j;
	struct nm_itr *itr;
	struct nm_tbl *tbl;
	for (i = 0; i < bmi->num_iters; i++) {
		itr = bmi->nm_itr_list[i];
		tbl = bmi->nm_tbl_list[i];
		if (tbl) nm_tbl_destroy(tbl);

		if (itr == NULL) continue;

		nm_query(itr, start, end);
		tbl = nm_tbl_create();
		nm_tbl_slurp(tbl, itr);

		bmi->nm_tbl_list[i] = tbl;

		// now look for iterators in bmi that use this normal metrics
		// file and set their table to `tbl`
		for (j = 0; j < bmi->num_iters; j++) {
			if (bmi->itr_list[j].nmi == itr) {
				bmi->itr_list[j].nmt = tbl;
			}
		}
	}
}

void bmi_destroy(struct bam_multi_itr *bmi)
{
	int i;

	bmi_cleanup(bmi);

	struct bam_single_itr *s;
	for (i = 0; i < bmi->num_iters; i++) {
		s = &bmi->itr_list[i];
		if (s->idx) hts_idx_destroy(s->idx);
		if (s->hdr) bam_hdr_destroy(s->hdr);
		if (s->f) sam_close(s->f);
		if (s->buf) bam_destroy1(s->buf);
		if (bmi->nm_itr_list[i]) nm_destroy(bmi->nm_itr_list[i]);
	}

	free(bmi->itr_list);
	free(bmi->nm_itr_list);
	free(bmi->nm_tbl_list);
	free(bmi);
}

static
struct nm_itr *load_normal_metrics(const char *base_path, const char *chromosome)
{
	char *normal_path = malloc(strlen(base_path) +
				strlen(chromosome) +
				strlen(".txt.gz") + 1);

	struct nm_itr *nmi;

	strcpy(normal_path, base_path);
	strcat(normal_path, chromosome);
	strcat(normal_path, ".txt.gz");

	nmi = nm_open(normal_path, 0);

	free(normal_path);

	return nmi;
}

static
struct nm_itr *bmi_nm_open(	struct bam_multi_itr *bmi,
				const char *base_path,
				const char *chromosome	)
{
	unsigned int i;
	struct nm_itr_tbl *entry;
	HASH_FIND_STR(bmi->nms, base_path, entry);

	if (entry == NULL) {
		entry = malloc(sizeof(struct nm_itr_tbl));
		entry->file_name = base_path;
		entry->itr = load_normal_metrics(base_path, chromosome);
		HASH_ADD_STR(bmi->nms, file_name, entry);

		for (i = 0; bmi->nm_itr_list[i]; i++);
		bmi->nm_itr_list[i] = entry->itr;
	}

	return entry->itr;
}

struct bam_multi_itr *bmi_create(	const char *bam_files_fn,
					const char *default_normal_base_path,
					const char *chromosome	)
{
	struct bam_multi_itr *bmi = NULL;
	FILE *blfp;
	bam_hdr_t *hdr;
	hts_idx_t *idx;

	char buf[1024];
	char *normal_base_path;

	htsFile *f;
	/* first find out how many files are there */

	int count = 0, i;

	blfp = fopen(bam_files_fn, "r");
	if (blfp == NULL) {
		err_printf("could not open \"%s\" for reading: ", bam_files_fn);
		perror("");
		goto fail_already_clean;
	}

	while ( fgets(buf, sizeof(buf), blfp) != NULL ) {
		count++;
	}
	fclose(blfp);

	if (count == 0) {
		err_printf("refusing to operate on empty bamlist\n");
		goto fail_already_clean;
	}

	bmi = malloc(sizeof(struct bam_multi_itr));
	if (bmi == NULL) {
		err_printf("failed to allocate memory\n");
		goto fail_already_clean;
	}
	
	bmi->num_iters = count;
	bmi->num_done = 0;
	bmi->itr_list = calloc(count, sizeof(struct bam_single_itr));
	bmi->nm_itr_list = calloc(count, sizeof(struct nm_itr));
	bmi->nm_tbl_list = calloc(count, sizeof(struct nm_tbl));
	bmi->nms = NULL;

	if (bmi->itr_list == NULL || bmi->nm_itr_list == NULL
			|| bmi->nm_tbl_list == NULL) {

		err_printf("failed to allocate memory\n");
		free(bmi->itr_list);
		free(bmi->nm_itr_list);
		free(bmi->nm_tbl_list);
		free(bmi);
		goto fail_already_clean;
	}

	i = 0;
	blfp = fopen(bam_files_fn, "r");
	while ( fgets(buf, sizeof(buf), blfp) != NULL ) {
		buf[strlen(buf)-1] = '\0';

		// If there's a comma, set it to \0 to split
		// the string
		normal_base_path = strchr(buf, ',');
		if (normal_base_path == NULL) {
			normal_base_path = (char *) default_normal_base_path;
		} else {
			*normal_base_path = '\0';
			normal_base_path++;
		}

		f = sam_open(buf, "rb");
		if (f == NULL) {
			err_printf("failed to open BAM file\n");
			goto fail;
		}

		hdr = bam_hdr_read(f->fp.bgzf);
		idx = sam_index_load(f, buf);
		if (idx == NULL) {
			err_printf("failed to load BAM index for %s\n", buf);
			goto fail;
		}

		bmi->itr_list[i].f = f;
		bmi->itr_list[i].idx = idx;
		bmi->itr_list[i].itr = NULL;
		bmi->itr_list[i].hdr = hdr;
		bmi->itr_list[i].buf = bam_init1();
		bmi->itr_list[i].vtable = NULL;
		bmi->itr_list[i].used = 1;
		bmi->itr_list[i].nmi = bmi_nm_open(bmi, normal_base_path, chromosome);

		i++;
	}
	fclose(blfp);

	return bmi;

fail:
	fclose(blfp);
	bmi_destroy(bmi);
fail_already_clean:
	return NULL;
}

inline
int bmi_get_tid(struct bam_multi_itr *bmi, const char *target_name)
{
	if (bmi->num_iters < 1) {
		return -1;
	}

	return bam_name2id(bmi->itr_list[0].hdr, target_name);

}

int bmi_query(struct bam_multi_itr *bmi, const char *target_name, uint32_t start, uint32_t end)
{

	int i, tid;
	struct bam_single_itr *s;

	tid = bmi_get_tid(bmi, target_name);

	for (i = 0; i < bmi->num_iters; i++) {
		s = &bmi->itr_list[i];
		s->itr = sam_itr_queryi(s->idx, tid, start, end);
		if (s->itr == NULL) {
			err_printf("could not query region %s-%d-%d in BAM file #%d",
					target_name, start, end, i);
			return 0;
		}
	}

	bmi_nm_slurp(bmi, start - 1, end);

	return 1;
}

// returns the minimum begin position among the iterators
int bmi_get_begin(struct bam_multi_itr *bmi)
{
	int min_beg = INT_MAX, i;
	for (i = 0; i < bmi->num_iters; i++) {
		if (bmi->itr_list[i].itr->beg < min_beg) {
			min_beg = bmi->itr_list[i].itr->beg;
		}
	}

	return min_beg + 1;

}

// returns the maximum end position among the iterators
int bmi_get_end(struct bam_multi_itr *bmi) {
	int max_end = 0, i;
	for (i = 0; i < bmi->num_iters; i++) {
		if (bmi->itr_list[i].itr->end > max_end) {
			max_end = bmi->itr_list[i].itr->end;
		}
	}

	return max_end;
}

static
void bmi_ensure(struct bam_multi_itr *bmi)
{
	int i, result = 0;
	struct bam_single_itr s;
	for (i = 0; i < bmi->num_iters; i++) {
		s = bmi->itr_list[i];
	
		if (s.itr->finished || !s.used) {
			continue;
		}

		result = sam_itr_next(s.f, s.itr, s.buf);

#ifdef HTS_HAS_ENHANCED_ENDPOS
		// this check is here to make sure the endpos cache is being cleared.
		// I am NOT checking that the read's end position is zero.
		assert(s.buf->core.endpos == 0);
#endif

		if (result == -1) {
			bmi->num_done++;
		} else {
			bmi->itr_list[i].used = 0;
		}

	}

}

int bmi_finished(struct bam_multi_itr *bmi)
{
	return bmi->num_done == bmi->num_iters;
}

int bmi_next(struct bam_multi_itr *bmi, bam1_t **out)
{
	int i, min_offset, min_index;
	struct bam_single_itr s;

	*out = NULL;
	if (bmi_finished(bmi)) {
		return -1;
	}

	bmi_ensure(bmi);
	min_offset = INT_MAX;
	min_index = -1;

	for (i = 0; i < bmi->num_iters; i++) {
		s = bmi->itr_list[i];
		if (s.used) continue;

		if (s.buf->core.pos < min_offset) {
			min_offset = s.buf->core.pos;
			min_index = i;
		}
	}

	if (min_index >= 0) {
		s = bmi->itr_list[min_index];
		*out = s.buf;
		bmi->itr_list[min_index].used = 1;
		return min_index;
	}

	return -2;
}

