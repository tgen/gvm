#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>


#include "htslib/sam.h"

#include "utils.h"
#include "variant_table.h"
#include "bam_multi_itr.h"

void bmi_cleanup(struct bam_multi_itr *bmi)
{
	int i;
	bmi->num_done = 0;
	for (i = 0; i < bmi->num_iters; i++) {
		hts_itr_destroy(bmi->itr_list[i].itr);
		bmi->itr_list[i].itr = NULL;
		destroy_variant_table(bmi->itr_list[i].vtable);
		bmi->itr_list[i].vtable = NULL;
		bmi->itr_list[i].used = 1;
	}

}

void bmi_destroy(struct bam_multi_itr *bmi)
{
	int i;

	bmi_cleanup(bmi);

	for (i = 0; i < bmi->num_iters; i++) {
		hts_idx_destroy(bmi->itr_list[i].idx);
		bam_hdr_destroy(bmi->itr_list[i].hdr);
		if (bmi->itr_list[i].f) sam_close(bmi->itr_list[i].f);
		bam_destroy1(bmi->itr_list[i].buf);
	}

	free(bmi->itr_list);
	free(bmi);
}

struct bam_multi_itr *bmi_create(const char *bam_files_fn)
{
	struct bam_multi_itr *bmi = NULL;
	FILE *blfp;
	bam_hdr_t *hdr;
	hts_idx_t *idx;

	char buf[1024];
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
		goto fail;
	}
	
	bmi->num_iters = count;
	bmi->num_done = 0;
	bmi->itr_list = calloc(count, sizeof(struct bam_single_itr));
	if (bmi->itr_list == NULL) {
		err_printf("failed to allocate memory\n");
		goto fail;
	}

	i = 0;
	blfp = fopen(bam_files_fn, "r");
	while ( fgets(buf, sizeof(buf), blfp) != NULL ) {
		buf[strlen(buf)-1] = '\0';
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

		// this check is here to make sure the endpos cache is being cleared.
		// I am NOT checking that the read's end position is zero.
		assert(s.buf->core.endpos == 0);

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

