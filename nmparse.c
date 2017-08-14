#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "nmparse.h"

struct nm_itr *nm_open(const char *fn, int tid)
{
	htsFile *htsfp = hts_open(fn, "r");
	if (htsfp == NULL) goto fail;

	tbx_t *t = tbx_index_load(fn);
	if (t == NULL) goto fail1;

	struct nm_itr *nmi = malloc(sizeof(*nmi));
	if (nmi == NULL) goto fail2;

	nmi->htsfp = htsfp;
	nmi->t = t;
	nmi->itr = NULL;
	nmi->tid = tid;

	return nmi;

fail2:
	tbx_destroy(t);
fail1:
	hts_close(htsfp);
fail:
	return NULL;

	
}

void nm_destroy(struct nm_itr *nmi)
{
	if (nmi == NULL) return;

	tbx_destroy(nmi->t);

	hts_close(nmi->htsfp);

	free(nmi);
}

void nm_cleanup(struct nm_itr *nmi)
{
	tbx_itr_destroy(nmi->itr);
}

int nm_query(struct nm_itr *nmi, uint32_t start, uint32_t end)
{
	assert(nmi != NULL /* nm_query called with null nm_itr */);

	if (nmi->t == NULL) return -1;
	nmi->itr = tbx_itr_queryi(nmi->t, nmi->tid, start, end);
	if (nmi->itr == NULL) return -2;

	return 1;
}

static
int nm_parse_line(const char *line, struct nm_entry *out)
{
	int res = sscanf(line, "%*s\t%u\t%lf\t%lf\t%lf\t%lf",
			&out->pos,
			&out->norm_read_depth,
			&out->prob_map_err,
			&out->read_pass,
			&out->a_or_b);
	if (isnan(out->norm_read_depth)) out->norm_read_depth = 0;
	if (isnan(out->prob_map_err)) out->prob_map_err = 0;
	if (isnan(out->read_pass)) out->read_pass = 0;
	if (isnan(out->a_or_b)) out->a_or_b = 0;

	return res;
}

int nm_parse_next(struct nm_itr *nmi, struct nm_entry *out)
{
	kstring_t str = {0,0,0};
	int res = tbx_itr_next(nmi->htsfp, nmi->t, nmi->itr, &str);

	if (res > 0) nm_parse_line(str.s, out);
	free(ks_release(&str));
	return res;
}

struct nm_tbl *nm_tbl_create()
{

	struct nm_tbl *tbl = calloc(1, sizeof(*tbl));
	if (tbl == NULL) return NULL;

	return tbl;
}

void nm_tbl_destroy(struct nm_tbl *nmt)
{
	struct nm_entry *ent, *tmp;
	HASH_ITER(hh, nmt->tbl, ent, tmp) {
		HASH_DEL(nmt->tbl, ent);
		free(ent);
	}

	free(nmt);

}

#define UPDATE_AVERAGE(CUR, NEXT, TOTAL) CUR = (CUR) * ((TOTAL) - 1) / (TOTAL) + (NEXT) / (TOTAL)
int nm_tbl_add(struct nm_tbl *nmt, struct nm_entry ent)
{
	struct nm_entry *heap_ent = malloc(sizeof(*heap_ent));
	if (heap_ent == NULL) return -1;

	memcpy(heap_ent, &ent, sizeof(ent));

	HASH_ADD_INT(nmt->tbl, pos, heap_ent);

	nmt->count++;

	// update averages
	UPDATE_AVERAGE(nmt->avgs.norm_read_depth, ent.norm_read_depth, nmt->count);
	UPDATE_AVERAGE(nmt->avgs.prob_map_err, ent.prob_map_err, nmt->count);
	UPDATE_AVERAGE(nmt->avgs.read_pass, ent.read_pass, nmt->count);
	UPDATE_AVERAGE(nmt->avgs.a_or_b, ent.a_or_b, nmt->count);

	return 1;

}
#undef UPDATE_AVERAGE

int nm_tbl_slurp(struct nm_tbl *nmt, struct nm_itr *nmi)
{
	assert(nmi->itr != NULL /* Need to call nm_query before nm_tbl_slurp */);
	struct nm_entry ent = {0};
	while (nm_parse_next(nmi, &ent) >= 0) {
		int res = nm_tbl_add(nmt, ent);
		if (res != 1) return res;
	}

	return 1;
}

struct nm_entry *nm_tbl_get(struct nm_tbl *nmt, uint32_t pos)
{
	assert(nmt != NULL /* nm_tbl_get with null nm_tbl */);

	struct nm_entry *out;
	HASH_FIND_INT(nmt->tbl, &pos, out);
	return out;
}

/* REMOVE THIS LINE TO BRING BACK MAIN
int main()
{


	const char *fn = "/scratch/rhalperin/TumorOnlyCaller/Controls/S5U_200Controls_051117_14.txt.gz";

	struct nm_itr *nmi = nm_open(fn, 0);
	if (!nmi) {
		perror("nm_open");
		return 1;
	}

	struct nm_tbl *nmt = nm_tbl_create();
	if (!nmt) {
		perror("nm_tbl_create");
		return 1;
	}

	nm_query(nmi, 0, 19434040);

	nm_tbl_slurp(nmt, nmi);


	nm_tbl_destroy(nmt);
	nm_destroy(nmi);

	return 0;

}
// */
