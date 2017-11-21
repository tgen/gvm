#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <uthash.h>

#include "utils.h"
#include "base_seq_repr.h"

#include "htslib/synced_bcf_reader.h"

struct cosm_count_table {
	// the following two fields form the key
	uint32_t pos;
	uint32_t b5alt;

	// value
	uint32_t count;

	UT_hash_handle hh;
};

static struct cosm_count_table *cct;

int32_t get_cosmic_count(uint32_t pos, uint32_t b5alt)
{
	uint32_t key_len;
	struct cosm_count_table key_holder;
	struct cosm_count_table *ent;

	memset(&key_holder, 0, sizeof(key_holder));

	key_holder.pos = pos;
	key_holder.b5alt = b5alt;

	key_len =   offsetof(struct cosm_count_table, b5alt)
		  + sizeof(b5alt)
		  - offsetof(struct cosm_count_table, pos);

	HASH_FIND(hh, cct, &key_holder.pos, key_len, ent);

	if (ent == NULL) return -1;

	return ent->count;
}

int load_cosmic_table(const char *vcf_fname, const char *chr)
{
	int res;
	uint32_t pos, b5alt, count;
	uint32_t key_len;
	struct cosm_count_table *ent;
	bcf1_t *line;
	bcf_hdr_t *header;
	bcf_info_t *info;
	bcf_srs_t *rs = bcf_sr_init();

	bcf_sr_set_regions(rs, chr, 0);

	res = bcf_sr_add_reader(rs, vcf_fname);
	header = rs->readers[0].header;

	cct = NULL;

	key_len =   offsetof(struct cosm_count_table, b5alt)
		  + sizeof(b5alt)
		  - offsetof(struct cosm_count_table, pos);

	if (res == 0) {
		err_printf("error: %s\n", bcf_sr_strerror(rs->errnum));
		return 0;
	}

	while (bcf_sr_next_line(rs)) {
		if (!bcf_sr_has_line(rs, 0)) continue;

		line = bcf_sr_get_line(rs, 0);
		info = bcf_get_info(header, line, "CNT");
		pos = line->pos + 1;
		if (!info) continue;

		// TODO: enforce the condition below more strictly!
		if (line->n_allele != 2) continue;

		b5alt = str_to_b5seq(line->d.allele[1]); // [1] is the ALT
		if (get_cosmic_count(pos, b5alt) > -1) continue;

		count = info->v1.i;
		ent = malloc(sizeof(struct cosm_count_table));
		ent->pos = pos;
		ent->b5alt = b5alt;
		ent->count = count;

		HASH_ADD(hh, cct, pos, key_len, ent);
	}


	bcf_sr_destroy(rs);

	return 1;

}
