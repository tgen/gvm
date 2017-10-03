// variant table data structures

#include <stdlib.h>
#include <stdint.h>

#include "uthash.h"

#include "variant_table.h"



void destroy_variant_counts(struct variant_counts *vc)
{
	struct variant_counts *el, *tmp;
	HASH_ITER(hh, vc, el, tmp) {
		HASH_DEL(vc, el);
		free(el);
	}
}

void destroy_variant_table(struct variant_table *vt)
{
	struct variant_table *el, *tmp;
	HASH_ITER(hh, vt, el, tmp) {
		destroy_variant_counts(el->counts);
		HASH_DEL(vt, el);
		free(el);
	}

}

