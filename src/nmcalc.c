#include "nmcalc.h"

void nmcalc(struct context *context, struct variant_table *v, struct nm_entry *out)
{
	(void) context;
	out->pos = v->offset;
	out->norm_read_depth = v->read_count_pass;

	out->prob_map_err = 0;
	out->read_pass = 0;
	out->a_or_b = 0;
}
