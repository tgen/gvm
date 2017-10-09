#ifndef GVM_NMCALC_H
#define GVM_NMCALC_H

#include "nmparse.h"
#include "variant_table.h"

struct context;


void nmcalc(struct context *context, struct variant_table *v, struct nm_entry *out);

#endif
