#ifndef GVM_BEDPARSE_H
#define GVM_BEDPARSE_H

#include <stdio.h>
#include <stdint.h>

#include "utils.h"

struct context;

typedef int (*region_handle_func)(void *, uint32_t index, uint32_t start, uint32_t end);

int bedf_read_single(	FILE *bedf,
			char *chr,
			size_t chr_size,
			uint32_t *start,
			uint32_t *end );

int bedf_forall_region_chr(	void *context,
				const char *bed_fn, 
				const char *chromosome,
				region_handle_func reg_func );

#endif
