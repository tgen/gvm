#include <stdio.h>
#include <string.h>

#include "bedparse.h"
#include "utils.h"

int bedf_read_single(	FILE *bedf,
			char *chr,
			size_t chr_size,
			uint32_t *start,
			uint32_t *end )
{
	int dummy;

	char line[4096];
	char format[128];

	if (!fgets(line, sizeof(line), bedf)) {
		return EOF;
	}

	snprintf(format, sizeof(format), "%%%ds\t%%d\t%%d\t%%d", (int) chr_size - 1);
	return sscanf(line, format, chr, start, end, &dummy);
}

int bedf_forall_region_chr(	void *context,
				const char *bed_fn, 
				const char *chromosome,
				region_handle_func reg_func )
{
	char chr[30];
	uint32_t start, end;
	uint32_t index;

	FILE *bedf = fopen(bed_fn, "r");
	if (!bedf) {
		err_printf("Could not open %s: ", bed_fn);
		perror("");
		return 0;
	}

	index = 0;
	while (bedf_read_single(bedf, chr, sizeof(chr), &start, &end) != EOF) {
		if (strcmp(chr, chromosome) == 0) {
			reg_func(context, index, start, end);
			index++;
		}
	}

	return 1;

}
