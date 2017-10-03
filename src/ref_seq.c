#include <stdint.h>
#include <stdlib.h>

#include "ref_seq.h"

char ref_seq_get(struct ref_seq ref_seq_info, uint32_t offset)
{
	if (ref_seq_info.data == NULL) {
		return 'N';
	}

	return ref_seq_info.data[offset - ref_seq_info.offset];
}
