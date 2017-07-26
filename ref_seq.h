#ifndef GVM_REF_SEQ_H
#define GVM_REF_SEQ_H

#include <stdint.h>

struct ref_seq {
	uint32_t offset;
	char *data;
	int len;
};

char ref_seq_get(struct ref_seq ref_seq_info, uint32_t offset);

#endif
