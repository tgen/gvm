#ifndef GVM_BASE_SEQ_REPR_H
#define GVM_BASE_SEQ_REPR_H

#include <stdint.h>

#include "ref_seq.h"

#define INDEL_REPR_LEN_CUTOFF 13
#define INDEL_REPR_CUTOFF 1220703125 - 1 /* 5^13 - 1 */

char base_to_char(uint8_t base);
char b5base_to_char(uint8_t base);
uint32_t int_to_b5base(int next);
uint32_t char_to_b5base(char next);
uint32_t seq_append(uint32_t seq, int base_val);
uint32_t seq_append2(uint32_t seq, uint32_t o_seq);
uint32_t seq_appendi(uint32_t seq, int next);
uint32_t seq_appendc(uint32_t seq, char next);
uint32_t collect_seqi(uint8_t *data, uint32_t start, uint32_t len, uint32_t current);
uint32_t collect_seqc(struct ref_seq ref_seq_info, uint32_t start, uint32_t len, uint32_t current);

#endif
