
#include <stdint.h>

#include "htslib/sam.h"

#include "base_seq_repr.h"
#include "ref_seq.h"

char base_to_char(uint8_t base)
{
	return seq_nt16_str[base];
}

char b5base_to_char(uint8_t base)
{
	if (base > 4) return 'M';
	return "NACGT"[base];
}

inline
uint32_t int_to_b5base(int next)
{
	switch(next) {
	case 1: return 1;
	case 2: return 2;
	case 4: return 3;
	case 8: return 4;
	case 15: return 0;
	default: return 0;
	}
}

inline
uint32_t char_to_b5base(char next)
{
	switch(next) {
	case 'N': return 0;
	case 'A': return 1;
	case 'C': return 2;
	case 'G': return 3;
	case 'T': return 4;
	default: return 0;
	}
}

uint32_t str_to_b5seq(char *bstr)
{
	uint32_t result = 0;
	char *c = bstr;
	while (*c) {
		result = seq_append(result, char_to_b5base(*c));
		c++;
	}

	return result;
}

inline
uint32_t seq_append(uint32_t seq, uint32_t base_val)
{
	if (seq > INDEL_REPR_CUTOFF) {
		return seq + 1;
	}

	uint64_t attempt = (uint64_t) seq * 5 + base_val;
	if (attempt > INDEL_REPR_CUTOFF) {
		return INDEL_REPR_CUTOFF + INDEL_REPR_LEN_CUTOFF;
	}

	return attempt;
}


uint32_t seq_append2(uint32_t seq, uint32_t o_seq)
// what a function!
// but does it work?
{
	uint32_t tmp, seq_len, o_seq_len;
	if (o_seq > INDEL_REPR_CUTOFF && seq > INDEL_REPR_CUTOFF) {
		return (o_seq - INDEL_REPR_CUTOFF) + seq;
	}

	if (o_seq <= INDEL_REPR_CUTOFF && seq <= INDEL_REPR_CUTOFF) {
		seq_len = 0;
		o_seq_len = 0;
		tmp = seq;
		while (tmp > 0) {
			tmp /= 5;
			seq_len++;
		}

		tmp = o_seq;
		while (tmp > 0) {
			tmp /= 5;
			o_seq_len++;
		}

		if (seq_len + o_seq_len > INDEL_REPR_LEN_CUTOFF) {
			return INDEL_REPR_CUTOFF + seq_len + o_seq_len - 1;
		} else {
			while (o_seq_len-- > 0) seq *= 5;

			return seq + o_seq;
		}
	}

	if (o_seq > INDEL_REPR_CUTOFF) {
		while (seq > 0) {
			seq /= 5;
			o_seq++;
		}

		return o_seq;
	} else {
		return seq_append2(o_seq, seq);
	}
}

inline
uint32_t seq_appendi(uint32_t seq, int next)
{
	int base_val = int_to_b5base(next);
	return seq_append(seq, base_val);
}

inline
uint32_t seq_appendc(uint32_t seq, char next)
{
	int base_val = char_to_b5base(next);
	return seq_append(seq, base_val);
}

uint32_t collect_seqi(uint8_t *data, uint32_t start, uint32_t len, uint32_t current)
{
	uint32_t i;
	for (i = start; i < start + len; i++) {
		current = seq_appendi(current, bam_seqi(data, i));
	}

	return current;
}


uint32_t collect_seqc(struct ref_seq ref_seq_info, uint32_t start, uint32_t len, uint32_t current)
{
	uint32_t i;
	for (i = start; i < start + len; i++) {
		current = seq_appendc(current, ref_seq_get(ref_seq_info, i));
	}

	return current;
}
