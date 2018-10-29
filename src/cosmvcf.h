#ifndef GVM_COSMVCF_H
#define GVM_COSMVCF_H
#include <stdint.h>

void destroy_cosmic_table();
uint32_t get_cosmic_count(uint32_t pos, uint32_t b5alt);
int load_cosmic_table(const char *vcf_fname, const char *chr);
#endif
