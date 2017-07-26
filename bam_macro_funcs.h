#ifndef _MACRO_FUNCS_H
#define _MACRO_FUNCS_H

/* These macros turn sam.h macros into functions so that I can use them from
 * the debugger. This should not be included in production. */

#define GEN_BAM_MACRO_FUNC(RET, M) RET m_##M(bam1_t *b) { return M(b); }

GEN_BAM_MACRO_FUNC(uint8_t *, bam_get_qual)
GEN_BAM_MACRO_FUNC(char *, bam_get_qname)
GEN_BAM_MACRO_FUNC(uint8_t *, bam_get_seq)

uint8_t m_bam_seqi(bam1_t *bam, int i) { return bam_seqi(bam_get_seq(bam), i); }

#endif

