bin_PROGRAMS = gvm
gvm_SOURCES = bam_macro_funcs.h bam_multi_itr.h bedparse.h nmparse.h utils.h bam_mate_table.c base_seq_repr.c config.h ref_seq.c variant_table.c bam_mate_table.h base_seq_repr.h gvm.c ref_seq.h variant_table.h bam_multi_itr.c bedparse.c nmparse.c report.h

BASE_FLAGS = -I./lib/uthash/src/ -Wall -Wextra -pedantic --std=c11  -DHASH_FUNCTION=HASH_BER

gvm_CPPFLAGS = $(BASE_FLAGS) -g -O3 -DNDEBUG

