bin_PROGRAMS = gvm
gvm_SOURCES = bam_macro_funcs.h bam_multi_itr.h bedparse.h nmparse.h utils.h bam_mate_table.c base_seq_repr.c config.h ref_seq.c variant_table.c bam_mate_table.h base_seq_repr.h gvm.c ref_seq.h variant_table.h bam_multi_itr.c bedparse.c nmparse.c report.h cigar.c cigar.h cmdline.c cmdline.h

cmdline.c cmdline.h: @abs_top_srcdir@/gengetopt/gvm.ggo
	gengetopt -i @abs_top_srcdir@/gengetopt/gvm.ggo --output-dir=@abs_top_srcdir@

BASE_FLAGS = -Wall -Wextra -pedantic --std=c11

gvm_CFLAGS = $(BASE_FLAGS) -g -O3 
gvm_CPPFLAGS = -DHASH_FUNCTION=HASH_BER -I@abs_top_srcdir@/lib/uthash/src/ 


