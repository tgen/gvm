SUBDIRS = gengetopt

bin_PROGRAMS = gvm
gvm_SOURCES = bam_macro_funcs.h bam_multi_itr.h bedparse.h nmparse.h nmcalc.h utils.h bam_mate_table.c base_seq_repr.c config.h ref_seq.c variant_table.c bam_mate_table.h base_seq_repr.h gvm.c ref_seq.h variant_table.h bam_multi_itr.c bedparse.c nmparse.c nmcalc.c report.h cigar.c cigar.h cosmvcf.c cosmvcf.h gengetopt/cmdline.c gengetopt/cmdline.h

BASE_FLAGS = -Wall -Wextra -pedantic --std=c11

gvm_CFLAGS = $(BASE_FLAGS)
gvm_CPPFLAGS = -DHASH_FUNCTION=HASH_BER -I@abs_top_srcdir@/lib/uthash/src/ 

gvm-cosmvcf.o: cosmvcf.c cosmvcf.h
	$(CC) $(CFLAGS:-Og=) $(gvm_CFLAGS:-Og=) -O3 $(gvm_CPPFLAGS) $< -c -o $@

clean-local:
	rm -f $(srcdir)/gengetopt/cmdline.[ch]
