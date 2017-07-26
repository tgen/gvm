CC=gcc
UNAME=$(shell uname)

UTHASH=lib/uthash
UTHASH_HFUNC=BER

INCS=-I$(UTHASH)/include
LIBS=-lhts -lyaml

HEADERS=utils.h report.h nmparse.h bedparse.h ref_seq.h base_seq_repr.c bam_multi_itr.h variant_table.h bam_mate_table.h
SRCS=gvm.c nmparse.c bedparse.c ref_seq.c base_seq_repr.c bam_multi_itr.c variant_table.c bam_mate_table.c

EXENAME=gvm

EXE=$(DIRPREFIX)/$(EXENAME)
OBJS=$(addprefix $(DIRPREFIX)/, $(SRCS:.c=.o))

ifeq ($(UNAME), Darwin)
  DIRPREFIX=local
else
  DIRPREFIX=prod
endif


CFLAGS = -Wall -Wextra -pedantic --std=c11
CFLAGS += -DHASH_FUNCTION=HASH_$(UTHASH_HFUNC)

all: release

debug: CFLAGS += -g -ggdb -DDEBUG -DCLEANUP -O0
debug: $(EXE)

# the NDEBUG definition is to eliminate calls to assert
release: CFLAGS += -O3 -g -DNDEBUG
release: $(EXE)

$(DIRPREFIX)/%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCS) -c -o $@ $<

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) $(INCS) -o $(EXE) $^

clean:
	rm -f $(EXE) $(OBJS)

.PHONY: clean
