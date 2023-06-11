IDIR =./include
CC=gcc
# note that usually use -O2 but for valgrind debugging use -O0 which is slower but more accurate
CFLAGS=-g -Wall -O2 -Wno-return-type -Wno-unused-variable -Wno-unused-function -I$(IDIR)
LFLAGS= -lz -lm 

ODIR=./obj/
SDIR=./src/
BDIR=./bin/

LIBS=-lm


#PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
    PREFIX := /usr/local/bin
endif

install: $(BDIR)fast-adapter-trimming 
	install -d $(DESTDIR)$(PREFIX)
	install -m 755 $^ $(DESTDIR)$(PREFIX)

BASE=seqs_to_ints rob_dna store-primers create-snps read_primers search-adapter-file hash primer-match-counts
FAT=$(BASE) paired_end_search fast_search
fatobj := $(addsuffix .o, $(addprefix $(ODIR), $(FAT)))
objects := $(fatobj)


$(objects): $(ODIR)%.o: $(SDIR)%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@ $(FLAGS)


$(BDIR)fast-adapter-trimming: $(fatobj)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

fast-adapter-trimming: $(BDIR)fast-adapter-trimming

EXEC=fast-adapter-trimming
all: $(addprefix $(BDIR), $(EXEC))



.PHONY: clean

clean:
	rm -fr bin/ obj/


