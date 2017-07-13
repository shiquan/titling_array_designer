PROG=    generate_oligos merge_oligos

all: mk $(PROG)

CC       = gcc
CFLAGS   = -Wall -Wc++-compat -O2
DFLAGS   = -lz -lhts
INCLUDES = -I . -I src/

all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION = 0.01
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
DOC_VERSION :=  $(shell git describe --always)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define OLIGOS_VERSION "$(PACKAGE_VERSION)"' > $@

mk:
	-mkdir bin

generate_oligos: version.h
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$@ $(DFLAGS) src/bed_utils.c src/generate_oligos.c src/faidx.c src/kstring.c src/bgzf.c

merge_oligos:
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$@ $(DFLAGS) src/merge_oligos.c  src/kstring.c src/bgzf.c

clean: 
	-rm -f gmon.out *.o *~ $(PROG) version.h  
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM *.bed bin

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
