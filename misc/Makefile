PROG=    generate_oligos

all: $(PROG)

CC       = gcc
CFLAGS   = -g -Wall -Wc++-compat -O2
DFLAGS   = -lz -lhts


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



generate_oligos: version.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(DFLAGS) bed_utils.c generate_oligos.c


clean: 
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM *.bed

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
