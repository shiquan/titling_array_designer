PROG=    generate_oligos merge_oligos

all: mk $(PROG)

CC       = gcc
CFLAGS   = -Wall -Wc++-compat -O2
CFLAGS_DEBUG   = -Wall -Wc++-compat -O0 -g
DFLAGS   = -lz -pthread
INCLUDES = -I . -I htslib-1.3.1/ -I src

all:$(PROG)

HTSDIR = htslib-1.3.1
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip
TABIX  = $(HTSDIR)/tabix
HTSLIB_LDFLAGS = $(HTSLIB_static_LDFLAGS)
HTSLIB_LIBS = $(HTSLIB_static_LIBS)

ifeq "$(shell uname -s)" "Darwin"
DYNAMIC_FLAGS = -Wl,-export_dynamic
else
DYNAMIC_FLAGS = -rdynamic
endif

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

mk: clean $(HTSLIB)
	-mkdir -p bin

generate_oligos: version.h
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$@ $(DFLAGS) src/bed_utils.c src/generate_oligos.c $(HTSLIB)

generate_oligos_debug: version.h
	$(CC) $(CFLAGS_DEBUG) $(INCLUDES) -o bin/$@ $(DFLAGS) src/bed_utils.c src/generate_oligos.c $(HTSLIB)

merge_oligos:
	$(CC) $(CFLAGS) $(INCLUDES) -o bin/$@ $(DFLAGS) src/merge_oligos.c  $(HTSLIB)

debug: mk generate_oligos_debug

clean: 
	-rm -f gmon.out *.o *~ $(PROG) version.h  
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM *.bed bin

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
