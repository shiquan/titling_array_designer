// This is a updated version of bedutil.h from bedutils program.
// Copyright shiquan@link.cuhk.edu.hk
//

#ifndef BED_UTILS_HEADER
#define BED_UTILS_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <kstring.h>
#include <khash.h>

#ifndef KSTRINT_INIT
#define KSTRINT_INIT { 0, 0, 0}
#endif

#define MEMPOOL_LINE 10000 // todo: memory management

// flag of bed_file struct
// bits offset rule : right first
//                           bed file is empty or not
//                          /
// uint8_t : | | | | | | | |
//                 | | |  \
//                 | | |   part of file is cached
//                 | | |_  has extra data (more than 3 cols in this bed file)
//                 | |___  sorted
//                 |_____  merged
#define BED_FILE_EMPTY  1
#define BED_FILE_CACHED 2 // when finished, should clear this bit
#define BED_FILE_EXTRA  4
#define BED_FILE_SORTD  (1<<3)
#define BED_FILE_MERGED (1<<4)

struct bed_line_core {
    int start; // 0-based
    int end; // 1-based
    kstring_t data; //data should be inited in very first early time
};

struct bed_chrom {
    uint8_t flag;
    int m, l;
    struct bed_line_core *a;
    int id; // name id, usually for chromosomes or contigs
    uint32_t length;    
};

typedef struct bed_chrom* bed_chrom_point;
KHASH_MAP_INIT_STR(reg, struct bed_chrom_point)
typedef kh_reg_t reghash_t;

struct bed_file {
    int errno;
    uint8_t flag;
    int l_names, m_names;
    char **names;
    reghash_t *hash;
    
    uint32_t regions_count; // how many gapped regions in this bed file
    uint32_t length_total; // total length of all these chromosomes
    gzFile fp; // Hold the handler of file. For very big files, read part of data into memory first and merge and read other parts
};

struct bed_readers {
    int n_files;
    char **file_names;
    //int sort_required;
    struct bed_file **files;    
};

int bed_readers_add_file(struct bed_readers *files, const char *fname);
void bed_readers_remove_file(struct bed_readers*files, int i);
struct bed_readers *bed_readers_init();
void bed_readers_destroy(struct bed_readers *files);

#endif
