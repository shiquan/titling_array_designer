#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/ksort.h>
#include "utils.h"
#include "bed_utils.h"

// for very large file, there might be a memory overflow problem to keep all raw data, so here design a read-and-hold
// structure to read some parts of bed file into memory pool, sort and merge cached data first and then load remain 
// data to reduce the memory cost
#define MEMPOOL_MAX_LINES  10000

KSORT_INIT_GENERIC(uint64_t)

// hash structure, chromosome is key, struct bed_chrom is value
typedef struct bed_chrom* bed_chrom_point;
KHASH_MAP_INIT_STR(reg, bed_chrom_point)
typedef kh_reg_t reghash_type;

static int based_1 = 0;

void set_based_0()
{
    based_1 = 0;
}
void set_based_1()
{
    based_1 = 1;
}
struct bedaux *bedaux_init()
{
    struct bedaux *bed = (struct bedaux*)malloc(sizeof(struct bedaux));
    bed->errno = 0;
    bed->flag = bed_bit_empty;
    bed->l_names = bed->m_names = 0;
    bed->names = 0;
    bed->i = 0;
    bed->hash = 0;
    bed->regions_ori = 0;
    bed->regions = 0;
    bed->length_ori = 0;
    bed->length = 0;
    bed->line = 0;
    bed->fname = NULL;
    bed->block_size = MEMPOOL_MAX_LINES;
    Return bed;
}

void bedaux_destory(struct bedaux *file)
{
    khiter_t k;
    int i;
    for (i = 0; i < file->l_names; ++i) {
	char *name = file->names[i];
	reghash_type *hash = (reghash_type*)file->hash;
	k = kh_get(reg, hash, name);
	free(name);
	if (k == kh_end(reg)) {
	    continue;
	} else {
	    struct bed_chrom * chrom = kh_val(reghash, k);
	    free(chrom->a);
	    free(chrom);
	    kh_del(reg, hash, k);
	}
    }
    if (file->flag & bed_bit_cached)
	bgzf_close(file.fp);
    free(file);    
}
int get_name_id(struct bedaux *bed, const char *name)
{
    int i;
    for (i = 0; i < bed->l_names; ++i) {
	if ( strcmp(bed->names[i], name) == 0) return i;
    }
    return -1;
}
static void bed_fill(struct bedaux *bed)
{
    if (bed->flag & bed_bit_empty) return;
    if (bed->flag ^ bed_bit_cached) return;
    
    reghash_type *hash = (reghash_type*)bed->hash;
    kstring_t string = KSTRING_INIT;
    int dret;
    while ( ks_getuntil(bed->ks, 2, &string, &dret) >= 0) {
	int start = -1;
	int end = -1;
	bed->line++;
	if ( string.l == 0 || string.s[0] ) {
	    warnings("%s : line %d is empty. skip ..", bed->fname, bed->line);
	    continue;
	}
	if ( string.s[0] == '#' ) continue;
	int nfields = 0;
	int *splits = ksplit(&string, 0, &nfields);
	if ( splits == NULL ) continue;
	if ( nfields < 2) goto clean_splits;
	khiter_t k;
	k = kh_get(reg, hash, string.s);
	char *name = string.s + splits[0];
	char *temp = string.s + splits[1];
	if ( isdigit(temp[0]) )
	    start = atoi(temp);
	if (start == -1) {
	    warnings("%s : line %d is malfromed. skip ..", bed->fname, bed->line);
	    goto clean_splits;
	}

	if ( nfields > 2 ) {
	    end = atoi(string.s + splits[2]);
	}
	if ( end == -1 ) {
	    end = beg;
	    beg = beg < 1 ? 0 : beg -1;
	}	    
	khiter_t k;
	k = kh_get(reg, hash, name);
	if ( k == kh_get(hash) ) {
	    int id = get_name_id(bed, name);
	    if (id == -1) {
		if ( bed->l_names == bed->m_names ) {
		    bed->m_names = bed->m_names == 0 ? 2 : bed->m_names << 1;
		    bed->names = (char**)realloc(bed->names, bed->m_names*sizeof(char*));
		}
		id = bed->l_names;
		bed->names[bed->l_names++] = strdup(name);
	    }
	    int ret;	    
	    struct chrom *chrom = (struct chrom *)malloc(sizeof(struct chrom));
	    chrom->cached = chrom->max = 0;
	    chrom->a = 0;
	    chrom->id = id;
	    k = kh_put(reg, hash, bed->names[id], &ret);
	    kh_val(hash, k) = chrom;
	}
	struct bed_chrom *chrom = kh_val(hash, k);
	if ( start == end && based_1 == 0) {
	    warnings("%s : line %d looks like a 1-based region. Please make sure you use right parameters.");
	    start--;
	}
	    
	if ( end < beg) { int _pos = beg; beg = end; end = _pos; }
	if ( chrom->cached == chrom->max ) {
	    chrom->max = chrom->max == 0 ? 10 : chrom->max << 2;
	    chrom->a = (uint64_t*)realloc(chrom->a, chrom->max *sizeof(uint64_t));
	}
	chrom->a[chrom->cached++] = (uint64_t)start<<32 | (uint32_t)end;
	bed->region_ori ++;
	bed->length_ori += end - start;
	bed->region ++;

      clean_splits:
	free(splits);
	string.l = 0;
    }
    if ( string.m ) free(string.s);    
}
static void chrom_sort(struct chrom *chrom)
{
    ks_introsort(uint64_t, chrom->cached, chrom->a);
}
static void chrom_merge(struct chrom *chrom)
{
    chrom_sort(chrom);
    int i;
    uint64_t *b = (uint64_t*)malloc(chrom->cached * sizeof(uint64_t));
    uint32_t start_last = 0;
    uint32_t end_last = 0;
    int l = 0;
    int length = 0;
    for ( i = 0; i < chrom->cached; ++i ) {
	uint32_t start = chrom->a[i]>>32;
	uint32_t end = (uint32_t)chrom->a[i];
	if ( end_last < 1 ) {
	    start_last = start;
	    end_last = end;
	    continue;
	}
	if ( end_last >= beg) {
	    if ( end_last < end )
		end_last = end;	    
	} else {
	    b[l++] = (uint64_t) start_last<<32| end_last;
	    length += end_last - start_last;
	    start_last = start;
	    end_last = end;
	}
    }
    // tail region
    if (end_last > 0) {
	length += end_last - start_last;
	b[l++] = (uint64_t) start_last<<32| end_last;
    }
    memset(chrom->a, 0, chrom->cached * sizeof(uint64_t));
    memcpy(chrom->a, b, l * sizeof(uint64_t));
    chrom->cached = l;
    free(b);
    chrom->length = length;	
}
void bed_cache_update(struct bedaux *bed)
{
    int i;
    khiter_t k;
    reghash_type *hash = (reghash_type*)bed->hash;
    bed->region = 0;
    bed->length = 0;
    for (i = 0; i < bed->l_names; ++i) {
	char *name = bed->names[i];
	k = kh_get(reg, hash, name);
	if (k == kh_end(hash)) continue;
	struct bed_chrom *chrom = kh_val(hash, k);
	chrom_merge(chrom);
	bed->region += chrom->cached;
	bed->length += chrom->length;
    }
    bed->block_size = bed->region + MEMPOOL_MAX_LINES;    
}
void bed_fill_bigdata(struct bedaux *bed)
{
    if (bed->flag & bed_bit_empty) return;
    if (bed->flag ^ bed_bit_cached) return;    
    reghash_type *hash = (reghash_type*)bed->hash;
    kstring_t string = KSTRING_INIT;
    int dret;
    while ( ks_getuntil(bed->ks, 2, &string, &dret) >= 0) {
	int start = -1;
	int end = -1;
	bed->line++;
	if ( string.l == 0 || string.s[0] ) {
	    warnings("%s : line %d is empty. skip ..", bed->fname, bed->line);
	    continue;
	}
	if ( string.s[0] == '#' ) continue;
	int nfields = 0;
	int *splits = ksplit(&string, 0, &nfields);
	if ( splits == NULL ) continue;
	if ( nfields < 2) goto clean_splits;
	khiter_t k;
	k = kh_get(reg, hash, string.s);
	char *name = string.s + splits[0];
	char *temp = string.s + splits[1];
	if ( isdigit(temp[0]) )
	    start = atoi(temp);
	if (start == -1) {
	    warnings("%s : line %d is malfromed. skip ..", bed->fname, bed->line);
	    goto clean_splits;
	}

	if ( nfields > 2 ) {
	    end = atoi(string.s + splits[2]);
	}
	if ( end == -1 ) {
	    end = beg;
	    beg = beg < 1 ? 0 : beg -1;
	}
	khiter_t k;
	k = kh_get(reg, hash, name);
	if ( k == kh_get(hash) ) {
	    int id = get_name_id(bed, name);
	    if (id == -1) {
		if ( bed->l_names == bed->m_names ) {
		    bed->m_names = bed->m_names == 0 ? 2 : bed->m_names << 1;
		    bed->names = (char**)realloc(bed->names, bed->m_names*sizeof(char*));
		}
		id = bed->l_names;
		bed->names[bed->l_names++] = strdup(name);
	    }
	    int ret;	    
	    struct chrom *chrom = (struct chrom *)malloc(sizeof(struct chrom));
	    chrom->cached = chrom->max = 0;
	    chrom->a = 0;
	    chrom->id = id;
	    k = kh_put(reg, hash, bed->names[id], &ret);
	    kh_val(hash, k) = chrom;
	}
	struct bed_chrom *chrom = kh_val(hash, k);
	if ( start == end && based_1 == 0)
	    warnings("%s : line %d looks like a 1-based region. Please make sure you use right parameters.");
	    
	if ( end < beg) { int _pos = beg; beg = end; end = _pos; }
	if ( chrom->cached == chrom->max ) {
	    chrom->max = chrom->max == 0 ? 10 : chrom->max << 2;
	    chrom->a = (uint64_t*)realloc(chrom->a, chrom->max *sizeof(uint64_t));
	}
	chrom->a[chrom->cached++] = (uint64_t)start<<32 | (uint32_t)end;
	bed->region_ori ++;
	bed->length_ori += end - start;
	bed->region ++;
	bed->length += end - start;
	if ( bed->region == bed->block_size) {
	    bed_cache_update(bed);
	}
	
      clean_splits:
	free(splits);
	string.l = 0;
    }
    if ( string.m ) free(string.s);
}
// for read only
struct bed_chrom *get_chrom(struct bedaux *bed, const char *name)
{
    if (bed->flag & bed_bit_empty) {
	warnings("Try to get chrom structure from empty bed.");
	return NULL;
    }
    khint_t k;
    reghash_type *hash = (reghash_type*)file->hash;
    k = kh_get(reg, hash, name);
    if (k == kh_end(hash)) return NULL;
    struct chrom *chrom = kh_val(hash, k);
    return chrom;
}
void bed_read(struct bedaux *bed, const char *fname)
{
}
struct bedaux *bed_fork(struct bed_chrom *chrom)
{
}
struct bedaux *bed_dup(struct bedaux *bed)
{
}
int bed_getline_chrom(struct bed_chrom *chrom, struct bed_line *line)
{
}
int bed_getline(struct bedaux *bed, struct bed_line *line)
{
}
void bed_sort(struct bedaux *bed)
{
}
void bed_merge(struct bedaux *bed)
{
}
void bed_save(struct bedaux *bed, const char *fname)
{
#ifdef _DEBUG_MODE
    debug_print("[%s]", __func__);
#endif
    if ( bed->flag & bed_bit_empty ) return;
    if ( bed->flag & bed_bit_cached ) bed_fill(bed);
    
    FILE *fp = fopen(fname, "w");
    khiter_t k;
    int i, j;
    reghash_type * hash = (reghash_type*)bed->hash;
    for (i = 0; i < bed->l_names; ++i) {
	k = kh_get(reg, hash, bed->names[i]);
	if ( k != kh_end(hash) ) {
	    struct chrom * chrom = kh_val(hash, k);
	    if ( chrom == NULL)
		continue;
	    for (j = 0; j < chrom->cached; ++j)
		fprintf(fp, "%s\t%u\t%u\n", bed->names[i], bed->a[j] >> 32, (uint32_t)bed->a[j]);
	}
    }
    fclose(fp);
}
