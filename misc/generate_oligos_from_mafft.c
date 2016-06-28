/*  generate oligos for most common regions and alternative locus from mafft alignment results
 *  INPUT: mafft.fa (preindexed by samtools faidx)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/bgzf.h"
#include "utils.h"
#include "ksw.h"

typedef struct {
    int sa, sb, gapo, gape, forward_only;
    int minsc, xtra;
    int length;
    char *input;
    char *output;
    faidx_t *fai;
} arg_t;

arg_t args = {
    .sa = 1,
    .sb = 3,
    .gapo = 3,
    .gape = 1,
    .length = 50,
    .forward_only = 0,
    .minsc = 0,
    .xtra = KSW_XSTART,
    .input = NULL,
    .output = NULL,
};
unsigned char seq_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void release_args()
{
    if (args.input) free(args.input);
    if (args.output) free(args.output);
    fai_destroy(args.fai);    
}
typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};


typedef struct {
    int length;
    int n, m;
    kstring_t *str;
    uint32_t *offsets;
} multi_aligns_t;

multi_aligns_t *multi_aligns_cache_init() {
    int n = args.fai->n;
    multi_aligns_t *aligns = (multi_aligns_t*)malloc(sizeof(multi_aligns_t));
    aligns->length = 0;
    aligns->m = n;
    aligns->n = 0;
    aligns->str = (kstring_t*)calloc(n, sizeof(kstring_t));
    aligns->offsets = (uint32_t*)calloc(n, sizeof(uint32_t));
    int i;
    for (i=0; i<n; ++i)
	aligns->offsets[i] = 0;
    
    return aligns;
}
int check_seqs_length(const faidx_t *idx)
{
    int i;
    int length = 0;
    int n = faidx_nseq(idx);
    for (i=0; i<n; ++i) {
	const char *name = faidx_iseq(idx, i);
	khiter_t l = kh_get(s, idx->hash, name);
	faidx1_t x = kh_value(idx->hash, l);
	if (length == 0) length = x.len;
	if (length != x.len) return 1;	
    }
    return 0;
}
int fill_aligns_cache(multi_aligns_t *mcache) {
    int i;
    const faidx_t *idx = args.fai;
    for (i=0; i<mcache->m; ++i) {
	mcache->str[i].l = 0;
	const char *name = faidx_iseq(idx, i);
	int l =0;
	uint32_t start = mcache->offsets[i];
	uint32_t stop = start + mcache->length;
	char *seq = faidx_fetch_seq(idx, name, start, stop, &l);
	if (l == 0) return 1;
	kputsn(seq, l, &mcache->str[i]);	
    }
    return 0;  
}
void reconstruct_seqs(multi_aligns_t *mcache)
{
    int i, j;
    for (i=0; i<mcache->m; ++i) {
	kstring_t * str = &mcache->str[i];
	for (j=0; j<str->l;) {
	    unsigned char c = seq_nt4_table[str->s[j]];
	    if ( c== 4) {
		memmove(str->s +j, str->s + j + 1, str->l -j -1);
		--str->l;
	    } else {
		str->s[j] = c;
		++j;
	    }
	}
	// skip it if sequence is too short
	if (str->l < mcache->length/2) str->l = 0;
	
	str->s[str->l] = 0;
	debug_print("%s , %zu", str->s, str->l);
    }
}
void calculate_sw(multi_aligns_t *mcache)
{
    
}
int generate_common_oligos(multi_aligns_t *mcache)
{
    do {	
	if ( fill_aligns_cache(mcache)) break;
	calculate_sw(mcache);
	int i;
	if (mcache->n == 0) return 1; // end
	for (i=0; i<mcache->n; ++i) {
	    printf("%s\n", mcache->str[i].s);
	}
    } while(1);    
    return 0;
}
int usage()
{
    return 1;
}
void destroy_mcache(multi_aligns_t *mcache)    
{
    assert(mcache);
    int i;
    for (i=0; i<mcache->m; ++i) 
	if (mcache->str[i].m) free(mcache->str[i].s);

    free(mcache->str);
    free(mcache);
}
int main(int argc, char **argv)
{
    int c;
    
    while ((c = getopt(argc, argv, "a:b:q:r:t:o:h?l:")) >= 0) {
	switch(c) {
	    
	    case 'a':
		args.sa = atoi(optarg);
		break;
		
	    case 'b':
		args.sb = atoi(optarg);
		break;
		
	    case 'q':
		args.gapo = atoi(optarg);
		break;
		
	    case 'r':
		args.gape = atoi(optarg);
		break;
		
	    case 't':
		args.minsc = atoi(optarg);
		break;
		
	    case 'o':		
		args.output = strdup(optarg);
		break;
		
	    case 'l':
		args.length = atoi(optarg);
		break;
		
	    case 'h':
	    case '?':		
		return usage();
		
	    default:
		error("Unknown parameter %c", c);
		
	}
    }

    if (argc == optind) {
	error("Usage: %s <seqs.mafft.fa>", argv[0]);
    } else {
	args.input = strdup(argv[optind]);
    }
    args.fai = fai_load(args.input);
    if (args.fai == NULL)
	error("Failed to load index.");

    if ( check_seqs_length(args.fai) == 1) {
	release_args();
	error("Different sequence length ! use mafft align %s first.", args.input);
    }
    multi_aligns_t *mcache = multi_aligns_cache_init();
    debug_print("Start generate oligos ...");
    int ret = generate_common_oligos(mcache);
    debug_print("Start clean memory ...");
    destroy_mcache(mcache);
    release_args();
    
    return ret;
}
