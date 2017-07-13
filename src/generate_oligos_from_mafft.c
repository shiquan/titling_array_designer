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
    int misc, xtra;
    int skip;
    int length;
    int depth;
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
    .skip = 35,
    .forward_only = 0,
    .depth = 1,
    .misc = 35,
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
    int fai_length;
    int n, m;
    kstring_t *str;
    uint64_t last_offset;
    uint64_t offset;
} multi_aligns_t;

multi_aligns_t *multi_aligns_cache_init() {
    int n = args.fai->n;
    multi_aligns_t *aligns = (multi_aligns_t*)malloc(sizeof(multi_aligns_t));
    aligns->length = 0;
    aligns->m = n;
    aligns->n = 0;
    aligns->str = (kstring_t*)calloc(n, sizeof(kstring_t));
    aligns->offset =0;
    return aligns;
}
int check_seqs_length(const faidx_t *idx, uint64_t *fai_length)
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
    *fai_length = length;
    return 0;
}
int fill_aligns_cache(multi_aligns_t *mcache) {
    int i;
    const faidx_t *idx = args.fai;
    mcache->length = 0;
    int depth_offset = args.length - args.length/args.depth;
    uint64_t start = mcache->offset+1 - depth_offset;
    uint64_t stop = start + args.length -1;
    for (i=0; i<mcache->m; ++i) {
	mcache->str[i].l = 0;
	const char *name = faidx_iseq(idx, i);
	int l =0;
	if (mcache->offset > mcache->fai_length) return 1;

	char *seq = faidx_fetch_seq(idx, name, start, stop, &l);
	if (l == 0 && seq == 0) return 1;
	if(mcache->length < l) mcache->length = l;
	kputsn(seq, l, &mcache->str[i]);

    }
    mcache->last_offset = mcache->offset;
    mcache->offset = stop;
    return 0;  
}
void reconstruct_seqs(multi_aligns_t *mcache)
{
    int i, j;
    mcache->n = 0;
    for (i=0; i<mcache->m; ++i) {
	kstring_t * str = &mcache->str[i];
	for (j=0; j<str->l;) {
	    unsigned char c = seq_nt4_table[str->s[j]];
	    if ( c== 4 && str->l > 0) {
		memmove(str->s +j, str->s + j + 1, str->l -j -1);
		--str->l;
	    } else {
		str->s[j] = c;
		++j;
	    }
	}
	// skip it if sequence is too short
	if (str->l) mcache->n ++;
	if (str->l < args.skip) str->l = 0;	
	str->s[str->l] = 0;

    }
}
void calculate_sw(multi_aligns_t *mcache)
{
    //args.misc = args.length * args.sa/2;
    int8_t mat[25];
    int i, j, k;
    for (i=k=0; i <4; ++i) {
	for(j=0; j<4; ++j)
	    mat[k++] = i==j ? args.sa : -args.sb;
	mat[k++] = 0;
    }
    for (j=0; j<5; ++j) mat[k++] = 0;
    kswr_t r;

    
    for (i=0; i<mcache->m-1; ++i) {
	int count = 0;
	kswq_t *q = 0;
	kstring_t *ref = &mcache->str[i];
	
	if (ref->l == 0) continue;
	for (j=i+1; j<mcache->m; ++j) {
	    kstring_t *query = &mcache->str[j];
	    if (query->l ==0) continue;
	    r = ksw_align(ref->l, (uint8_t*)ref->s, query->l, (uint8_t*)query->s, 5, mat, args.gapo, args.gape, args.xtra, &q);
	    int misc = ref->l > query->l ? ref->l *args.sa : query->l *args.sa;
	    misc = args.misc > misc ?  misc - misc/args.sa*0.1*args.sb : args.misc;
	    if (r.score > misc) {
		if (ref->l > query->l) {
		    query->l = 0;
		} else {
		    ref->l = 0;
		    break;
		}
	    	//mcache->str[j].l= 0;
	    }

	    //debug_print("score : %d, misc : %d, misc : %d, count : %d", r.score, args.misc, misc, count);
	    count ++;
	}
	free(q);
    }

}
int generate_common_oligos(multi_aligns_t *mcache)
{
    do {	
	if ( fill_aligns_cache(mcache)) break;
	reconstruct_seqs(mcache);
	if (mcache->n == 0) return 1; // end
	calculate_sw(mcache);
	int i, j;
	//debug_print("offset : %llu", mcache->offset);
	for (i=0; i<mcache->m; ++i) {

	    if (mcache->str[i].l ) {
		fprintf(stdout, "%s\t%llu\t%llu\t",args.fai->name[i], mcache->last_offset, mcache->offset);
		for (j=0; j<mcache->str[i].l; ++j)
		    putc("ACGTN"[(uint8_t)mcache->str[i].s[j]], stdout);
		putc('\n', stdout);
	    }	    
	}
    } while(1);    
    return 0;
}
int usage()
{
    fprintf(stderr,"-a  match score");
    fprintf(stderr,"-b  minmatch score");
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
    
    while ((c = getopt(argc, argv, "a:b:q:r:t:o:h?l:s:d:")) >= 0) {
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
		args.misc = atoi(optarg);
		break;
		
	    case 'o':		
		args.output = strdup(optarg);
		break;
		
	    case 'l':
		args.length = atoi(optarg);
		break;
	    case 'd':
		args.depth = atoi(optarg);
		break;
	    case 's':
		args.skip = atoi(optarg);
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
    

    uint64_t length = 0;
    if ( check_seqs_length(args.fai, &length) == 1) {
	release_args();
	error("Different sequence length ! use mafft align %s first.", args.input);
    }

    multi_aligns_t *mcache = multi_aligns_cache_init();
    mcache->fai_length = length;
    //debug_print("Start generate oligos ...");
    int ret = generate_common_oligos(mcache);
    //debug_print("Start clean memory ...");
    destroy_mcache(mcache);
    release_args();
    
    return ret;
}
