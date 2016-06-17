/* calculate the fragements covered on the genome and the distribution 
   INPUT: bam file
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

#define check_map_same_strand(a) ((((a) & BAM_FREVERSE) << 1) ^((a) & BAM_FMREVERSE))

enum ori_type {
    left,
    right,
    unknown
};
/* output */
char *out = NULL;

/* check the bam file is coordinated or not, this func was adapt from  bam_mate.c */
/* int check_sam_is_sorted(bam_hdr_t *header) */
/* { */
/*     if (header == NULL) */
/*         return 1; */
/*     // Accept coordinate sorted. */
/*     fprintf(stderr, "[debug] %s\n", header->text); */
/*     if ((header->l_text > 3) && (strncmp(header->text, "@HD", 3) == 0)) { */
/*         char *p; */
/*         p = strstr(header->text, "\tSO:coordinate"); */
/* 	if (p == 0) */
/*             return 1; */
/*     } */
/*     return 0; */
/* } */

int cal_fra_core(samFile *in, char *out)
{
    bam_hdr_t *header = NULL;
    bam1_t *b = NULL;
    header = sam_hdr_read(in);
    int ret = 0;
    if (header == NULL) {
	fprintf(stderr, "[%s] ERROR: couldn't read header\n", __func__);
	return 1;
    }
    /* if ( check_sam_is_sorted(header) == -1) { */
    /* 	fprintf(stderr, "[%s] ERROR: coordinate is not sorted. use `samtools sort` first.\n", __func__); */
    /* 	bam_hdr_destroy(header);	 */
    /* 	return 1; */
    /* } */
    
    FILE *fp = NULL;
    int export_log = 1;
    if (out == NULL || !strcmp(out,"-")) {
	fp = stdout;
	export_log = 0;
    } else {
	fp = fopen(out, "w");
    }
    if (fp == NULL) {
	fprintf(stderr, "%s : %s\n", out, strerror(errno));
	return 1;
    }

    //kstring_t str = {0, 0, 0};
    b = bam_init1();
    int flag = ~(BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREVERSE | BAM_FMREVERSE | BAM_FDUP | BAM_FREAD1);
    int fmask = ~(1<<12);
    fprintf(fp, "#chr\tstart\tstop\tlength\tstrand\tseq_orientation\tdistribution\tname\n");
    while (sam_read1(in, header, b)>=0) {

	/* only count paired read, and just stat read1 for simple */
	if ( (b->core.flag & fmask) & flag)
	    continue;

	/* check read1 and read2 map on teh same strand */
	if ( ! check_map_same_strand(b->core.flag))
	    continue;
	
	/* check read1 and read2 map on the same chromosome */
	if ( b->core.tid != b->core.mtid)
	    continue;

	int32_t isize = b->core.isize;
	
	if(isize == 0) {
	    fprintf(stderr,"ERROR: [%s] unrecongnised map type : %s\n", __func__, bam_get_qname(b));
	    ret = 1;
	    goto clean_end;
	}
	

	int32_t start = b->core.pos;
	int32_t stop = 0;
	enum ori_type ori = unknown;
	if (isize > 0) {
	    stop = start + isize;
	    ori = left;
	} else {
	    stop = start;
	    start = start - isize;
	    ori = right;
	}
	fprintf(fp, "%s\t%d\t%d\t%d\t%c\t%s\t%s\t%s\n", header->target_name[b->core.tid], start, stop, isize > 0 ? isize : -isize, b->core.flag & BAM_FREVERSE ? '-' : '+', ori == left ? "left" : "right", b->core.flag & BAM_FPROPER_PAIR ? "properly" : "unproperly", bam_get_qname(b));
    }
  clean_end:
    
    bam_destroy1(b);
    bam_hdr_destroy(header);
    return ret;
}

int usage()
{
    
    return 1;
}
int main(int argc, char **argv)
{
    int c;
    char *out = NULL;
    while ((c=getopt(argc, argv, "o:h")) != -1) {
	switch (c) {
	    case 'o':
		out = strdup(optarg);
		break;
	    case 'h':
	    case '?':
	    default:
		return usage();
	}
    }
    char *fn_in = argc == optind ? "-" : argv[optind];
    htsFormat in_fmt;
    samFile *in = NULL;
    struct hFILE *hfile = hopen(fn_in, "r");
    if (hfile == NULL) {
	fprintf(stderr, "failed to open %s : %s\n", fn_in, strerror(errno));
	return -1;
    }
    
    if ( hts_detect_format(hfile, &in_fmt) || in_fmt.format == unknown_format) {
	fprintf(stderr, "failed to detect format : %s\n", fn_in);
	return -1;
    }

    if (hclose(hfile))
	return -1;
    if (in_fmt.category != sequence_data) {
	fprintf(stderr, "%s is not sequence data\n", fn_in);
	return -1;
    }

    if ((in = hts_open_format(fn_in, "r", &in_fmt)) == NULL) {
	fprintf(stderr, "failed to open %s : %s\n", fn_in, strerror(errno));
	return -1;
    }
    int ret;
    ret = cal_fra_core(in, out);
    hts_close(in);
    
    return ret;
}
