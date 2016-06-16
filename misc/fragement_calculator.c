/* calculate the fragements covered on the genome and the distribution 
   INPUT: bam file
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

/* output */
char *out = NULL;

/* check the bam file is coordinated or not */
int check_sam_is_sorted(bam_hdr_t *h)
{
    
}

int cal_fra_core(int argc, char **argv)
{
}

int usage()
{
    
    return 1;
}
int main(int argc, char **argv)
{
    int c;
    while ((c=getopt(argc, argv, "o:h")) != -1) {
	switch (c) {
	    case 'o':
		args.out = strdup(optarg);
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
    bam_hdr_t *header = NULL;
    bam1_t *line = NULL;
    hFILE *hfile = hts_hopen(fn_in, "r");
    if (hfile == NULL) {
	fprintf(stderr, "failed to open %s : %s\n", fn_in, strerror(errno));
	return -1;
    }
    
    if ( hts_detect_format(hfile, &in_fmt) ) {
	fprintf(stderr, "failed to detect format : %s\n", fn_in);
	return -1;
    }

    hclose(hfile);
    if (in_fmt.category != sequence_data) {
	fprintf(stderr, "%s is not sequence data\n", fn_in);
	return -1;
    }

    if ((in = hts_open_format(fn_in, "r", &in_fmt)) == NULL) {
	fprintf(stderr, "failed to open %s : %s\n", fn_in, strerror(errno));
	return -1;
    }
    
     
}
