//
//
//
// Here, I defined probe is a set of oligonucletides. One oligo is a DNA sequence.
#include <stdio.h>
#include <stdlib.h>
//#include <sys/type.h>
//#include <sys/stat.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/faidx.h"
#include "htslib/kseq.h"
#include "utils.h"
#include "bed_utils.h"

#ifndef KSTRING_INIT
#define KSTRING_INIT { 0, 0, 0 }
#endif
struct args {
    // species reference genome, retrieve oligos from this reference
    const char *fasta_fname;
    // user upload input bed file
    const char *input_bed_fname;
    // serious uniq regions for design oligos
    const char *uniq_bed_fname;
    // for different projects, may tolerant some repeats
    // const char *tolerant_bed_fname;
    // output directary for keep designs and summary file
    const char *output_dir;
    const char *common_variants_fname;
    // design for each region as much as possible, 0 for default
    int must_design; 
    // defaulf oligo length is 50 now, if set to 0 dynamic mode will enabled
    // 0 for dynamic design (50b - 90b)
    int oligo_length;
    // pesudo sequences to fill short oligos, because machine can only synthesis equal length oligos
    //const char *fork_sequence;
    // origal target bed file auxiliary
    struct bedaux *target_regions;
    // intersection of input bed and uniq bed
    struct bedaux *design_regions;
    // prediction captured regions
    struct bedaux *predict_regions;
    // required a uniq database, data_required == 0 if uniq_bed_fname and tolerant_bed_fname are empty
    int data_required;
    int variants_skip_required;
    //tbx_t *uniq_data_tbx;
    int gap_size;
    // oligo coverage 
    int depth;
    // Memory cache for output probes
    kstring_t string;
    // last chromosome id, default is -1
    int last_chrom_id;
    // last region start
    int last_start;
    // last region end
    int last_end;
    // last region designed or not
    int last_is_empty;
    // current line cache
    struct bed_line line;
    faidx_t *fai;
};

struct args args = {
    .data_required = 0,
    .variants_skip_required = 0,
    .fasta_fname = 0,
    .input_bed_fname = 0,
    .uniq_bed_fname = 0,
    .common_variants_fname = 0,
    .output_dir = 0,
    .oligo_length = 50,
    .target_regions = 0,
    .design_regions = 0,
    //.uniq_data_tbx = 0,
    .gap_size = 200,
    .must_design = 0,
    .depth = 2,
    .last_start = 0,
    .last_end = 0,
    .last_is_empty = 0,
    .string = KSTRING_INIT,
    .line = BED_LINE_INIT,
    .fai = 0,
};

#define OLIGO_LENGTH_MIN 50
#define OLIGO_LENGTH_MAX 90

#define BUBBLE_GAP_MAX 30
#define BUBBLE_GAP_MIN 5

#define SMALL_REGION 200
// in case two regions come very close
const int flank_region_length = 50;
const int trim_region_length = -50;
// in case small gaps in the uniq regions
const int flank_uniq_length = 10;
const int trim_uniq_length = -10;

int usage()
{
    fprintf(stderr,
	    "generate_probes - generate DNA titling array probes from target regions and reference genome sequences.\n"
	    "Usage: \n"
	    "generate_probes [options] -fasta hg19.fa -target target.bed -uniq_regions database.bed.gz -outdir output_dir\n"
	    "options:\n"
	    "       -r, -fasta   \n"
	    "       -t, -target  \n"
	);
    return 1;
}
// quiet mode, 0 for default, will export logs
static int quiet_mode = 0;

int prase_args(int argc, char **argv)
{
    int i;
    kstring_t buff = KSTRING_INIT;
    for (i = 0; i < argc; ++i) {
	if ( i ) kputc(' ', &buff);
	kputs(argv[i], &buff);
    }
    const char *length = 0;
    const char *depth = 0;
    for (i = 0; i < argc; ) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 || strcmp(a, "-help") == 0 )
	    return usage();
	if ( strcmp(a, "-quiet") == 0 ) {
	    quiet_mode = 1;
	    continue;
	}
	const char **var = 0;
	if ( (strcmp(a, "-r") == 0 || strcmp(a, "-fasta") == 0) && args.fasta_fname == 0 )
	    var = &args.fasta_fname;
	else if ( (strcmp(a, "-t") == 0 || strcmp(a, "-target") == 0) && args.input_bed_fname == 0 )
	    var = &args.input_bed_fname;
	else if ( (strcmp(a, "-u") == 0 || strcmp(a, "-uniq_regions") == 0) && args.uniq_bed_fname == 0 )
	    var = &args.uniq_bed_fname;
	/* else if ( (strcmp(a, "-d2") == 0 || strcmp(a, "-tolerant_regions") == 0) && args.tolerant_bed_fname == 0) */
	/*     var = &args.tolerant_bed_fname; */
	else if ( (strcmp(a, "-o") == 0 || strcmp(a, "-outdir") == 0) && args.output_dir == 0 )
	    var = &args.output_dir;
	else if ( (strcmp(a, "-l") == 0 || strcmp(a, "-length") == 0) && length == 0)
	    var = &length;
	else if ( (strcmp(a, "-d") == 0 || strcmp(a, "-depth") == 0) && depth == 0 )
	    var = &depth;
	
	if ( var != 0 ) {
	    if (i == argc) {
		error_print("Missing an argument after %s.", a);
		return -2;
	    }
	    *var = argv[i++];
	    continue;
	}
	if ( strcmp(a, "-must_design") == 0) {
	    args.must_design = 1;
	}
	error_print("Unknown parameter : %s. Use -h to for more help.", a);
	return 1;
    }

    if (quiet_mode == 0) {
	LOG_print("The program was compiled at %s %s by %s.", __DATE__, __TIME__, getenv("USER"));
	LOG_print("Args: %s", buff.s);
    }
    if (args.output_dir != 0) {
	struct stat s = { 0 };
	if ( stat(args.output_dir, &s) == -1 ) {
	    if (mkdir(args.output_dir, 0755) ) {
		error_print("Failed to create directory %s.", args.output_dir);
		args.output_dir = "./";
	    } else {
		if (quiet_mode == 0)
		    LOG_print("Directory %s created.", args.output_dir);
	    }
	}
    }
    if (args.fasta_fname == 0)
	error("Required a reference genome sequence. Use -r/-fasta to specify it.");

    if (args.input_bed_fname == 0)
	error("Required a target bed file. Use -t/-target to specify it.");

    if (args.uniq_bed_fname == 0) {
	if (quiet_mode == 0)
	    LOG_print("No uniq regions datasets specified.");
    } else {
	args.data_required = 1;
	//args.uniq_data_tbx = tbx_index_load(args.uniq_bed_fname);	
	//if ( args.uniq_data_tbx == 0) {
	//    error_print("Failed to load tabix index of %s.", args.uniq_bed_fname);
	//    return 1;
	//}
    }

    if (args.common_variants_fname == 0) {
	if (quiet_mode == 0)
	    LOG_print("No common variantions datasets specified.");
    } else {
	args.variants_skip_required = 1;
    }
    if (length != 0) {
	args.oligo_length = atoi(length);
	if (args.oligo_length < 40 && args.oligo_length > 0) {
	    if (quiet_mode == 0)
		LOG_print("Oligo length is too short, <40. Force set to %db.", OLIGO_LENGTH_MIN);
	    args.oligo_length = OLIGO_LENGTH_MIN;
	} else if (args.oligo_length > 100) {
	    if (quiet_mode == 0)
		LOG_print("Oligo length is too long, > 100. Force set to %db.", OLIGO_LENGTH_MAX);
	    args.oligo_length = OLIGO_LENGTH_MAX;
	}
    }
    if (depth != 0) {
	args.depth = atoi(depth);
	if (args.depth == 0) {
	    args.depth = 2;
	    if (quiet_mode == 0) LOG_print("Depth should greater than 0. Force set to 2.");
	} else if (args.depth > 10) {
	    args.depth = 4;
	    if (quiet_mode == 0) LOG_print("Depth should not greater than 10. Force set to 4. Usually 2 ~ 4x.");
	}	
    }
    if (args.oligo_length == 0 && quiet_mode == 0) {
	LOG_print("Use dynamic design mode, the length of oligos will set from 50b to 90b.");
    }

    args.target_regions = bedaux_init();

    // assume input is 0 based bed file.
    set_based_0();
    bed_read(args.target_regions, args.input_bed_fname);

    // bed will merge auto.
    // for some target regions, usually short than 100, expand to 100 size.
    // the reason we define the expand size to 100 is our oligo length usually smaller than 100. so
    // for a single nucletide variantion, there should be at least two different oligos cover it at any depth.
    bed_round(args.target_regions, 100);

    // merge two near regions
    struct bedaux *bed = bed_dup(args.target_regions);

    // the target regions usually generated from raw exome or user specified regions. some of them are very close,
    // if two nearby regions are close enough, we could just treat them as one continuous region, and design tiled oligos.
    bed_flktrim(bed, flank_region_length, flank_region_length);
    bed_flktrim(bed, trim_region_length, trim_region_length);

    if ( args.data_required == 1) {
	htsFile *fp = hts_open(args.uniq_bed_fname, "r");
	tbx_t *tbx = tbx_index_load(args.uniq_bed_fname);
	// this function will find the overlap regions of target and uniq dataset for design. And more, for exactly
	// non-overlaped regions, means hang regions without any overlap with uniq database, will find the most nearest
	// uniq regions for design if possible
	args.design_regions = bed_find_rough_bigfile(bed, fp, tbx, args.gap_size, OLIGO_LENGTH_MAX);
	hts_close(fp);
	tbx_destroy(tbx);
    } else {
	args.design_regions = bed_dup(bed);
    }

    // sometimes, uniq regions in database are very small and will break a contine region into several small regions.
    // merge these small regions into one piece if the gap between them is shorter than flank_uniq_length*2.
    // defined the flank_uniq_length based on the insert gap size of bubble oligos.
    bed_flktrim(args.design_regions, flank_uniq_length, flank_uniq_length);
    bed_flktrim(args.design_regions, trim_uniq_length, trim_uniq_length);
    bed_destroy(bed);
    return 0;
}
float calculate_GC(const char *seq, int length)
{
    int i, j;
    for (i = 0, j = 0; i < length; ++i)
	if (seq[i] == 'G' || seq[i] == 'C' || seq[i] == 'g' || seq[i] == 'c') j++;
    return (float)j/length;
}
// for much design regions, usually very short, try to use short oligos for better oligos
void must_design(int cid, int start, int end)
{
    if (args.must_design == 0)
	return;
    // expand the small regions into longer one, the size of new region should consider of length of oligo and depth.
    // the algrithm here to generate oligos based on depth is by set oligo start from the 1/n part of previous oligos
    titling_design(cid, start, end);
}
// bubble design is one oligo cover two regions within a tolerant gap. There will be some fork sequence in the gap to make
// sure the bubble structure is constructed. but here only use blocks to remember to stand for bubbles. fork sequences will
// not be added in this program.
// for all bubble design, only consider to cover the last region and current region.
// And also should put the gaps in the middle of oligos as much as possible for better performance.
// based on different situtions of last region, if last region is designed and the purpose of bubble oligo is to cover
// current region, the main body of bubble oligos should located in the current region, but of course is current region
// is very short, we should not expand it for better performace, because the last region is very close to current one, the
// oligo of last region will capture big enough fragement to enhance the coverage of current one. so it is OK even there are
// not any designed oligos located in current region.
void bubble_design()
{
}
// rough design, not consider of common variants
void titling_design(int cid, int start, int end)
{
    int length = end - start;
    int oligo_length = args.oligo_length == 0 ? length < SMALL_REGION ? OLIGO_LENGTH_MIN : OLIGO_LENGTH_MAX : args.oligo_length;
    int n_parts = length/oligo_length == 0 ? args.depth : length/oligo_length * args.depth;
    int part = length/n_parts;
    int offset =  part > oligo_length ? 0 : (oligo_length - part)/2;
    float mid = (float)n_parts/2;

    int i;
    for (i = 0; i < n_parts; ) {
	int rank = 1;
	int start_pos = start + i *part;
	start_pos = start_pos - offset;
	if (mid > i && start_pos < start) {
	    start_pos = start;
	    if (start_pos + oligo_length > end) rank = 0;
	}
	if (mid > i && start_pos + oligo_length > end) {
	    start_pos = end - oligo_length;
	    if (start_pos < start) rank = 0; 
	}
	int l = 0;
	char *seq = faidx_fetch_seq(args.fai, args.design_regions->names[cid], start_pos, start_pos+oligo_length-1, &l);
	if (l < oligo_length) {
	    if (seq) free(seq);
	    continue;
	}
	ksprintf(&args.string, "%s\t%d\t%d\t%d\t%s\t%d\t%d,\t%d,\t%.2f\t%d\n", args.design_regions->names[cid], start_pos, start_pos + oligo_length, oligo_length, seq, 1, start_pos, start_pos+oligo_length, calculate_GC(seq, oligo_length), rank);
	
    }
    
    
}
// format of oligos file.
// chr, start(0-based), end, seq_length, sequences, n_blocks, blocks(seperated by commas, sometime the sequences are consist of different parts from reference sequences), gc percent, type, rank, score
int generate_oligos_core()
{
    struct bed_line *line = &args.line;
    // first line
    if ( args.last_chrom_id == -1 ) 
	goto design;

    // check last region is empty
    // first line of next chromosome
    if ( args.last_chrom_id != line->chrom_id ) {
	if ( args.last_is_empty == 1) {
	    must_design( args.last_chrom_id, args.last_start, args.last_end);
	}
	goto print_line;
    }
    int gap = line->start - args.last_end;
    if ( args.last_is_empty == 1) {
	if ( gap > BUBBLE_GAP_MAX ) {
	    must_design( args.last_chrom_id, args.last_start, args.last_end);
	} else {
	    // if gap size is short, design bubble oligos, create two blocks.
	    // remember, because we have expand and merge nearby regions after generate uniq design regions, so there should 
	    // not be more than two blocks in the downstream design.
	    bubble_design();
	}
    }
  design:	
    if ( bed_getline(args.design_regions, line) )
	return 1;
    
    int length = line->end - line->start;
    if (length < args.oligo_length) {
	if ( gap > BUBBLE_GAP_MAX ) {
	    args.last_is_empty = 1;
	} else {
	    bubble_design();
	}
    } else {
	titling_design(line->chrom_id, line->start, line->end);
    }
    
  print_line:	
    args.last_chrom_id = line->chrom_id;
    args.last_start = line->start;
    args.last_end = line->end;

    return 0;
}
void generate_oligos()
{
    args.fai = fai_load(args.fasta_fname);
    if (args.fai == NULL ) {
	if (fai_build(args.fasta_fname) == -1)
	    error("Failed to build the index of %s.", args.fasta_fname);
	args.fai = fai_load(args.fasta_fname);
    }
    struct bed_line line = BED_LINE_INIT;
    int oligo_length = args.oligo_length == 0 ? OLIGO_LENGTH_MIN : args.oligo_length;

    // create probe file in the out directary
    kstring_t probe_path = KSTRING_INIT;
    if ( args.output_dir)
	kputs(args.output_dir, &probe_path);
    if (probe_path.l && probe_path.s[probe_path.l-1] != '/')
	kputc('/', &probe_path);
    kputs("probes.txt.gz", &probe_path);    
    BGZF *fp = bgzf_open(probe_path.s, "w");
    free(probe_path.s);    
    if (fp == NULL)
	error("Failed to write %s : %s.", probe_path.s, strerror(errno));

    // write header to probe file
    kstring_t header = KSTRING_INIT;
    kputs("##filetype=probe\n", &header);    
    ksprintf(&header,"##length=%d\n", oligo_length);
    kputs("#chrom\tstart\tend\tseq_length\tsequence\tn_block\tstarts\tends\tGC_content\trank\n", &header);
    bgzf_write(fp, &header, header.l);
    free(header.s);
    
    while (1) {	
	if ( generate_oligos_core() ) break;	
	if ( args.string.l ) {
	    bgzf_write(fp , args.string.s, args.string.l);
	    args.string.l = 0;
	}	
    }    

    bgzf_close(fp);
}
void export_summary_reports()
{
    kstring_t path = KSTRING_INIT;
    if ( args.output_dir )
	kputs(args.output_dir, &path);
    if (path.l && path.s[path.l-1] != '/')
	kputc('/', &path);

    int l = path.l;
    kputs("target_regions.bed", &path);
    bed_save(args.target_regions, path.s);

    path.l = l;
    kputs("design_regions.bed", &path);
    bed_save(args.design_regions, path.s);

    path.l = l;
    kputs("predict_regions.bed", &path);
    bed_save(args.predict_regions, path.s);

    free(path.s);
}
void clean_memory()
{
    bed_destroy(args.target_regions);    
    bed_destroy(args.design_regions);
    bed_destroy(args.predict_regions);
    //if (args.data_required ) tbx_destroy(args.uniq_data_tbx);
    fai_destroy(args.fai);
    free(args.string.s);
}
int main(int argc, char **argv)
{
    if (argc == 0)
	return usage();
    
    if ( prase_args(argc, argv) != 0 ) 
	return 1;
    
    generate_oligos();
    
    export_summary_reports();
    
    clean_memory();
    
    return 0;
}

