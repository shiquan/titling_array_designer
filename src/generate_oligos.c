// generate_oligo.c - program to generate oligonucletides to capture target genome.
// 
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/faidx.h"
#include "htslib/kseq.h"
#include "utils.h"
#include "bed_utils.h"
#include "version.h"

#define ROUND_SIZE  100
#define DEPTH_LIMIT  20

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
    // project name
    const char *project_name;
    int debug_mode;
    // design for each region as much as possible, 0 for default
    int must_design; 
    // defaulf oligo length is 50 now, if set to 0 dynamic mode will enabled
    // 0 for dynamic design (50b - 90b)
    int oligo_length;
    // temp parameters for dynamic design mode
    int min_oligo_length;
    int max_oligo_length;
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
    float depth;
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
    kstring_t commands;
    uint32_t probes_number;
};

struct args args = {
    .data_required = 0,
    .variants_skip_required = 0,
    .fasta_fname = 0,
    .input_bed_fname = 0,
    .uniq_bed_fname = 0,
    .project_name = 0,
    .debug_mode = 0,
    .common_variants_fname = 0,
    .output_dir = 0,
    .oligo_length = 50,
    .min_oligo_length = 0,
    .max_oligo_length = 0,
    .target_regions = 0,
    .design_regions = 0,
    //.uniq_data_tbx = 0,
    .gap_size = 200,
    .must_design = 0,
    .depth = 2,
    .last_chrom_id = -1,
    .last_start = 0,
    .last_end = 0,
    .last_is_empty = 0,
    .string = KSTRING_INIT,
    .line = BED_LINE_INIT,
    .commands = KSTRING_INIT,
    .fai = 0,
    .probes_number = 0,
};

static int oligo_length_minimal = 50;
static int oligo_length_maxmal = 120;
// static int oligo_length_usual  = 90;
static char *get_version(void)
{
    return OLIGOS_VERSION;
}
void set_oligo_length_min(int length)
{
    oligo_length_minimal = length;
}
void set_oligo_length_max(int length)
{
    oligo_length_maxmal = length;
}

#define BUBBLE_GAP_MAX 30
#define BUBBLE_GAP_MIN 0

#define SMALL_REGION 200
// in case two regions come very close
static int flank_region_length = 50;
static int trim_region_length = -50;
void set_flank_trim_regions_length(int length)
{
    if ( length < 0 )
	length = -1 * length;
    flank_region_length = length;
    trim_region_length = -1 * length;    
}

// in case small gaps in the uniq regions
static int flank_uniq_length = 10;
static int trim_uniq_length = -10;
void set_flank_trim_uniq_length(int length)
{
    if ( length < 0 )
	length = -1 * length;
    flank_uniq_length = length;
    trim_uniq_length = -1 * length;
}

void titling_design(int cid, int start, int end);

int usage()
{
    fprintf(stderr,
	    "generate_probes - generate DNA titling array probes from target regions and reference genome sequences.\n"
	    "Usage: \n"
	    "generate_probes [options] -fasta hg19.fa -target target.bed -uniq_regions database.bed.gz -outdir output_dir\n"
	    "Options:\n"
	    "  -r, -fasta [fasta file]\n"
	    "            reference genome sequences, in fasta format.\n"
	    "  -t, -target [bed file]\n"
	    "            target bed file to design.\n"
	    "  -o, -outdir [dir]\n"
	    "            output directary, will create it if not exists.\n"
	    "  -u, -database [tabix-indexed bed file]\n"
	    "            database of non-repeats, or user pre-defined designable regions.\n"
	    "  -l, -length [50]\n"
	    "            pre-defined oligo length, usually from 50 base to 90 base, set 0 for dynamic design.\n"
            "  -min INT \n"
            "            minimal oligo length for dynamic design mode\n"
            "  -max INT \n"
            "            maximal oligo length for dynamic design mode\n"
	    "  -d, -depth [2]\n"
	    "            oligo depths pre base, increase this value will increase the dense of oligos.\n"
	    "  -p, -project [string]\n"
	    "            project id\n"
	    "  -must_design\n"
	    "            if no uniq regions around small target, must design it no matter repeat regions.\n"
	    "  -h, -help\n"
	    "            for help information.\n"
	    "Version: %s\n"
	    "Bugs report: shiquan@genomics.cn\n"
	    "Homepage: https://github.com/shiquan/titling_array_designer\n", OLIGOS_VERSION
	);
    return 1;
}
// quiet mode, 0 for default, will export logs
static int quiet_mode = 0;

int parse_args(int argc, char **argv)
{
    int i;

    for (i = 0; i < argc; ++i) {
	if ( i ) kputc(' ', &args.commands);
	kputs(argv[i], &args.commands);
    }
    const char *length = 0;
    const char *depth = 0;
    const char *max_oligo_length = 0;
    const char *min_oligo_length = 0;
    for (i = 0; i < argc; ) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 || strcmp(a, "-help") == 0 )
	    return usage();
	if ( strcmp(a, "-quiet") == 0 ) {
	    quiet_mode = 1;
	    continue;
	}
        if ( strcmp(a, "-debug") == 0 ) {
            args.debug_mode = 1;
            continue;
        }
	const char **var = 0;
	if ( (strcmp(a, "-r") == 0 || strcmp(a, "-fasta") == 0) && args.fasta_fname == 0 )
	    var = &args.fasta_fname;
	if ( (strcmp(a, "-p") == 0 || strcmp(a, "-project") == 0) && args.project_name == 0)
	    var = &args.project_name;
	else if ( (strcmp(a, "-t") == 0 || strcmp(a, "-target") == 0) && args.input_bed_fname == 0 )
	    var = &args.input_bed_fname;
	else if ( (strcmp(a, "-u") == 0 || strcmp(a, "-database") == 0) && args.uniq_bed_fname == 0 )
	    var = &args.uniq_bed_fname;
	else if ( (strcmp(a, "-o") == 0 || strcmp(a, "-outdir") == 0) && args.output_dir == 0 )
	    var = &args.output_dir;
	else if ( (strcmp(a, "-l") == 0 || strcmp(a, "-length") == 0) && length == 0)
	    var = &length;
	else if ( (strcmp(a, "-d") == 0 || strcmp(a, "-depth") == 0) && depth == 0 )
	    var = &depth;
        else if ( (strcmp(a, "-min") == 0) && min_oligo_length == NULL )
            var = &min_oligo_length;
        else if ( (strcmp(a, "-max") == 0 ) && max_oligo_length == NULL )
            var = &max_oligo_length;                  
	
	if ( var != 0 ) {
	    if (i == argc) {
		error_print("Miss an argument after %s.", a);
		return -2;
	    }
	    *var = argv[i++];
	    continue;
	}
	if ( strcmp(a, "-must_design") == 0) {
	    args.must_design = 1;
	    continue;
	}
	error_print("Unknown parameter : %s. Use -h to for more help.", a);
	return 1;
    }

    if (quiet_mode == 0) {
	LOG_print("Version %s.", OLIGOS_VERSION);
	LOG_print("Args: %s", args.commands.s);
    }
    if (args.project_name == 0) {
	error("No project name. Please use -p or -project to spectify.");
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
	error("Required a reference genome sequence. Use -r or -fasta to specify.");

    if (args.input_bed_fname == 0)
	error("Required a target bed file. Use -t or -target to specify.");

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
    if ( min_oligo_length ) {
        args.min_oligo_length = str2int(min_oligo_length);
        set_oligo_length_min(args.min_oligo_length);
    }

    if ( max_oligo_length ) {
        args.max_oligo_length = str2int(max_oligo_length);
        set_oligo_length_max(args.max_oligo_length);
    }
    if (length != 0) {
	args.oligo_length = atoi(length);
	if (args.oligo_length < 40 && args.oligo_length > 0) {
	    if (quiet_mode == 0)
		LOG_print("Oligo length is too short, < 40. Force set to %db.", oligo_length_minimal);
	    args.oligo_length = oligo_length_minimal;
	} else if (args.oligo_length > oligo_length_maxmal) {
	    if (quiet_mode == 0)
		LOG_print("Oligo length is too long, > 100. Force set to %db.", oligo_length_maxmal);
	    args.oligo_length = oligo_length_maxmal;
	}
    }
    if (depth != 0) {
	args.depth = atof(depth);
	if (args.depth == 0) {
	    args.depth = 2;
	    if (quiet_mode == 0) {
		LOG_print("Depth is unknown. Force set to 2. %s", depth);
	    }
	} else if (args.depth > DEPTH_LIMIT) {
	    args.depth = DEPTH_LIMIT;
	    if (quiet_mode == 0)
                LOG_print("Depth capped to %d. Force set to %d. Usually 2 ~ 4x.", DEPTH_LIMIT, DEPTH_LIMIT);
	}	
    }
    if (args.oligo_length == 0 && quiet_mode == 0) {
	LOG_print("Use dynamic design mode, the length of oligos will set from %dnt to %dnt.", oligo_length_minimal, oligo_length_maxmal);
    }

    args.target_regions = bedaux_init();

    // assume input is 0 based bed file.
    set_based_0();
    
    if ( bed_read(args.target_regions, args.input_bed_fname) )
        error("Empty file, %s", args.input_bed_fname);

    // bed will merge auto.
    // for some target regions, usually short than ROUND_SIZE, expand to ROUND_SIZE.
    // the reason we define the expand size to ROUND_SIZE is our oligo length usually smaller than ROUND_SIZE. so
    // for a single nucletide variantion, there should be at least two different oligos cover it at any depth.
    bed_round(args.target_regions, ROUND_SIZE);

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
	args.design_regions = bed_find_rough_bigfile(bed, fp, tbx, args.gap_size, oligo_length_maxmal);
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
float repeat_ratio(char *seq, int length)
{
    int i, j;
    for (i = 0, j = 0; i < length; ++i) {
	switch(seq[i]) {
	    case 'A' :
		break;
	    case 'T':
		break;
	    case 'G':
		break;
	    case 'C':
		break;
	    case 'a':
		seq[i] = 'A', j++;
		break;
	    case 't':
		seq[i] = 'T', j++;
		break;
	    case 'g':
		seq[i] = 'G', j++;
		break;
	    case 'c':
		seq[i] = 'C', j++;
		break;
	    case 'N':
	    default:
		error("There is a N is seq %s.", seq);
		break;
	}
    }
    return (float)j/length;
}
// for much design regions, usually very short, try to use short oligos for better oligos
void must_design(int cid, int start, int end)
{
    if (args.must_design == 1) {	
	// expand the small regions into longer one, the size of new region should consider of length of oligo and depth.
	// the algrithm here to generate oligos based on depth is by set oligo start from the 1/n part of previous oligos
	titling_design(cid, start, end);
    }
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
int bubble_design(int cid, int last_start, int last_end, int start, int end)
{
    int length = end - start + last_end - last_start;
    int oligo_length = args.oligo_length == 0 ?
        length < SMALL_REGION ? oligo_length_minimal : oligo_length_maxmal
        : args.oligo_length;
    float n_parts = (float)length/oligo_length < 1 ? args.depth : (float)length/oligo_length * args.depth;
    int part = length/n_parts;
    int offset = part > oligo_length ? 0 : (oligo_length - part)/2;
    int i;
    int head_length = last_end - last_start;
    int tail_length = end - start;
    // only works when head length smaller than oligo length, for longer region use titling_design() instead.
    assert(oligo_length > head_length);
    if ( args.debug_mode ) {
        debug_print("last empty: %d\tlast: %d-%d\t%s:%d-%d\t%d\tn_part: %f\tpart: %d\toffset: %d\thead length: %d\ttail length: %d",
                    args.last_is_empty, args.last_start, args.last_end,
                    args.design_regions->names[cid], start, end, oligo_length, n_parts, part, offset, head_length, tail_length);
    }

    // if length of regions shorter than oligo length, skip the tail.
    if ( head_length + tail_length  < oligo_length )
        return 1;
    kstring_t string = KSTRING_INIT;    
    for (i = 0; i < n_parts; ++i ) {
        int rank = 1;
        int offset_l = i * part;
        int start_pos = offset_l > head_length ? start + offset_l - head_length -1: last_start + offset_l-1;
        // start_pos = start_pos - offset;
        if (start_pos < last_start) {
            start_pos = last_start;
            // rank = 0;
        }
        int end_pos = start_pos >= start ? start_pos + oligo_length : start + oligo_length - (last_end-start_pos);
        if (end_pos > end) {
            end_pos = end;
            start_pos = end_pos - oligo_length >= start ? end_pos - oligo_length : last_end - (oligo_length - (end_pos - start));
        }
        int l = 0;
        //debug_print("%d\t%d\t%d\t%d\t%d\t%d\n", start_pos, end_pos, last_start, last_end, start, end);
        if ( start_pos < start) {
            char *head = faidx_fetch_seq(args.fai, args.design_regions->names[cid], start_pos+1, last_end, &l);
            char *tail = faidx_fetch_seq(args.fai, args.design_regions->names[cid], start+1, end_pos, &l);            
            kputs(head, &string);                
            kputs(tail, &string);
            free(head);
            free(tail);
            float repeat = repeat_ratio(string.s, string.l);
            float gc = calculate_GC(string.s, string.l);
            ksprintf(&args.string, "%s\t%d\t%d\t%d\t%s\t%d\t%d,%d,\t%d,%d,\t%.2f\t%.2f\t%d\n", args.design_regions->names[cid], start_pos, end_pos, oligo_length, string.s, 2, start_pos, start, last_end, end_pos, repeat, gc, rank);
            args.probes_number ++;
        
        } else {
            char *seq = faidx_fetch_seq(args.fai, args.design_regions->names[cid], start_pos+1, end_pos, &l);
            kputs(seq, &string);
            free(seq);
            float repeat = repeat_ratio(string.s, string.l);
            float gc = calculate_GC(string.s, string.l);
            ksprintf(&args.string, "%s\t%d\t%d\t%d\t%s\t%d\t%d,\t%d,\t%.2f\t%.2f\t%d\n", args.design_regions->names[cid], start_pos, end_pos, oligo_length, string.s, 1, start_pos, end_pos, repeat, gc, rank);
            args.probes_number ++;
        }
        if(string.l != oligo_length) {
            fprintf(stderr, "%s\t%d\t%d\t%d\t%s\t%d\t%d,%d,\t%d,%d,\n", args.design_regions->names[cid], start_pos, end_pos, oligo_length, string.s, 1, start_pos, start, last_end, end_pos);
            exit(1);
        }
        string.l = 0;
    }
    free(string.s);
    return 0;
}
// rough design, not consider of common variants
void titling_design(int cid, int start, int end)
{
    int length = end - start;
    int oligo_length = args.oligo_length == 0 ?
        length < SMALL_REGION ? oligo_length_minimal : oligo_length_maxmal
        : args.oligo_length;
    float n_parts = (float)length/oligo_length < 1 ? args.depth : (float)length/oligo_length * args.depth;
    int part = length/n_parts;
    int offset =  part > oligo_length ? 0 : (oligo_length - part)/2;
    float mid = (float)n_parts/2;
    if ( args.debug_mode ) {
        debug_print("last empty:%d\tlast:%d-%d\t%s:%d-%d\t%d\tn_parts: %f\tpart: %d\toffset: %d",
                    args.last_is_empty, args.last_start, args.last_end,args.design_regions->names[cid], start, end, oligo_length, n_parts, part, offset);
    }
    int i;
    for (i = 0; i < n_parts; ++i) {
	int rank = 1;
	int start_pos = start + i *part;
	start_pos = start_pos - offset;
	if (mid > i && start_pos < start) {
	    start_pos = start;
	    if (start_pos + oligo_length > end) rank = 0;
	}
	if (mid < i && start_pos + oligo_length > end) {
	    start_pos = end - oligo_length;
	    if (start_pos < start) rank = 0; 
	}
	int l = 0;
	char *seq = faidx_fetch_seq(args.fai, args.design_regions->names[cid], start_pos+1, start_pos+oligo_length, &l);
	if (l < oligo_length) {
	    if (seq) free(seq);
	    continue;
	}
	float repeat = repeat_ratio(seq, oligo_length);
	float gc = calculate_GC(seq, oligo_length);
	ksprintf(&args.string, "%s\t%d\t%d\t%d\t%s\t%d\t%d,\t%d,\t%.2f\t%.2f\t%d\n", args.design_regions->names[cid], start_pos, start_pos + oligo_length, oligo_length, seq, 1, start_pos, start_pos+oligo_length, repeat, gc, rank);
	args.probes_number ++;
        free(seq);
    }
}
// format of oligos file.
// chr, start(0-based), end, seq_length, sequences, n_blocks, blocks(seperated by commas, sometime the sequences are consist of different parts from reference sequences), gc percent, type, rank, score
int generate_oligos_core()
{
    struct bed_line *line = &args.line;
    if ( bed_getline(args.design_regions, line) ) {
        if (args.last_is_empty == 1) {
            must_design( args.last_chrom_id, args.last_start, args.last_end);
        }
	return 1;
    }
    // if databases is not merged properly    
    if ( args.last_chrom_id == line->chrom_id )  {
        // totally overlapped
        if ( args.last_end > line->end)
            return 0;
        // trim new region
        if ( args.last_end > line->start ) {
            line->start = args.last_end;
        }
    }

    int length = line->end - line->start;        
    // first line
    if ( args.last_chrom_id == -1 )
	goto design;
    
    // check last region is empty
    // first line of next chromosome
    if ( args.last_chrom_id != line->chrom_id ) {
	if ( args.last_is_empty == 1) {
	    must_design(args.last_chrom_id, args.last_start, args.last_end);
	}
        args.last_is_empty = 0;
	goto design;
    }
    
    if ( args.last_is_empty == 1) {
        int gap = line->start - args.last_end;
        if ( gap > BUBBLE_GAP_MAX ) {
	    must_design(args.last_chrom_id, args.last_start, args.last_end);
	} else {
            if (args.oligo_length && args.last_end - args.last_start > args.oligo_length) {
                error("%d\t%d\t%d\t%d", args.last_end, args.last_start, args.oligo_length, length);            
            }
	    // if gap size is short, design bubble oligos, create two blocks.
	    // remember, because we have expand and merge nearby regions after generate uniq design regions, so there should 
	    // not be more than two blocks in the downstream design.
	    if ( bubble_design(args.last_chrom_id, args.last_start, args.last_end, line->start, line->end) )
                return 0;
	}
        args.last_is_empty = 0;
    }

    if ( length > args.oligo_length) {
        goto design;
    } 
    
  design:
    assert(args.last_is_empty == 0);
    if (length < args.oligo_length || (args.oligo_length == 0 && length < oligo_length_minimal)) {
        args.last_is_empty = 1;
        goto print_line;
    } else {
	titling_design(line->chrom_id, line->start, line->end);
        args.last_is_empty = 0;
        goto print_line;
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
    // struct bed_line line = BED_LINE_INIT;
    
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

    int oligo_length = args.oligo_length == 0 ? oligo_length_maxmal : args.oligo_length;
    // write header to probe file
    kstring_t header = KSTRING_INIT;
    kputs("##filetype=probe\n", &header);
    ksprintf(&header, "##generate_oligos Version=%s\n", get_version());
    ksprintf(&header, "##project_name=%s\n", args.project_name);
    ksprintf(&header, "##max_length=%d\n", oligo_length);
    // ksprintf(&header, "##oligo_number=%u\n", args.probes_number); // should always be 0
    ksprintf(&header, "##Command=%s\n", args.commands.s);    
    kputs("#chrom\tstart\tend\tseq_length\tsequence\tn_block\tstarts\tends\trepeat_ratio\tGC_content\trank\n", &header);
    if ( bgzf_write(fp, header.s, header.l) != header.l )
        error ( "Write error : %d.", fp->errcode);
    free(header.s);
    
    while (1) {	
	if ( generate_oligos_core() ) break;	
	if ( args.string.l ) {
	    if ( bgzf_write(fp , args.string.s, args.string.l) != args.string.l)
                error("Writer error : %d.", fp->errcode);
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
void clean_memory(void)
{
    bed_destroy(args.target_regions);
    bed_destroy(args.design_regions);
    bed_destroy(args.predict_regions);    
    fai_destroy(args.fai);
    free(args.string.s);
}
int main(int argc, char **argv)
{
    if (argc == 0)
	return usage();
    
    if ( parse_args(--argc, ++argv) != 0 ) 
	return 1;
    
    generate_oligos();
    
    export_summary_reports();
    
    clean_memory();
    LOG_print("%u oligos were generated. See probe.txt.gz for details.", args.probes_number);
    LOG_print("Sucess.");
    return 0;
}
