#include "utils.h"
#include "file.h"
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "zlib.h"
#include "htslib/kseq.h"

#define KSTRING_INIT {0, 0, 0}

int usage()
{
    fprintf(stderr, 
"- Details: Merge probes files into one.\n"
"- Usage: merge_probes  probes1.txt.gz probes2.txt.gz [ probes3.txt.gz ... ]\n"
"- Author: Shi Quan (shiquan@genomics.cn)\n"
);
    return 1;
}
struct args {
    int n_files;
    int m_files;
    htsFile **fp;
    kstring_t string;
    int oligo_length;
    // 0 on default, 1 on header only, 2 on no header
    int header_flag;
} args = {
    .n_files = 0,
    .m_files = 0,
    .fp = NULL,
    .string = {0, 0, 0},
    .oligo_length = 0,
    .header_flag = 0,
};

int read_comment_line(htsFile *fp, kstring_t *string)
{
    for ( ;; ) {
        if ( hts_getline(fp, KS_SEP_LINE, string) < 0 )
            return 1;
        if ( string->l == 0 )
            continue;
        if ( string->s[0] == '#')
            break;
        return 1;
    }
    return 0;
}
int read_line(htsFile *fp, kstring_t *string)
{
    for ( ;; ) {
        if ( hts_getline(fp, KS_SEP_LINE, string) < 0 )
            return 1;
        if ( string->l == 0 )
            continue;
        if ( string->s[0] == '#' || string->s[0] == '/')
            continue;
        break;
    }
    return 0;
}
void release_memory()
{
    int i;
    for ( i = 0; i < args.n_files; ++i )
        hts_close(args.fp[i]);
    free(args.fp);
    if ( args.string.l )
        free(args.string.s);
}
int parse_length(htsFile *fp)
{
    kstring_t string = KSTRING_INIT;
    int length = 0;
    if ( file_seek(fp, 0, SEEK_SET) < 0 )
        return 0;
    
    for ( ;; ) {
        if ( read_comment_line(fp, &string) )
            return length;
        if ( strncmp(string.s, "##max_length=", 13) == 0 ) {
            sscanf(string.s+13, "%d", &length);
            return length;
        }
    }
    return 0;
}

int parse_args(int argc, char **argv)
{
    
    int i;    
    kstring_t command = KSTRING_INIT;
    for ( i = 0; i < argc; ++i ) {
        if ( i )
            kputc(' ', &command);
        kputs(argv[i], &command);
    }
    if ( argc == 0 )
        return usage();
    if ( command.l )
        free(command.s);
    
    for ( i = 0; i < argc; ) {    
        const char *a = argv[i++];
        if ( strcmp(a, "-h") == 0 ) {
            args.header_flag = 1;
            continue;
        } else if ( strcmp(a, "-H") == 0 ) {
            args.header_flag = 2;
            continue;
        }
        if ( args.m_files == args.n_files ) {
            args.m_files += 2;
            args.fp = (htsFile**)realloc(args.fp, args.m_files * sizeof(htsFile*));
        }
        args.fp[args.n_files] = hts_open(a, "r");
        if ( args.fp[args.n_files] == NULL ) {
            warnings("%s : %s", a, strerror(errno));
            continue;
        }
        args.n_files++;
    }
    if ( args.n_files < 2)
        error("Must merge at least two probes files. %d", args.n_files);

    for ( i = 0; i < args.n_files; ++i ) {
        int length = 0;
        length = parse_length(args.fp[i]);
        if (length > args.oligo_length )
            args.oligo_length = length;
    }
    
    return 0;
}
void merge_probes()
{
    int i;
    kstring_t string = KSTRING_INIT;
    // Generate header.
    for ( i = 0; i < args.n_files; ++i ) {
        if ( file_seek(args.fp[i], 0, SEEK_SET) < 0 )
            continue;
        for ( ;; ) {
            if ( read_comment_line(args.fp[i], &string) )
                break;
            if ( i == 0 ) {
                if ( string.s[1] == '#') {
                    kputc('\n', &string);
                    printf("%s", string.s);
                }
            } else if ( strncmp(string.s, "##Command=", 10 ) == 0 ) {
                kputc('\n', &string);
                printf("%s", string.s);
            } else if ( i == args.n_files-1 && string.s[1] != '#') {
                kputc('\n', &string);
                printf("%s", string.s);
            }
        }        
    }
    // Generate body.
    for ( i = 0; i < args.n_files; ++i ) {
        if ( file_seek(args.fp[i], 0, SEEK_SET) < 0 )
            continue;
        for ( ;; ) {
            if ( read_line(args.fp[i], &string) )
                break;
            kputc('\n', &string);
            printf("%s", string.s);
        }
    }    
}
int main(int argc, char **argv)
{
    if ( parse_args(--argc, ++argv) )
        return 1;

    merge_probes();
    release_memory();
    return 0;
}
