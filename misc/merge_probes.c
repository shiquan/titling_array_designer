#include "utils.h"
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "htslib/string.h"

int usage()
{
    fprintf(stderr, "
- Details: Merge probes files into one.
- Usage: merge_probes  probes1.txt.gz probes2.txt.gz [ probes3.txt.gz ... ]
- Author: Shi Quan (shiquan@genomics.cn)
");
    return 1;
}
struct args {
    int n_files;
    char **fnames;
    
};
int parse_args(int argc, char **argv)
{
    return 0;
}
void merge_probes()
{
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    merge_probes();
    return 0;
}
