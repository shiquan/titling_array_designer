#include "bed_utils.h"
#include <errno.h>

struct bed_readers *bed_readers_init()
{
    struct bed_readers *files = (struct bed_readers*)malloc(sizeof(struct bed_readers));
    return files;
}
int bed_readers_add_file(struct bed_readers *files, const char *fname)
{
    struct bed_file *file = (struct bed_file*)malloc(sizeof(struct bed_file));
    file->fp = gzopen(fname, "r");
    file->flag = BED_FILE_EMPTY;
    file->l_name = file->m_names = 0;
    file->names = 0;
    file->regions_count = file->length_total = 0;
    file->errno = 0;
    if (file->fp == NULL) {
	file->errno = errno;
	return 1;
    }
    file->hash = kh_init(reg);
    files->files = (struct bed_file**)realloc(files->files, (files->n_files+1) * sizeof(struct bed_file*));
    files->file_names = (char**)realloc(files->file_names, (files->n_files+1)* sizeof(char *));
    files->file_names[file->n_files] = strdup(fname);
    files->files[files->n_files] = file;
    files->n_files++;
    return 0;
}
void bed_readers_remove_file(struct bed_readers *files, int i)
{
    
}
void bed_chrom_destory(struct bed_chrom *chrom)
{
    int i;
    for (i = 0; i < chrom->l; ++i) {
	if (chrom->a[i].data.m) free(chrom->a[i].data.s);
    }
    free(chrom->a);
    free(chrom);
}
void bed_file_destory(struct bed_file *file)
{
    khiter_t k;
    int i;
    for (i = 0; i < file->l_names; ++i) {
	char *name = file->names[i];
	k = kh_get(reg, file->reghash, name);
	free(name);
	if (k == kh_end(reg)) {
	    continue;
	} else {
	    struct bed_chrom * chrom = kh_val(reghash, k);
	    bed_chrom_destory(chrom);
	    kh_del(reg, reghash, k);
	}
    }
}
void bed_readers_destroy(struct bed_readers *files)
{
    int i;
    for (i = 0; i < files->n_files; ++i) {
	free(files->file_names[i]);
	bed_file_destory(files->files[i]);
    }
    free(files->file_names);
    free(files->files);
    free(files);
}
int bed_readers_next_line(struct bed_readers *files, struct bed_line_core *line)
{
    
}
