# titling array designer - a lightweight oligonucletide designer for target resequencing panels



## generate_oligos

**generate_oligos** is the core program to generate oligos with defined target regions (in bed format) and genome reference.

```
generate_probes - generate DNA titling array probes from target regions and reference genome sequences.
Usage:
generate_probes [options] -fasta hg19.fa -target target.bed -uniq_regions database.bed.gz -outdir output_dir
Options:
  -r, -fasta [fasta file]
            reference genome sequences, in fasta format.
  -t, -target [bed file]
            target bed file to design.
  -o, -outdir [dir]
            output directary, will create it if not exists.
  -u, -database [tabix-indexed bed file]
            database of non-repeats, or user pre-defined designable regions.
  -l, -length [50]
            pre-defined oligo length, usually from 50 base to 90 base, set 0 for dynamic design.
  -d, -depth [2]
            oligo depths pre base, increase this value will increase the dense of oligos.
  -p, -project [string]
            project id
  -must_design
            if no uniq regions around small target, must design it no matter repeat regions.
  -h, -help
            for help information.
Version: 6c9c2d4-dirty
Bugs report: shiquan@genomics.cn
Homepage: https://github.com/shiquan/titling_array_designer
```

About the parameters:
* **-p**, project id, this is mandatory for user to record the poject information;
* **-t**, specify target regions, all the regions should be formated in BED, please notice that all the start coordinate is 0 based and end coordinate is 1 based in BED file.
* **-r**, specify the reference genome in FASTA format. And the reference database should be indexed with `samtools faidx` , to make sure your data are properly indexed, please check the `*.fai` file in the same directory.
* **-u**, designable region database. This database tell program what exactly regions could be *designed*, please notice that it is not a mandatory database, but if you set, all the oligos should be covered by the regions in the database.


Output files include:
* **design_regions.bed**, oligos covered regions in BED format.
* **target_regions.bed**, target regions to design, this file may be slightly different with your specified target regions, because program will round any small region to 100nt for better performance. And this file is the final target regions after region-check.
* **


About the format of probe file, please refer to [format of probe file](https://github.com/shiquan/titling_array_designer/blob/master/document/format.md).

## merge_oligos

