# Short-read RNA-seq

| Script | Step | Main inputs |
| --- | --- | --- |
| `1_fastp.sh` | Read filtering | Paired FASTQ files |
| `2_star_index.sh` | STAR index generation | Reference FASTA and GTF |
| `3_star_align.sh` | Paired-end alignment | STAR index and clean reads |
| `4_mark_duplicates.sh` | Duplicate marking | Sorted BAM |
| `5_stringtie_quant.sh` | Reference-guided quantification | BAM and annotation GTF |
| `6_prepDE_counts.sh` | Gene/transcript count matrices | prepDE.py3 and sample GTF list |

Run scripts without arguments to display usage. Reference files use GRCg7b and Ensembl release 110.

The STAR index uses `sjdbOverhang=149` and `genomeSAindexNbases=13`. Alignment retains the splice-junction settings in script 3, a maximum intron/mate gap of 500,000 bp and a maximum mismatch/read-length ratio of 0.07. Duplicate marking retains the reads; apply this step according to the BAM preprocessing used for quantification. StringTie uses `-e` with the reference annotation. prepDE requires the average read length as an argument.

Version records: fastp 0.23.2, STAR 2.7.10a and SAMtools 1.11. StringTie is recorded as 2.1.7 in the script path and 2.2 in the manuscript table; its exact version remains unconfirmed.

The downstream [differential-expression analysis](../differential_expression/README.md) starts from the processed matrices in [results/gene_counts](../../results/gene_counts).
