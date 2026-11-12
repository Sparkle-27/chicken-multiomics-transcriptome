# Whole-genome sequencing

| Script | Step | Main inputs |
| --- | --- | --- |
| `1_fastp.sh` | Read filtering | Paired FASTQ files |
| `2_bwa_align.sh` | BWA-MEM alignment and BAM sorting | Indexed reference FASTA and clean reads |
| `3_haplotypecaller.sh` | Per-sample variant calling | Reference FASTA and BAM |
| `4_snp_filter.sh` | SNP hard filtering | Reference FASTA and SNP-only VCF |

Run scripts without arguments to display usage. Example:

```bash
bash scripts/wgs/1_fastp.sh sample_R1.fq.gz sample_R2.fq.gz work/sample
bash scripts/wgs/2_bwa_align.sh reference.fa work/sample_R1.fq.gz work/sample_R2.fq.gz sample 16 work/sample
bash scripts/wgs/3_haplotypecaller.sh reference.fa work/sample.sort.bam work/sample.vcf.gz
```

BWA-MEM uses `-M` and a sample read group. VariantFiltration marks SNPs with QD < 5, FS > 60, SOR > 3, MQ < 40, MQRankSum < −12.5 or ReadPosRankSum < −8. FILTER annotations are retained in the output; cohort-level merging, joint genotyping and PASS selection are separate steps.

Version records: fastp 0.23.2, BWA 0.7.17, SAMtools 1.11 and GATK 4.4.0.0. Earlier variant-calling commands referenced GRCg6a, while the transcriptome analysis uses GRCg7b. Reference FASTA files and BAM/VCF coordinates must use the same assembly; these scripts do not convert coordinates between assemblies.
