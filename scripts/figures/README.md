# Figures and QC scripts

1. `1_sample_qc.R`: blind VST, global and within-tissue top500 PCA, sample distances and tissue-centroid diagnostics.
2. `2_plot_gene_figures.R`: 29 gene-expression/QC figures from fitted results and QC tables.
3. `3_combine_figures.py`: combined PDF with page index and bookmarks.
4. `4_catalog_summary.R`: supplied final-catalog classifications, exon and length summaries.
5. `5_length_summary.R`: supplied positive read/transcript lengths, median/IQR and N50.
6. `6_aqp9_structure.R`: reference/PacBio AQP9 exon tracks from supplied GTFs.
7. `7_aqp9_expression.R`: descriptive plot of supplied AQP9 TPM data.

Scripts 1–3 are executed by the gene-analysis runner. Scripts 4–7 are standalone helpers requiring their stated inputs; they do not create or validate a study-wide unified catalog.

PCA uses the original UMAP group colors and symbols matched by sample ID, with coordinates recomputed from VST counts. Choosing 500 variable genes is for PCA only, not a restriction of gene testing to 500 features.

Volcano plots retain the original 8×8-inch layout, size 2 points, blue #5CA7BA/grey/red #FE4365 colors, typography and bottom legend. X=log2(TS/CB), Y=-log10(nominal P); significance colors use BH padj<0.05. Red means TS higher, blue CB higher. Stable denotes not FDR-significant, not demonstrated equality.

AQP9 total expression, cloned intron-retention structures and differential isoform usage are separate questions. The helpers do not establish isoform switching.
