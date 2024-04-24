# Guideline for RNAseq downstream analysis

# DropBox

- In `reports/QC`
  - [ ] copy `bcbio-se.rds` and `tximport-counts.csv`
  - [ ] copy QC `Rmd/R/html/figures`
- In `reports/DE`
  -	[ ] Normalized counts for all genes x all samples (csv format)
- In `reports/DE`, for *each analysis*:
  - **Note** For multiple comparisons/analysis, do a single report/template if possible in the parent folder using parameters whenever possible. 
  - Create a folder with the comparison names in the files. Numbering by comparison (`01.1_DE_comp1`, `01.2_DE_comp2`, etc.). If you’re running multiple models for the same comparison, append `_M#`. Add the following files under each folder:
  - [ ] Normalized count table with the samples used in this analysis/comparison.
  -	[ ] Full results `DESeq2` for all genes (csv format) with annotation columns appended. 
  -	[ ] Significant genes results file (subset of annotated full results by chosen p-value and LFC). Separate files will be created for each individual contrast.
  -	[ ] Significant genes results file as described above, but additionally append columns containing normalized count values for each sample.
  -	Make sure to append the gene symbols to these tables so the researcher can interpret the results.

# GitHub

- [ ] Push all `*Rmd` `*R` files used for the *QC* and *DE* analysis respecting folder structure.

Please, ignore `*html/figures/csv` and any output of the code.
