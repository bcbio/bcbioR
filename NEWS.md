# bcbioR 0.4.4


# bcbioR 0.4.3

migrate rnaseq templates to its own repo #120
add HD visium QC #101
add HD visum clustering #102
include code to generate files for DEGpattern #116
double code for single cell #130
adapt to QC to ATACseq #133
fix error in names #132
fix GSEA rnaseq error and styling code bcbio/rnaseq-reports#2

# bcbioR 0.4.2

Filter genes according biotype or gene name #117
Add Visium HD QC #101 #115
Add degPattern template #97
Collapse multiple params files #65
Add rmaseq library links #60
Fix version comparison #104
Save log object in sc data #105
More docs for clusters ids #106

# bcbioR 0.4.1

- in DEG.Rmd, write expression table only once	https://github.com/bcbio/bcbioR/issues/63
- add methods to RNA QC and DE templates	https://github.com/bcbio/bcbioR/issues/61
- sc DE pseudobulk
- sc DE MAST
- compositional analysis scRNA
- entrezid for enrichment analysesS	v0.4.1	P1
- subset data based on metadata file.	https://github.com/bcbio/bcbioR/issues/64
- ATACseq QC report
- Add WGCNA analysis to RNAseq
- bug, should be !is.na	https://github.com/bcbio/bcbioR/issues/73
- line of code in volcano plot causes odd figures	https://github.com/bcbio/bcbioR/issues/70
- volcano plot colors	https://github.com/bcbio/bcbioR/issues/74
- remove inline option in RMD	https://github.com/bcbio/bcbioR/issues/78
- change x lab to reads	https://github.com/bcbio/bcbioR/issues/75
- Visium RMD template	
- fix mouse genome for nfcore and templates
- adapt to mouse and another genome in the DE RNAseq template	https://github.com/bcbio/bcbioR/issues/80
- Change to Annotation Hub	https://github.com/bcbio/bcbioR/issues/94

# bcbioR 0.3.1

- Fix bugs in RNAseq
- Add function to set up library

# bcbioR 0.3.0

* re-structure templates
* Add text with best practices
* Reproducibility:
  * test data for RNAseq, singlecell, CHIPseq
* Base project:
  * Guidelines to create repo easily
  * Example Rmd with headers and aesthetics
* RNASEQ
  * Use provenance for FA in DE report
  * Support multiple comparisons
* New templates:
  * methylation - draft
  * singcell cell QC and Inegration - stable
  * scQC shiny app - stable
  * chipseq QC and Diffbind - beta
  * COSMX - draft

# bcbioR 0.1.3

* fix duplicated gene names
* Add combatseq
* Add draft templates for TEAseq, COSMX
* Adapt templates to nf-core rnaseq
* Fix when sample start by number
* Fix when rRNA rate is missing
* Add by sample plots in QC
* Add function to check nfcore samplesheet
* Add PCA with variance analysis
* Minor details fixed in QC RNAseq report
* Correct metric for rRNA and tRNA

# bcbioR 0.1.1

* Add page to github
* Make CHECK to pass
* Add vignette
