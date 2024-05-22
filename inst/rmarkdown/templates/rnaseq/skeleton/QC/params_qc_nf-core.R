# info params

# Your data
# This is the file used to run nf-core or compatible to that
metadata_fn='/Path/to/metadata/meta.csv'
# This file is inside star_salmon/ folder
se_object='/path/to/nf-core/output/star_salmon/salmon.merged.gene_counts.rds'
# This folder called "multiqc_report_data" is inside the output directory star_salmon inside multiqc folder
multiqc_data_dir='/path/to/nf-core/output/star_salmon/multiqc_report_data'
# This file is inside the genome folder in the output directory, use this only non-model organism
# gtf_fn='/path/to/nf-core/output/genome/hg38.filtered.gtf'

# Example data: COMMENT THESE LINE IF YOU ARE USING YOUR DATA
metadata_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/coldata.csv'
se_object=url('https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/star_salmon/salmon.merged.gene_counts.rds')
# This folder is in the output directory inside multiqc folder
multiqc_data_dir='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/multiqc/star_salmon/multiqc-report-data/'
# This file is inside the genome folder in the output directory
gtf_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/devel/nf-core/genome/genome.filtered.gtf.gz'
