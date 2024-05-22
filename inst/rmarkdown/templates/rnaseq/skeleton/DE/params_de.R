# project params
date = "YYYYMMDD"
basedir <- './' # where to write down output files

# params for bcbio
# coldata_fn = "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/coldata.csv"
# counts_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/tximport-counts.csv'
# se_object=url("https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/bcbio-se.rds")
#

# params for nfcore
# Your data
# This is the file used to run nf-core or compatible to that
coldata_fn='/Path/to/metadata/meta.csv'
# This file is inside star_salmon/ folder
counts_fn='/path/to/nf-core/output/star_salmon/salmon.merged.gene_counts.tsv'
# This folder called "multiqc_report_data" is inside the output directory star_salmon inside multiqc folder
multiqc_data_dir='/path/to/nf-core/output/star_salmon/multiqc_report_data'
# This file is inside the genome folder in the output directory, use this only non-model organism
# gtf_fn='/path/to/nf-core/output/genome/hg38.filtered.gtf'


# Example data: COMMENT THESE LINE IF YOU ARE USING YOUR DATA
coldata_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/coldata.csv'
counts_fn=url('https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/star_salmon/salmon.merged.gene_counts.tsv')
# This folder is in the output directory inside multiqc folder
multiqc_data_dir='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/multiqc/star_salmon/multiqc-report-data/'
# This file is inside the genome folder in the output directory
gtf_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/devel/nf-core/genome/genome.filtered.gtf.gz'
se_object = NA
