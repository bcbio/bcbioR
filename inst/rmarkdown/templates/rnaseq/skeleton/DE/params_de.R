# project params
date = "YYYYMMDD"
basedir <- './' # where to write down output files

# params for bcbio
# coldata_fn = "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/coldata.csv"
# counts_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/tximport-counts.csv'
# se_object=url("https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/bcbio-se.rds")
#

# params for nfcore
# This is the file used to run nf-core or compatible to that
coldata_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/coldata.csv'
# This file is inside star_salmon/ folder
counts_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/salmon.merged.gene_counts.tsv'
# This folder is in the output directory inside multiqc folder
multiqc_data_dir = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/multiqc-report-data/'
# This file is inside the genome folder in the output directory
gtf_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/genome.filtered.gtf.gz'
se_object = NA
