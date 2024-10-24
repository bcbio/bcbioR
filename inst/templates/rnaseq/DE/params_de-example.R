# project params
date = "YYYYMMDD"
basedir <- './' # where to write down output files

# params for bcbio
# coldata_fn = "https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/coldata.csv"
# counts_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/tximport-counts.csv'
# se_object=url("https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/bcbio/bcbio-se.rds")
#

# Example data
coldata_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/coldata.csv'
counts_fn=url('https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/star_salmon/salmon.merged.gene_counts.tsv')
# This folder is in the output directory inside multiqc folder
multiqc_data_dir='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/devel/rnaseq/nf-core/multiqc/star_salmon/multiqc-report-data/'
# This file is inside the genome folder in the output directory
gtf_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/devel/nf-core/genome/genome.filtered.gtf.gz'
se_object = NA
