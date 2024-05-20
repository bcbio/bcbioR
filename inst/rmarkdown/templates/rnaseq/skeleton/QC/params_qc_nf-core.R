# info params

# This is the file used to run nf-core or compatible to that
metadata_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/coldata.csv'
# This file is inside star_salmon/ folder
# remove URL below to point to local files
se_object=url('https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/salmon.merged.gene_counts.rds')
# This folder is in the output directory inside multiqc folder
multiqc_data_dir='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/multiqc-report-data/'
# This file is inside the genome folder in the output directory
gtf_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/genome.filtered.gtf.gz'
