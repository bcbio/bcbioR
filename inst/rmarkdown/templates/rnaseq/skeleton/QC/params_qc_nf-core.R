# info params

# This is the file used to run nf-core or compatible to that
metadata_fn='/Path/to/metadata/meta.csv'
# This file is inside star_salmon/ folder
se_object=('/path/to/nf-core/output/salmon/salmon.merged.gene_counts.rds')
# This folder called "multiqc_report_data" is inside the output directory star_salmon inside multiqc folder
multiqc_data_dir='/path/to/nf-core/output/star_salmon/multiqc_report_data'
# This file is inside the genome folder in the output directory
gtf_fn='/path/to/nf-core/output/genome/hg38.filtered.gtf'
# Put hg38, mm10, mm39, or other
genome='hg38'
