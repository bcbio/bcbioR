# info params


coldata_fn='/path/to/nf-core/samplesheet.csv'
# This folder is in the nf-core output directory inside multiqc folder
multiqc_data_dir='/path/to/nf-core/output/multiqc/narrowPeak/multiqc_data/'
# This folder is in the nf-core output directory, maybe is broadPeak instead of narrowPeak
peaks_dir = '/path/to/nf-core/output/bowtie2/mergedLibrary/macs2/narrowPeak/'
# This folder is in the nf-core output directory, maybe is broadPeak instead of narrowPeak, also includes antibody name
counts_fn = '/path/to/nf-core/output/bowtie2/mergedLibrary/macs2/narrowPeak/consensus/antibody/deseq2/antibody.consensus_peaks.rds'
