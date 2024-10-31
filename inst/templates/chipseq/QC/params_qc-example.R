# info params
# Example data

# this is the samplesheet used to run nf-core, with additional columns containing covariates of interest
coldata_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/refs/heads/main/chipseq/chipseq_peakanalysis_H3K27Ac.csv'

# This folder is in the output directory inside multiqc folder
multiqc_data_dir='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/refs/heads/main/chipseq/multiqc/narrowPeak/multiqc_data/'

# This folder is in the output directory
# peaks_dir = "https://api.github.com/repos/bcbio/bcbioR-test-data/contents/chipseq/bowtie2/mergedLibrary/macs2/narrowPeak"

# This folder is in the output directory
counts_fn = url('https://raw.githubusercontent.com/bcbio/bcbioR-test-data/refs/heads/main/chipseq/bowtie2/mergedLibrary/macs2/narrowPeak/consensus/H3K27ac/deseq2/H3K27ac.consensus_peaks_small.rds')
