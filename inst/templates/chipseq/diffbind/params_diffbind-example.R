# info params
# Example data

# this is the samplesheet used to run nf-core, with additional columns containing covariates of interest
coldata_fn='https://raw.githubusercontent.com/bcbio/bcbioR-test-data/refs/heads/main/chipseq/chipseq_peakanalysis_H3K27Ac.csv'

# example data doesn't need this but this folder is in the nf-core output directory, maybe is broadPeak instead of narrowPeak
# peaks_dir = '/path/to/nf-core/output/bowtie2/mergedLibrary/macs2/narrowPeak/'

# example data doesn't need this but this folder is in the nf-core output directory
# bam_dir = '/path/to/nf-core/output/bowtie2/mergedLibrary/'

diffbind_samplesheet_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/refs/heads/main/chipseq/diffbind_samplesheet.csv'
diffbind_counts_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/refs/heads/main/chipseq/diffbind_counts.qs'
