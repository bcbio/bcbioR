# info params

coldata_fn='/path/to/nf-core/samplesheet.csv'

# This folder is in the nf-core output directory, maybe is broadPeak instead of narrowPeak
peaks_dir = '/path/to/nf-core/output/bowtie2/mergedLibrary/macs2/narrowPeak/'

# This folder is in the nf-core output directory
bam_dir = '/path/to/nf-core/output/bowtie2/mergedLibrary/'

# this will be the file that the diffbind samplesheet is eventually saved in 
diffbind_samplesheet_fn = 'diffbind_samplesheet.csv'

# This will be the file that the diffbind counts matrix is eventually saved in
diffbind_counts_fn = 'diffbind_counts.qs'