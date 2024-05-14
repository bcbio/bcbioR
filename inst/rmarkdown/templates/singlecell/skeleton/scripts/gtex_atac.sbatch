#!/bin/bash

#SBATCH --job-name=Cellranger-Arc      # Job name
#SBATCH --partition=priority            # Partition name
#SBATCH --time=1-23:59                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=16             # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

module load cellranger-ARC/2.0.0

localcores=$SLURM_CPUS_PER_TASK
localmem=$SLURM_MEM_PER_NODE

cellranger-arc count --id=cellranger-arc_output \
                       --reference=/n/shared_db/GRCh38/uk/cellranger-ARC/2.0.0/2.0.0/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/n/data1/cores/bcbio/PIs/athe_tsibris/10x_TEA-Seq_Tsibris_human_HIV_ex-vivoCD4_hbc04730_ZZ/meta/libraries_3rdbatch.csv \
                       --localcores=16 \
                       --localmem=128
                       