#!/bin/env bash

#SBATCH --mem=40GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=Ryan.Johnson@nyulangone.org
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --job-name=nextflow_chip

module load nextflow/23.04.4
module load r/4.3.2
# just checking if the pe_reads parameter can be seen form the workflow section.
# it can which means that this will work for any pair end reads for any project
# ues -entry parameter followed by the name of the workflow to run
nextflow run chip_seq_nf_pipeline.nf -resume --se_reads '../chip_fastqs/chip_*.fastq'
