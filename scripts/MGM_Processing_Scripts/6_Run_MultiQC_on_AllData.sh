#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G  #
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --output=MultiQC_on_AllData.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Run_MultiQC_on_AllData"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

ABSprojdir=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project
workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project/ABS_MGMs

# run MultiQC on entire data set

multiqc ${ABSprojdir}/FastQC_Results/MGM_PreTrim_FastQC/ ${ABSprojdir}/FastQC_Results/MGM_PostTrim_FastQC/ ${workplace}/Trimmed_Seqs ${workplace}/Bowtie2_NoHumanReads_Results ${workplace}/Humann3_Results . -o ${workplace}

conda deactivate

