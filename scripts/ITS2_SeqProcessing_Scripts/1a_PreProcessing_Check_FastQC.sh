#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G  #
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --output=QualityCheck_16SV3V4_12.6.24.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="16S V3V4 FastQC quality check"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate AmpliconAnalysis

# workplace will be the directory where your sequence data are stored
workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_ITS2_11.29.2024/ITS2_NoAdapters

if [[ ! -d ./FastQC_Results ]]; then
    mkdir FastQC_Results
    mkdir FastQC_Results/16S_FastQC
    mkdir FastQC_Results/ITS2_FastQC
fi


for FILE in 16S_Seqs/*.fastq.gz;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    
    fastqc $FILE --outdir=./FastQC_Results/16S_FastQC
    
done

for FILE in ITS2_Seqs/*.fastq.gz;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    fastqc $FILE --outdir=./FastQC_Results/ITS2_FastQC
    
done

conda deactivate
