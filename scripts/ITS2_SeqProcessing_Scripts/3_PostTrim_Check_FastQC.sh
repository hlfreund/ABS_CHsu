#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G  #
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --output=CheckQuality_Trimmed_ITS2Seqs_12.2.24.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="CheckQuality_Trimmed_ITS2Seqs_12.2.24"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

source /tscc/nfs/home/lfreund/miniforge3/bin/activate AmpliconAnalysis

# workplace will be the directory where your sequence data are stored
workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_ITS2_11.29.2024/ITS2_NoAdapters

if [[ ! -d ${workplace}/Trim_FastQC_Results ]]; then
    mkdir ${workplace}/Trim_FastQC_Results
    #mkdir ${workplace}/Trim_FastQC_Results/16S_Trimmed_FastQC
    mkdir ${workplace}/Trim_FastQC_Results/ITS2_Trimmed_FastQC
fi

#for FILE in ${workplace}/Trimmed_Seqs/16S_Trimmed/*.fastq;
#do
#    f=$(basename $FILE)
#    SAMPLE=${f%.fastq*}
#
#    fastqc $FILE --outdir=${workplace}/Trim_FastQC_Results/16S_Trimmed_FastQC
#
#done

for FILE in ${workplace}/Trimmed_Seqs/ITS2_Trimmed/*.fastq;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    
    fastqc $FILE --outdir=${workplace}/Trim_FastQC_Results/ITS2_Trimmed_FastQC
    
done

conda deactivate
