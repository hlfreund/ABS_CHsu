#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G  #
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --output=ErrorTrim_ITS2_eestats_12.2.24.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="ErrorTrim_ITS2_eestats_12.2.24"
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

if [[ ! -d ${workplace}/EEstats_Results ]]; then
    mkdir ${workplace}/EEstats_Results
    #mkdir ${workplace}/EEstats_Results/16S_EEstats
    mkdir ${workplace}/EEstats_Results/ITS2_EEstats
fi

## * be sure to gunzip files before running script

#for FILE in ${workplace}/16S_Seqs/*.fastq;
#do
#    f=$(basename $FILE)
#    SAMPLE=${f%.fastq*}
#
#    #usearch -fastq_eestats $FILE -output ${SAMPLE}_eestats.txt
#    usearch -fastq_eestats2 $FILE -output ${workplace}/EEstats_Results/16S_EEstats/${SAMPLE}_eestats2.txt
#
#    #mv ${SAMPLE}_eestats2.txt EEstats_Results/16S_EEstats
#done

for FILE in ${workplace}/*.fastq;
do
    f=$(basename $FILE)
    SAMPLE=${f%.fastq*}
    
    #usearch -fastq_eestats $FILE -output ${SAMPLE}_eestats.txt
    usearch -fastq_eestats2 $FILE -output ${workplace}/EEstats_Results/ITS2_EEstats/${SAMPLE}_eestats2.txt
    
    #mv ${SAMPLE}_eestats2.txt EEstats_Results/ITS2_EEstats
    
done

# More info here: https://www.drive5.com/usearch/manual/cmd_fastq_eestats.html
# https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html

conda deactivate
