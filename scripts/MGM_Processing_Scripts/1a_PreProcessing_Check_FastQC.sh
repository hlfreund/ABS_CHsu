#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G  #
#SBATCH --time=2-00:00:00     # 2 day
#SBATCH --output=QualityCheck_ABS_MGMs.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="QualityCheck_ABS_MGMs"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow

# workplace will be the directory where your sequence data are stored
workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project

if [[ ! -d ${workplace}/FastQC_Results ]]; then
    mkdir ${workplace}/FastQC_Results
    mkdir ${workplace}/FastQC_Results/MGM_PreTrim_FastQC
fi

for FILE in ${workplace}/ABS_MGMs/*.fastq.gz;
do
    #f=$(basename $FILE)
    #SAMPLE=${f%.fastq}
    
    fastqc $FILE --outdir=${workplace}/FastQC_Results/MGM_PreTrim_FastQC/
    
done

# run MultiQC on these results

multiqc ${workplace}/ABS_MGMs ${workplace}/FastQC_Results/MGM_PreTrim_FastQC/ -o ${workplace}/FastQC_Results/MGM_PreTrim_FastQC/
# ^ format: multiqc <directories to include> <directory to store results>

conda deactivate
