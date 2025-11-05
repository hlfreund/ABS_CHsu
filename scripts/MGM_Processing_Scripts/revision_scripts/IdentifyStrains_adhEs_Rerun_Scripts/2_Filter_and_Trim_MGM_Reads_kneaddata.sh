#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=5-00:00:00     # 5 day
#SBATCH --output=Filter_and_Trim_Reads_Kneaddata.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Filter_and_Trim_Reads_Kneaddata"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs
knead_db=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Kneaddata_Resources

if [[ ! -d ${workplace}/Kneaddata_Results ]]; then
    mkdir ${workplace}/Kneaddata_Results

fi

gunzip ${workplace}/*.gz
#if [[ -f ${workplace}/*.fastq.gz ]]; then
#    gunzip ${workplace}/*.fastq.gz
#fi

# Kneaddata can trim the reads and remove host genome reads before taxonomic and functional annotation

for i in ${workplace}/*_R1_001.fastq;
do
    
    f=$(basename $i)
    SAMPLE=${f%_R1_001.fastq}
    
    # if results dir called ${SAMPLE}_knead_files does not exist, then run kneaddata
    if [[ ! -d ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files ]]; then
        
        kneaddata --input1 ${workplace}/${SAMPLE}_R1_001.fastq --input2 ${workplace}/${SAMPLE}_R2_001.fastq -db ${knead_db} -t 6 --reorder --store-temp-output --output-prefix ${SAMPLE} --output ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files
        
        mv ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_*.fastq ${workplace}/Kneaddata_Results/
    
    fi
   
    
done


conda deactivate

# More info on KneadData: https://github.com/biobakery/kneaddata?tab=readme-ov-file#paired-end-run

## Databases for Kneaddata - Downloading Notes - Important!
# Downloading Kneaddata database  (see comment here for details;
## tldr it’s a permissions issue and we have to download it this way on the university server)
# Must do the following before you use Humann:
# > wget https://g-227ca.190ebd.75bc.data.globus.org/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
# > tar -xvzf Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz

## Kneaddata Output Files explained:
# When performing quality filtering and trimming for paired end files, three things can happen:
## Both reads in the pair pass.
## The read in the first mate passes, and the one in the second does not pass.
## The read in the second mate passes, and the one in the first does not pass.

# seq_kneaddata.log: Log file containing statistics about the run.
# seq_kneaddata_paired_1.fastq: Reads from the first mate (if both reads in pair pass) identified as NOT belonging to any of the reference databases.
# seq_kneaddata_paired_2.fastq: Reads from the second mate (if both reads in pair pass) identified as NOT belonging to any of the reference databases.
# seq_kneaddata_unmatched_1.fastq: Reads from the first mate (read passed from R1, not R2) identified as NOT belonging to any of the reference databases.
# seq_kneaddata_unmatched_2.fastq: Reads from the second mate (read passed from R2, not R1)  identified as NOT belonging to any of the reference databases.

## One last important note...
## Kneaddata is particular about the format of FASTQ files and how it identifies the read pairs...
## More info and some solutions here: https://forum.biobakery.org/t/all-paired-end-read-unmatched/2895
