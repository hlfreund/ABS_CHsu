#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=5-00:00:00     # 5 day
#SBATCH --output=Assign_Taxonomy_RawMGMs_metaphlan_11.7.24.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Assign_Taxonomy_RawMGMs_metaphlan_11.7.24"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/Burcin_MGMData

echo "Building the dir for MetaPhlan results"

if [[ ! -d ${workplace}/MetaPhlan_Results ]]; then
    mkdir ${workplace}/MetaPhlan_Results
fi

# using trimmed, host-read-filtered sequences (aka the raw reads were first trimmed, then human host sequences were removed, now ready for MetaPhlan)

for i in ${workplace}/Bowtie2_NoHumanReads_Results/*_host_removed_R1.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_host_removed_R*}
    
    # if the fastq files are still gzipped, gunzip them
 #   if [[ -f ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1.fastq.gz ]]; then
 #
 #       gunzip ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1.fastq.gz
 #       gunzip ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R2.fastq.gz
 #
#    fi
    
    # if metaphlan output files do not exist, run metaphlan
    if [[ ! -f ${workplace}/MetaPhlan_Results/${SAMPLE}_metagenome.bowtie2.bz2 ]] && [[ ! -f ${workplace}/MetaPhlan_Results/${SAMPLE}_profiled_metagenome.txt ]]; then
        
        # run metaphlan with F & R trimmed mgm fastq fiiles
        echo "Running metaphlan on" ${SAMPLE}
        
        metaphlan ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1.fastq,${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R2.fastq --bowtie2out ${workplace}/MetaPhlan_Results/${SAMPLE}_metagenome.bowtie2.bz2 --nproc 6 --input_type fastq -o ${workplace}/MetaPhlan_Results/${SAMPLE}_profiled_metagenome.txt

        
    fi
    
done

# use metaphlan's scrit to merge metaphlan outputs into single file
# merge_metaphlan_tables.py *_profiled_metagenome.txt > ${workplace}/MetaPhlan_Results/BO_LF_MetaPhlan_merged_abundance_table_10.31.24.txt
# ^^ if this fails, go to 3a_Merge_MetaPhlan_Outputs.sh

conda deactivate

