#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=6-00:00:00     # 5 day
#SBATCH --output=Assign_Taxonomy_RawMGMs_ABS_and_HCs_strainphlan_part2.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Assign_Taxonomy_RawMGMs_ABS_and_HCs_strainphlan_part2"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs

echo "Building the dir for MetaPhlan results"

if [[ ! -d ${workplace}/MetaPhlan_Results ]]; then
    mkdir ${workplace}/MetaPhlan_Results
    mkdir ${workplace}/MetaPhlan_Results/metaphlan_sams
    mkdir ${workplace}/MetaPhlan_Results/consensus_markers
    mkdir ${workplace}/MetaPhlan_Results/db_markers
fi

# using trimmed, host-read-filtered sequences (aka the raw reads were first trimmed, then human host sequences were removed, now ready for MetaPhlan)
# filtering and trimming was performed by KneadData

for i in ${workplace}/Kneaddata_Results/*_paired_1.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_paired_1*}
    
    # if metaphlan output files do not exist, run metaphlan
    if [[ ! -f ${workplace}/MetaPhlan_Results/${SAMPLE}_metagenome.bowtie2.bz2 ]] && [[ ! -f ${workplace}/MetaPhlan_Results/${SAMPLE}_profiled_metagenome.txt ]]; then
        
        # run metaphlan with F & R trimmed mgm fastq fiiles
        echo "Running metaphlan on" ${SAMPLE}
        
        metaphlan ${workplace}/Kneaddata_Results/${SAMPLE}_paired_1.fastq,${workplace}/Kneaddata_Results/${SAMPLE}_paired_2.fastq -s ${workplace}/MetaPhlan_Results/metaphlan_sams/${SAMPLE}_mgm.sam.bz2 --bowtie2out ${workplace}/MetaPhlan_Results/${SAMPLE}_metagenome.bowtie2.bz2 --nproc 6 --input_type fastq -o ${workplace}/MetaPhlan_Results/${SAMPLE}_profiled_metagenome.txt

        # create consensus-marker files to be used eventually for StrainPhlan
        ## this step will reconstruct all species found in it and store them in a pickle file (*.pkl). These are "sample-reconstructed strains"
        
         sample2markers.py -i ${workplace}/MetaPhlan_Results/metaphlan_sams/${SAMPLE}_mgm.sam.bz2 -o ${workplace}/MetaPhlan_Results/consensus_markers --nproc 6
         

    fi
    
    
done


# must have Python module loaded for these scripts to work!

module load python # need this for strainphlan .py functions

# Extract strain-specific markers (Ecoli) from MetaPhlan db
#extract_markers.py -c t__SGB10068 -o ${workplace}/MetaPhlan_Results/db_markers/

# Extract strain-specific markers (Klebsiella pneumoniae) from MetaPhlan db
#extract_markers.py -c t__SGB10115_group -o ${workplace}/MetaPhlan_Results/db_markers/

if [[ ! -d ${workplace}/MetaPhlan_Results/StrainPhlan_Results_v2/ ]]; then
    mkdir ${workplace}/MetaPhlan_Results/StrainPhlan_Results_v2

fi

# StrainPhlan uses all samples to pull out clade markers to generate MSA and phylo tree, tracing each sample in the process

echo "Running StrainPhlan on all samples"

# first find Ecoli clade markers
strainphlan -s ${workplace}/MetaPhlan_Results/consensus_markers/*.json.bz2 -m ${workplace}/MetaPhlan_Results/db_markers/t__SGB10068.fna -o ${workplace}/MetaPhlan_Results/StrainPhlan_Results_v2 -n 6 -c t__SGB10068 --mutation_rates

# then find Klebsiella pneumoniae clade markers
strainphlan -s ${workplace}/MetaPhlan_Results/consensus_markers/*.json.bz2 -m ${workplace}/MetaPhlan_Results/db_markers/t__SGB10115_group.fna -o ${workplace}/MetaPhlan_Results/StrainPhlan_Results_v2 -n 6 -c t__SGB10115_group --mutation_rates

echo "./StrainPhlan_Results_v1 contains StrainPhlan results for the minimum 10% of markers kept after filtering to keep a sample" > ${workplace}/MetaPhlan_Results/StrainPhlan_Results_v1/README.md

conda deactivate

# use metaphlan's scrit to merge metaphlan outputs into single file
# merge_metaphlan_tables.py *_profiled_metagenome.txt > ${workplace}/MetaPhlan_Results/BO_LF_MetaPhlan_merged_abundance_table_10.31.24.txt
# ^^ if this fails, go to 3a_Merge_MetaPhlan_Outputs.sh

## More info here:
# https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4
# database for Metaphlan here: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/

