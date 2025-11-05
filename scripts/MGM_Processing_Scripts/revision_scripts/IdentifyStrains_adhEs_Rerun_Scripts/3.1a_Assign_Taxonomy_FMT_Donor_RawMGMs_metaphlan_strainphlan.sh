#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=6-00:00:00     # 5 day
#SBATCH --output=Assign_Taxonomy_RawMGMs_ABS_and_HCs_strainphlan_FMT.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Assign_Taxonomy_RawMGMs_ABS_and_HCs_strainphlan_FMT_Only"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs

# must have Python module loaded for these scripts to work!

module load python # need this for strainphlan .py functions

# Extract strain-specific markers (Ecoli) from MetaPhlan db
#extract_markers.py -c t__SGB10068 -o ${workplace}/MetaPhlan_Results/db_markers/

# Extract strain-specific markers (Klebsiella pneumoniae) from MetaPhlan db
#extract_markers.py -c t__SGB10115_group -o ${workplace}/MetaPhlan_Results/db_markers/

if [[ ! -d ${workplace}/MetaPhlan_Results/StrainPhlan_Results_FMT2/ ]]; then
    mkdir ${workplace}/MetaPhlan_Results/StrainPhlan_Results_FMT2

fi

# StrainPhlan uses all samples to pull out clade markers to generate MSA and phylo tree, tracing each sample in the process

echo "Running StrainPhlan on FMT samples & FMT-Donor sample only"

# first find Ecoli clade markers
strainphlan -s ${workplace}/MetaPhlan_Results/consensus_markers/FMT_Patient_Only/*.json.bz2 -m ${workplace}/MetaPhlan_Results/db_markers/t__SGB10068.fna -o ${workplace}/MetaPhlan_Results/StrainPhlan_Results_FMT2 -n 6 -c t__SGB10068 --sample_with_n_markers 10 --sample_with_n_markers_perc 10 --marker_in_n_samples_perc 25 --sample_with_n_markers_after_filt 10 --sample_with_n_markers_after_filt_perc 10 --mutation_rates

# then find Klebsiella pneumoniae clade markers
#strainphlan -s ${workplace}/MetaPhlan_Results/consensus_markers/FMT_Patient_Only/*.json.bz2 -m ${workplace}/MetaPhlan_Results/db_markers/t__SGB10115_group.fna -o ${workplace}/MetaPhlan_Results/StrainPhlan_Results_FMT -n 6 -c t__SGB10115_group --mutation_rates

#echo "Running StrainPhlan " > ${workplace}/MetaPhlan_Results/StrainPhlan_Results_FMT/README.md

conda deactivate

# use metaphlan's scrit to merge metaphlan outputs into single file
# merge_metaphlan_tables.py *_profiled_metagenome.txt > ${workplace}/MetaPhlan_Results/BO_LF_MetaPhlan_merged_abundance_table_10.31.24.txt
# ^^ if this fails, go to 3a_Merge_MetaPhlan_Outputs.sh

## More info here:
# https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4
# database for Metaphlan here: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/

