#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=7-00:00:00     # 7 day
#SBATCH --output=Functional_Annotation_RawMGMs_HCs_Humann3.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Functional_Annotation_RawMGMs_HCs_Humann3"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project/ABS_MGMs/ABS_Revisions_3.3.25


if [[ ! -d ${workplace}/Humann3_Results ]]; then
    mkdir ${workplace}/Humann3_Results
    mkdir ${workplace}/Humann3_Results/output_relab
    mkdir ${workplace}/Humann3_Results/output_regrouped
    #mkdir ${workplace}/Humann3_Results/output_regrouped
    mkdir ${workplace}/Humann3_Results/output_merged
    #mkdir ${workplace}/Humann3_Results/output_merged/output_merged_relab
fi

# using trimmed, host-read-filtered sequences (aka the raw reads were first checked for quality, then trimmed, then human host sequences were removed, now ready for Humann3)

# NOTE: you need to run the following code to get the databases beyond Chocophlan and Uniref90 and other databases before you run Humann3:
 # if there are open permissions: humann_databases --download utility_mapping full $DIR
 # $DIR is where you would store the databases
 # if not open permissions, follow the example below:
 # wget https://g-227ca.190ebd.75bc.data.globus.org/humann_data/full_mapping_v201901b.tar.gz
 # tar -xvzf --update database_folders /full/path/to/contents/of/full_mapping
 
for i in ${workplace}/Kneaddata_Results/*_knead_files/*_paired_1.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_paired_1*}
    
    echo "concatenating R1 and R2 files that have been trimmed and had host reads removed"
    
    # if concatenated R1R2 file doesn't exist, create it
    if [[ ! -f ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_12.fastq ]]; then
        # concatenate R1 and R2 fastqs into single fastq file
        cat ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_1.fastq ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_2.fastq > ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_12.fastq
    
    fi
    
    # if the concatenated R1R2 file exists, run Humann3
    if [[ -f ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_12.fastq ]]; then
        # if concatenated R1R2.fastq exists, run Humann using concatenated R1 + R2 fastq file
        echo "running Humann3 on concatenated R1 and R2 mgm files"
        humann -i ${workplace}/Kneaddata_Results/${SAMPLE}_knead_files/${SAMPLE}_paired_12.fastq -o ${workplace}/Humann3_Results --threads 6 --input-format fastq --remove-temp-output --output-basename ${SAMPLE}
    fi
    
done

#normalize each output (from Humann3) to relative abundance - * see Normalization notes below!

for i in ${workplace}/Humann3_Results/*.tsv
do
        echo "Renormalizing: ${i} by relative abundance"
        f=$(basename $i)
        SAMPLE=${f%*.tsv}
                
        # relative abundance normalization
        humann_renorm_table --input $i --output ${workplace}/Humann3_Results/output_relab/${SAMPLE}_relab.tsv --units relab --update-snames
        
        # copies per million normalization aka CoPM
        # humann_renorm_table --input $i --output ${workplace}/Humann3_Results/output_CoPM/${SAMPLE}_CoPM.tsv --units cpm --update-snames
done

#regroup all genefamilies_relab files by sample first
for i in ${workplace}/Humann3_Results/output_relab/*_genefamilies_relab.tsv
do
        echo "Regrouping RelAb: ${i}"
        f=$(basename $i)
        SAMPLE=${f%_genefamilies_*}
        humann_regroup_table --input ${i} --groups uniref90_rxn --output ${workplace}/Humann3_Results/output_regrouped/${SAMPLE}_metacyc_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_go --output ${workplace}/Humann3_Results/output_regrouped/${SAMPLE}_go_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_ko --output ${workplace}/Humann3_Results/output_regrouped/${SAMPLE}_ko_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_pfam --output ${workplace}/Humann3_Results/output_regrouped/${SAMPLE}_pfam_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_eggnog --output ${workplace}/Humann3_Results/output_regrouped/${SAMPLE}_eggnog_genefamilies_relab.tsv
done

#join genefamilies, pathabundance, pathcoverage, and respective regrouped genefamilies files

### first for relative abundance data, then cpm data

## first join the *_relab files
echo "first join the pathway *_relab output"
#humann_join_tables --input ${workplace}/Humann3_Results/output_relab --output ${workplace}/Humann3_Results/output_merged/all_genefamilies_relab.tsv --file_name genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_relab --output ${workplace}/Humann3_Results/output_merged/all_pathcoverage_relab.tsv --file_name pathcoverage_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_relab --output ${workplace}/Humann3_Results/output_merged/all_pathabundance_relab.tsv --file_name pathabundance_relab

## then join the *_genefamilies_relab files
echo "then join the *_genefamilies relab files output"

humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped --output ${workplace}/Humann3_Results/output_merged/all_metacyc_genefamilies_relab.tsv --file_name metacyc_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped --output ${workplace}/Humann3_Results/output_merged/all_go_genefamilies_relab.tsv --file_name go_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped --output ${workplace}/Humann3_Results/output_merged/all_ko_genefamilies_relab.tsv --file_name ko_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped --output ${workplace}/Humann3_Results/output_merged/all_pfam_genefamilies_relab.tsv --file_name pfam_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped --output ${workplace}/Humann3_Results/output_merged/all_eggnog_genefamilies_relab.tsv --file_name eggnog_genefamilies_relab

conda deactivate

# More info on input files here: https://github.com/biobakery/humann?tab=readme-ov-file#workflow-by-input-file-type

## Databases for Humann3 - Downloading Notes - Important!
# Downloading Chocophlan database and Uniref90 and other databases for Humann3 manually (see comment here for details;
## tldr it’s a permissions issue and we have to download it this way on the university server)
# Must do the following before you use Humann:
# > wget  https://g-227ca.190ebd.75bc.data.globus.org/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
# > tar -xvzf full_chocophlan.v201901_v31.tar.gz
# > humann_config --update database_folders nucleotide /tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Humann_Resources/chocophlan

### Normalization notes
# This script provides the choice to normalize to relative abundance or copies per million (CPM) units.
## Both of these represent "total sum scaling (TSS)"-style normalization: in the former case, each sample is constrained to sum to 1, whereas in the latter case (CPMs) samples are constrained to sum to 1 million. Units out of 1 million are often more convenient for tables with many, many features (such as genefamilies.tsv tables).
## Note: CPM as used here does not refer to unnormalized COUNTS per million, but rather copies per million. CPMs as used here are a generic analog of the TPM (transcript per million) unit in RNA-seq. You may wish to use the abbreviation CoPM for added clarity.
## By default, this script normalizes all stratification levels to the sum of all community feature totals, but other options (such as level-wise normalization) are supported. "Special" features (such as UNMAPPED) can be included or excluded in the normalization process.

