#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=7-00:00:00     # 7 day
#SBATCH --output=Functional_Annotation_RawMGMs_Humann3.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Functional_Annotation_RawMGMs_Humann3"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project/ABS_MGMs


if [[ ! -d ${workplace}/Humann3_Results ]]; then
    mkdir ${workplace}/Humann3_Results
    mkdir ${workplace}/Humann3_Results/output_relab
    mkdir ${workplace}/Humann3_Results/output_CoPM
    mkdir ${workplace}/Humann3_Results/output_regrouped
    mkdir ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab
    mkdir ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM
    mkdir ${workplace}/Humann3_Results/output_merged
    mkdir ${workplace}/Humann3_Results/output_merged/output_merged_relab
    mkdir ${workplace}/Humann3_Results/output_merged/output_merged_CoPM
fi

# using trimmed, host-read-filtered sequences (aka the raw reads were first checked for quality, then trimmed, then human host sequences were removed, now ready for Humann3)

for i in ${workplace}/Bowtie2_NoHumanReads_Results/*_host_removed_R1.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_host_removed_R*}
    
    echo "concatenating R1 and R2 files that have been trimmed and had host reads removed"
    # if concatenated R1R2 file doesn't exist, create it
    if [[ ! -f ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1R2.fastq ]]; then
        # concatenate R1 and R2 fastqs into single fastq file
        cat ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1.fastq ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R2.fastq > ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1R2.fastq
    
    fi
    
    # if the concatenated R1R2 file exists, run Humann3
    if [[ -f ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1R2.fastq ]]; then
        # if concatenated R1R2.fastq exists, run Humann using concatenated R1 + R2 fastq file
        echo "running Humann3 on concatenated R1 and R2 mgm files"
        humann -i ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1R2.fastq -o ${workplace}/Humann3_Results --threads 6 --input-format fastq --remove-temp-output --output-basename ${SAMPLE}
    fi
    
done

#normalize each output to relative abundance - * see Normalization notes below!

for i in ${workplace}/Humann3_Results/*.tsv
do
        echo "Renormalizing: ${i} by relative abundance and copies per million (CoPM)"
        f=$(basename $i)
        SAMPLE=${f%*.tsv}
                
        # relative abundance normalization
        humann_renorm_table --input $i --output ${workplace}/Humann3_Results/output_relab/${SAMPLE}_relab.tsv --units relab --update-snames
        
        # copies per million normalization aka CoPM
        humann_renorm_table --input $i --output ${workplace}/Humann3_Results/output_CoPM/${SAMPLE}_CoPM.tsv --units cpm --update-snames
done

# NOTE: you need to run the following code to get the databases beyond Chocophlan and Uniref90 and other databases:
 # if there are open permissions: humann_databases --download utility_mapping full $DIR
 # $DIR is where you would store the databases
 # if not open permissions, follow the example below:
 # wget https://g-227ca.190ebd.75bc.data.globus.org/humann_data/full_mapping_v201901b.tar.gz
 # tar -xvzf --update database_folders /full/path/to/contents/of/full_mapping
 
#regroup all genefamilies_relab files by sample first
for i in ${workplace}/Humann3_Results/output_relab/*_genefamilies_relab.tsv
do
        echo "Regrouping RelAb: ${i}"
        f=$(basename $i)
        SAMPLE=${f%_genefamilies_*}
        humann_regroup_table --input ${i} --groups uniref90_rxn --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab/${SAMPLE}_metacyc_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_go --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab/${SAMPLE}_go_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_ko --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab/${SAMPLE}_ko_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_pfam --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab/${SAMPLE}_pfam_genefamilies_relab.tsv
        humann_regroup_table --input ${i} --groups uniref90_eggnog --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab/${SAMPLE}_eggnog_genefamilies_relab.tsv
done

#then regroup all genefamilies_CoPM files by sample
for i in ${workplace}/Humann3_Results/output_CoPM/*_genefamilies_CoPM.tsv
do
        echo "Regrouping CoPM: ${i}"
        f=$(basename $i)
        SAMPLE=${f%_genefamilies_*}

        humann_regroup_table --input ${i} --groups uniref90_rxn --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM/${SAMPLE}_metacyc_genefamilies_CoPM.tsv
        humann_regroup_table --input ${i} --groups uniref90_go --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM/${SAMPLE}_go_genefamilies_CoPM.tsv
        humann_regroup_table --input ${i} --groups uniref90_ko --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM/${SAMPLE}_ko_genefamilies_CoPM.tsv
        humann_regroup_table --input ${i} --groups uniref90_pfam --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM/${SAMPLE}_pfam_genefamilies_CoPM.tsv
        humann_regroup_table --input ${i} --groups uniref90_eggnog --output ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM/${SAMPLE}_eggnog_genefamilies_CoPM.tsv
done


#join genefamilies, pathabundance, pathcoverage, and respective regrouped genefamilies files

### first for relative abundance data, then cpm data

## first join the *_relab files
echo "first join the pathway *_relab output"
#humann_join_tables --input ${workplace}/Humann3_Results/output_relab --output ${workplace}/Humann3_Results/output_merged/all_genefamilies_relab.tsv --file_name genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_pathcoverage_relab.tsv --file_name pathcoverage_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_pathabundance_relab.tsv --file_name pathabundance_relab

## then join the *_genefamilies_relab files
echo "then join the *_genefamilies relab files output"
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_metacyc_genefamilies_relab.tsv --file_name metacyc_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_go_genefamilies_relab.tsv --file_name go_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_ko_genefamilies_relab.tsv --file_name ko_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_pfam_genefamilies_relab.tsv --file_name pfam_genefamilies_relab
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_relab --output ${workplace}/Humann3_Results/output_merged/output_merged_relab/all_eggnog_genefamilies_relab.tsv --file_name eggnog_genefamilies_relab

## first join the *_CoPM files
echo "then join the *_CoPM files output"
#humann_join_tables --input ${workplace}/Humann3_Results/output_CoPM --output ${workplace}/Humann3_Results/output_merged/all_genefamilies_CoPM.tsv --file_name genefamilies_CoPM
humann_join_tables --input ${workplace}/Humann3_Results/output_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_pathcoverage_CoPM.tsv --file_name pathcoverage_CoPM
humann_join_tables --input ${workplace}/Humann3_Results/output_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_pathabundance_CoPM.tsv --file_name pathabundance_CoPM

## then join the *_genefamilies_CoPM files
echo "then join the *_genefamilies_CoPM files output"
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_metacyc_genefamilies_CoPM.tsv --file_name metacyc_genefamilies_CoPM
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_go_genefamilies_CoPM.tsv --file_name go_genefamilies_CoPM
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_ko_genefamilies_CoPM.tsv --file_name ko_genefamilies_CoPM
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_pfam_genefamilies_CoPM.tsv --file_name pfam_genefamilies_CoPM
humann_join_tables --input ${workplace}/Humann3_Results/output_regrouped/output_regrouped_CoPM --output ${workplace}/Humann3_Results/output_merged/output_merged_CoPM/all_eggnog_genefamilies_CoPM.tsv --file_name eggnog_genefamilies_CoPM

conda deactivate

# More info on input files here: https://github.com/biobakery/humann?tab=readme-ov-file#workflow-by-input-file-type

## Databases for Humann3 - Downloading Notes - Important!
# Downloading Chocophlan database and Uniref90 and other databases for Humann3 manually (see comment here for details;
## tldr it’s a permissions issue and we have to download it this way on the university server)
# Must do the following before you use Humann:
# > wget  https://g-227ca.190ebd.75bc.data.globus.org/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
# > tar -xvzf full_chocophlan.v201901_v31.tar.gz
# > humann_config --update database_folders nucleotide /tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Humann_Resources/chocophlan

https://g-227ca.190ebd.75bc.data.globus.org/kneadData_databases/Homo_sapiens_hg39_T2D_Bowtie2_v0.1.tar.gz

### Normalization notes
# This script provides the choice to normalize to relative abundance or copies per million (CPM) units.
## Both of these represent "total sum scaling (TSS)"-style normalization: in the former case, each sample is constrained to sum to 1, whereas in the latter case (CPMs) samples are constrained to sum to 1 million. Units out of 1 million are often more convenient for tables with many, many features (such as genefamilies.tsv tables).
## Note: CPM as used here does not refer to unnormalized COUNTS per million, but rather copies per million. CPMs as used here are a generic analog of the TPM (transcript per million) unit in RNA-seq. You may wish to use the abbreviation CoPM for added clarity.
## By default, this script normalizes all stratification levels to the sum of all community feature totals, but other options (such as level-wise normalization) are supported. "Special" features (such as UNMAPPED) can be included or excluded in the normalization process.

