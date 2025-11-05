#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G  #
#SBATCH --time=6-00:00:00     # 6 day
#SBATCH --output=Assign_Taxonomy_StrainLevel_RawMGMs_ABS_HCs_SameSTR_07.03.25.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Assign_Taxonomy_StrainLevel_RawMGMs_ABS_HCs_SameSTR"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

source /tscc/nfs/home/lfreund/miniforge3/bin/activate FindStrainsMGMs # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs
db_location=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/StrainPhlan_Resources

#if [[ ! -d ${db_location}/samestr_db ]]; then
#    samestr db --markers-info ${db_location}/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl --markers-fasta ${db_location}/mpa_vJun23_CHOCOPhlAnSGB_202403.fna.bz2 --db-version ${db_location}/mpa_latest --output-dir ${db_location}/samestr_db
#fi

if [[ ! -d ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR ]]; then
    mkdir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR
fi

# using trimmed, host-read-filtered sequences (aka the raw reads were first trimmed, then human host sequences were removed, now ready for MetaPhlan)
# filtering and trimming was performed by KneadData

samestr convert --input-files ${workplace}/MetaPhlan_Results/metaphlan_sams/*.sam.bz2 --marker-dir ${db_location}/samestr_db --tax-profiles-dir ${workplace}/MetaPhlan_Results --tax-profiles-extension _profiled_metagenome.txt --min-vcov 5 --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_convert

samestr merge --input-files ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_convert/*/*.npz --marker-dir ${db_location}/samestr_db --nprocs 6 --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_merge
    
samestr filter --input-files ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_merge/*.npz --input-names ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_merge/*.names.txt --marker-dir ${db_location}/samestr_db/ --clade-min-samples 2 --marker-trunc-len 20 --global-pos-min-n-vcov 2 --sample-pos-min-n-vcov 5 --sample-var-min-f-vcov 0.1 --nprocs 6 --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_filter/

samestr stats --input-files ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_filter/*.npz --input-names ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_filter/*.names.txt --marker-dir ${db_location}/samestr_db/ --nprocs 6 --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_stats/

samestr compare --input-files ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_filter/*.npz --input-names ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_filter/*.names.txt --marker-dir ${db_location}/samestr_db/ --nprocs 6 --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_compare/

samestr summarize --input-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_compare/ --tax-profiles-dir ${workplace}/MetaPhlan_Results --tax-profiles-extension _profiled_metagenome.txt --marker-dir ${db_location}/samestr_db --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_summarize/

#for d in ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_convert/*_L001;
#do
#    SAMPLE=$(basename $d)
#    #SAMPLE=${f%.sam.*}
#
#    samestr merge --input-files ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_convert/${SAMPLE}/*.npz --marker-dir ${db_location}/samestr_db --nprocs 6 --output-dir ${workplace}/MetaPhlan_Results/ABS_HCs_SameSTR/samestr_out_merge
#
#
#done


conda deactivate

# use metaphlan's scrit to merge metaphlan outputs into single file
# merge_metaphlan_tables.py *_profiled_metagenome.txt > ${workplace}/MetaPhlan_Results/BO_LF_MetaPhlan_merged_abundance_table_10.31.24.txt
# ^^ if this fails, go to 3a_Merge_MetaPhlan_Outputs.sh

## More info here:
# https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4
# database for Metaphlan here: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/

## SameStr functions:
#db: regenerate species marker db from MetaPhlAn markers markers-pkl
#convert: convert MetaPhlAn alignments sam to SNV profiles npy
#extract: extract SNV profiles npy from reference genomes fasta
#merge: merge SNV profiles npy + npy from multiple sources
#filter: filter SNV profiles npy
#stats: report coverage stats tsv for SNV profiles npy
#compare: compare SNV profiles npy to determine their similarity and overlapping coverage tsv
#summarize: summarize shared strains and strain co-occurrence tsv
