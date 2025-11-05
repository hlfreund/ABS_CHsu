#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G  #
#SBATCH --time=5-00:00:00     # 5 day
#SBATCH --output=Filter_for_Ecoli_Reads_EcoliMarkers_Kneaddata.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Filter_for_Ecoli_Reads_EcoliMarkers_Kneaddata"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

## must have prepared database for KneadData first -- see notes below @ bottom of script!

workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs/Kneaddata_Results/
knead_db=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Kneaddata_Resources
ecoli_ref=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Reference_Genomes/Ecoli_Ref_Genome
ecoli_markers=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Reference_Genomes/Ecoli_Ref_Genome/Ecoli_MarkerGenes

# we are using human-read-removed FMT patient + Donor samples only as input to pull out only Ecoli reference reads
## we used KneadData to remove human-host reads, now we will run KneadData again to pull out only reads that map to Ecoli reference

if [[ ! -d ${workplace}/FMT_Donor_Seqs_Only/Kneaddata_EcoliMarkers ]]; then
    mkdir ${workplace}/FMT_Donor_Seqs_Only/Kneaddata_EcoliMarkers

fi

# Kneaddata can trim the reads and remove host genome reads before taxonomic and functional annotation

for i in ${workplace}/FMT_Donor_Seqs_Only/*_1.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_paired_*}
    
    # if results dir called ${SAMPLE}_knead_files does not exist, then run kneaddata
    if [[ ! -d ${workplace}/FMT_Donor_Seqs_Only/Kneaddata_EcoliMarkers/${SAMPLE}_ecoli_markers_knead ]]; then
        
        kneaddata --input1 ${workplace}/FMT_Donor_Seqs_Only/${SAMPLE}_paired_1.fastq --input2 ${workplace}/FMT_Donor_Seqs_Only/${SAMPLE}_paired_2.fastq -db ${ecoli_markers} -t 4 --reorder --store-temp-output --output-prefix ${SAMPLE}_ecoli --output ${workplace}/FMT_Donor_Seqs_Only/Kneaddata_EcoliMarkers/${SAMPLE}_ecoli_markers_knead
        
        mv ${workplace}/FMT_Donor_Seqs_Only/Kneaddata_EcoliMarkers/${SAMPLE}_ecoli_markers_knead/${SAMPLE}_*_contam.fastq ${workplace}/FMT_Donor_Seqs_Only/Kneaddata_EcoliMarkers
    
    fi
   
    
done

# NOTE -- because we are looking for sequences that are found in the Ecoli reference genome, we want to use the files that have _paired_contam_1.fastq or _paired_contam_2.fastq for downstream analyses!
## contam means that they were found in the reference database we used!

conda deactivate

# More info on KneadData: https://github.com/biobakery/kneaddata?tab=readme-ov-file#paired-end-run
# output information also here: https://github.com/biobakery/biobakery/wiki/kneaddata#paired-end-reads

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
