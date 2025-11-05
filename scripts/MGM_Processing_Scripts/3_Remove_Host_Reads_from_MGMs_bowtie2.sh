#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G  #
#SBATCH --time=3-00:00:00     # 3 day
#SBATCH --output=Remove_HostDNA_from_MGMs_bowtie2.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Remove_HostDNA_from_MGMs_bowtie2"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project/ABS_MGMs # workplace has been updated to be specifically the ABS_MGMs directory only
bowtie2resources=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Bowtie2_Resources

if [[ ! -d ${workplace}/Bowtie2_NoHumanReads_Results ]]; then
    mkdir ${workplace}/Bowtie2_NoHumanReads_Results
fi


for i in ${workplace}/Trimmed_Seqs/*_R1_trim.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    
    # here we are using Bowtie2 to remove human host reads from the metagenomes
    ## if the host is mouse, need to use a different genome than the one included here!
    bowtie2 -p 4 -x ${bowtie2resources}/GRCh38_noalt_as/GRCh38_noalt_as -1 ${workplace}/Trimmed_Seqs/${SAMPLE}_R1_trim.fastq -2 ${workplace}/Trimmed_Seqs/${SAMPLE}_R2_trim.fastq --very-sensitive-local --un-conc ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed -S ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removal.sam
    
    samtools view -b ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removal.sam > ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removal.bam
    
    mv ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed.1 ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1.fastq
    mv ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed.2 ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R2.fastq
    
    gunzip ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R1.fastq
    gunzip ${workplace}/Bowtie2_NoHumanReads_Results/${SAMPLE}_host_removed_R2.fastq
        
done

conda deactivate

# Notes:
# Bowtie2 Manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-end-example
# Host Seq Removal tutorial: https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
## -x is to specify the basename of the index, so you have to specify the directory of the index & the basename of the *.bt2 files --> -x ./index_dirname/basename (and index_dirname and basename are the same)

# --un-conc-gz is the key command from bowtie2 (more here: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#output-options & search for -un-conc
## --un-conc-gz -- "Write paired-end reads that fail to align concordantly to file(s) at <path>. These reads correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set (depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to distinguish which file contains mate #1 and mate #2."

#workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project
#bowtie2resources=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Bowtie2Resources
#bowtie2 -p 4 -x GRCh38_noalt_as/GRCh38_noalt_as -1 O_k1-74-22_S179_L005_R1_trim.fastq -2 O_k1-74-22_S179_L005_R2_trim.fastq --very-sensitive-local --un-conc-gz O_k1-74-22_S179_L005_host_removed -S O_k1-74-22_S179_L005_host_removal.sam
