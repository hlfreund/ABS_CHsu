#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G  #
#SBATCH --time=3-00:00:00     # 3 day
#SBATCH --output=Trim_MGM_Seqs_bbduk.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Trim MGM Sequences bbduk"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate RawMetagenomeWorkflow

workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/ABS_Project
bbduksources=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/BBdukResources

if [[ ! -d ${workplace}/ABS_MGMs/Trimmed_Seqs ]]; then
    mkdir ${workplace}/ABS_MGMs/Trimmed_Seqs
fi

# first, gunzip the FASTQ files
gunzip ${workplace}/ABS_MGMs/*.fastq.gz

# then let's trim and decontaminate!
for i in ${workplace}/ABS_MGMs/*_R1_001.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    
    if [[ ! -f ${workplace}/ABS_MGMs/Trimmed_Seqs/${SAMPLE}_R1_trim.fastq ]] && [[ ! -f ${workplace}/ABS_MGMs/Trimmed_Seqs/${SAMPLE}_R2_trim.fastq ]]; then
        bbduk.sh -Xmx10g in1=${workplace}/ABS_MGMs/${SAMPLE}_R1_001.fastq in2=${workplace}/ABS_MGMs/${SAMPLE}_R2_001.fastq out1=${workplace}/ABS_MGMs/Trimmed_Seqs/${SAMPLE}_R1_trim.fastq out2=${workplace}/ABS_MGMs/Trimmed_Seqs/${SAMPLE}_R2_trim.fastq  ref=${bbduksources}/adapters.fa ftl=10 ftr=147 rcomp=t ktrim=r k=23 maq=10 minlength=120 trimpolyg=6 trimpolya=6 mink=11 hdist=1 tpe tbo
    fi
    
done

conda deactivate

# Adapter Sequences identified by FastQC can be found here: https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/AdapterSequencesIntro.htm;
## http://docs.blast2go.com/user-manual/tools-(pro-feature)/fastq-quality-check/
## use literal= flag to remove these adapter seqs ^^

#/Volumes/HLF_SSD/Aronson_Lab_Data/bbmap/bbduk.sh in1=/Volumes/HLF_SSD/Aronson_Lab_Data/SaltonSea_POW/EA_Pool-POW_1-1a_S28_L001_R1_001_clean1.fastq in2=EA_Pool-POW_1-1a_S28_L001_R2_001_clean1.fastq out1=/Volumes/HLF_SSD/Aronson_Lab_Data/SaltonSea_POW/EA_Pool-POW_1-1a_S28_L001_R1_001_clean.fq out2=/Volumes/HLF_SSD/Aronson_Lab_Data/SaltonSea_POW/EA_Pool-POW_1-1a_S28_L001_R2_001_clean.fq  literal=ACTGCGAA,TTCGCAGT rcomp=t ktrim=r k=23 mink=11 hdist=1 maq=10 minlen=51 trimq=20 tpe tbo

# More info here: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
# also http://seqanswers.com/forums/showthread.php?t=42776
# also also: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk.sh
# ref ---> file provided by bbduk that holds collection of Illumina TruSeq adapters
# literal=(sequence here) ---> literal adapter sequences to remove; "N" represents any base -- in this case, they are indexes within the adapters
# rcomp=t ---> Rcomp looks for kmers and their reverse-complements, rather than just forward kmer, if set to true
# ktrim=r ---> “ktrim=r” is for right-trimming (3′ adapters);
## In ktrim=r mode, once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; this is the normal mode for adapter trimming.
# k=23 ---> look for kmer that is 23 bp long
# trimpolyg=6 --> trim poly G sequences by 6 bps; can also trim polyA or polyT tails
# mink=11 ---> in addition to kmers of x length, look for shorter kmers with lengths 23 to 11 (in this case)
# maq=10 ---> This will discard reads with average quality below 10
# hdist=1 ---> hamming distance of 1 (for identifying kmers)
# mlf=50 ---> (minlengthfraction=50) would discard reads under 50% of their original length after trimming
# trimq=10 ---> quality-trim to Q10 using the Phred algorithm, which is more accurate than naive trimming.
# qtrim=r ---> means it will quality trim the right side only [happens after all base kmer operations]
# tpe ---> which specifies to trim both reads to the same length
# tbo ---> which specifies to also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences)
# mm ----> Maskmiddle ignores the middle base of a kmer, can turn off with mm=f
# ftl ----> (force) trim the left most 5 bases

# more on hdist...
# A 4.6Mbp genome like E.coli contains around 4.6 million unique kmers. If a hamming distance is used, such as hdist=1, then the number of kmers stored will be
# multiplied by 1+(3*k)^hdist. So, for E.coli with K=31 and hdist=0, there are 4554207 kmers stored, using around 140MB, taking about 0.5 seconds; with
# hdist=1, there are 427998710 kmers stored (94 times as many), using 15GB, and taking 104 seconds
# BBDuk needs around 20 bytes per kmer

# BBDuk supports kmers of length 1-31. The longer a kmer, the high the specificity
# Note that it is impossible to do kmer matching for reference sequences that are shorter than K.
# When kmer filtering, you can specify a kmer length greater than 31 (such as k=40) for even higher specificity.
# For k=40, BBDuk will consider a read to match the reference if there are 10 consecutive 31-mer matches. This is not exactly the same as actually matching 40-mers, but is very similar.
# Example Usage
## bbduk.sh -Xmx10g in1=${workplace}/16S_Seqs/${SAMPLE}_R1.fastq in2=${workplace}/16S_Seqs/${SAMPLE}_R2.fastq out1=${workplace}/Trimmed_Seqs/16S_Trimmed/${SAMPLE}_R1_clean.fastq out2=${workplace}/Trimmed_Seqs/16S_Trimmed/${SAMPLE}_R2_clean.fastq literal=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ref=/bigdata/aronsonlab/shared/bbmap_resources/adapters.fa ftl=5 ftr=265 rcomp=t ktrim=r k=23 maq=10 minlength=250 mink=11 hdist=1 tpe tbo

