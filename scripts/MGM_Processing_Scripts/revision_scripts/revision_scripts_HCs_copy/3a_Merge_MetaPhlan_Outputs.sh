#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G  #
#SBATCH --time=2-00:00:00     # 5 day
#SBATCH --output=Merge_MetaPhlan_Results.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Merge_MetaPhlan_Results"
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

# Input file: *_profiled_metagenome.txt; output from metaphlan

# first let's add column with sample ID so we can keep track of that
for i in ${workplace}/MetaPhlan_Results/*_profiled_metagenome.txt;
do
    f=$(basename $i)
    SAMPLE=${f%_profiled*}
    
    awk -F '\t' -v A="$SAMPLE" '{OFS="\t"} NR>4 { if ( NR == 5 ) print "clade_name", "NCBI_tax_id" ,"relative_abundance", "additional_species", "SampleID"; else print ($0,$(NF+1)=A) } ' ${workplace}/MetaPhlan_Results/${SAMPLE}_profiled_metagenome.txt > ${workplace}/MetaPhlan_Results/${SAMPLE}_profiled_mgm_labeled.txt
done

# merge all updated metaphlan outputs together into single file
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {if (!a[$0]++) print}' ${workplace}/MetaPhlan_Results/*_profiled_mgm_labeled.txt > ${workplace}/MetaPhlan_Results/ABS_HC_Proj_PRJEB32762_MetaPhlan_merged_abundance_table.txt

conda deactivate


# Understanding the output from MetaPhlan: https://github.com/biobakery/biobakery/wiki/metaphlan4#151-links-to-profile-output-files

# Code Notes
## Add sample ID column awk code:
# 'BEGIN --> tells awk we what to do when reading in the file (at beginning)
# {OFS="\t"} --> output field separator (OFS) is \t
# -v A="$OGSample" --> -v allows you to create variables in awk that you can refer to later
## ^ A is a variable we created earlier in the loop
# FNR refers to the record number (typically the line number) in the current file.
# NR refers to the total record number aka how many total lines its processed in all files read by awk
# The operator == is a comparison operator, which returns true when the two surrounding operands are equal.
## ^ Example: This means that the condition NR==FNR is normally only true for the first file, as FNR resets back to 1 for the first line of each file but NR keeps on increasing. We are now just resetting this every time so each file is treated as the first file in the loop
# NR>3 after processing third line in input file
# $0,$(NF+1)=A) --> print whole line, then add column (NF+1) and put $A there, which here is SampleID
## ^ NF --> number of columns in file
## $0 means the whole line in the file; $0=$0 and more text separated by quotes are just adding the texts and tabes to line $0

## Merge awk code:
# awk 'FNR==1 && NR!=1 --> when the line (aka record) # in current file being read is 1, but the total # of lines that have been read is not equal to 1
## FNR refers to the record number (typically the line number) in the current file.
## NR refers to the total record number aka how many total lines its processed in all files read by awk
# { while (/^<header>/) getline; } 1 --> while reading file, get header once (?)
# awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 --> when line being read in current file is 1 but the total lines read into awk is not 1, grab the header once
# {if (!a[$0]++) print} --> only print lines into new file that are NOT duplicated lines (aka lines of data that appear more than once are ignored)
# more on find command: https://math2001.github.io/article/bashs-find-command/
