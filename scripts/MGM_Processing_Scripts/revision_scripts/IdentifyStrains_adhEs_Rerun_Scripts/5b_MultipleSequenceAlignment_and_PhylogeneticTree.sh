#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G  #
#SBATCH --time=2-00:00:00     # 2 day
#SBATCH --output=Build_16SV3V4_PhyloTree.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Build_16SV3V4_PhyloTree"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

#source activate /tscc/nfs/home/lfreund/miniforge3/envs/MGM_Analysis
source /tscc/nfs/home/lfreund/miniforge3/bin/activate PhyloTreeConda

# workplace will be the directory where your sequence data are stored
workplace=/tscc/nfs/home/lfreund/ps-hsulab/data/FT_16S_V3V4

# You must have generated a clean, updated FASTA file only containing the ASVs you are studying, aka no contaminant or eukaryotic ASVs
## if you do not have this yet, go back and run 5a_Prep_Data_for_MSA_PhylogeneticTree.R

# First load all modules you will need to generate multiple sequence alignment (MSA) and build phylogenetic tree!
module load ssu-align hmmer fasttree iqtree

# Step1: run the multiple sequence alignment of 16S V3V4 rRNA with ssu-align
ssu-align ${workplace}/DADA2_Results/FT_16SV3V4_ASVs_Updated.fa FT_16SV3V4_MSA_ssualign

# Step 2: convert Stockholm aligned file (.stk) into a multifasta using esl-reformat which is part of the easel toolkit within hmmer
esl-reformat afa ${workplace}/FT_16SV3V4_MSA_ssualign/FT_16SV3V4_MSA_ssualign.bacteria.stk > ${workplace}/FT_16SV3V4_MSA_ssualign/FT_16SV3V4_MSA_ssualign.bacteria.fasaln

# Step 3: run FastTree to build phylogenetic tree from 16S V3V4 MSA
FastTreeMP -nt -gtr -gamma < ${workplace}/FT_16SV3V4_MSA_ssualign/FT_16SV3V4_MSA_ssualign.bacteria.fasaln > ${workplace}/FT_16SV3V4_MSA_ssualign/FT_16SV3V4_MSA_ssualign.bacteria.bacteria.FastTree.tre

# Step 4 (Optional) : run iqtree to build phylogenetic tree
iqtree2 -s ${workplace}/FT_16SV3V4_MSA_ssualign/FT_16SV3V4_MSA_ssualign.bacteria.fasaln -nt AUTO -B 1000 --alrt 1000

# After you have finished running this, you can now generate Unifrac distance!
