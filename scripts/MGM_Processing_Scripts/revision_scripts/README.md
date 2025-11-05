# Revision Scripts - ABS MGMs

Data used for these revisions came from https://www.ebi.ac.uk/ena/browser/view/PRJEB32762
Project PRJEB32762
Include healthy controls that were selected via propensity matching by Dr. Cynthia Hsu

The order to run these scripts is provided by their #s, starting with 0a

The workflow:

0. download sequence data from ENA
-- get download script (ena-file-downlad-read_run-*.sh)
-- 0a. pull out accessions for samples from metadata csv (0a_Prep_ENA_Data_for_Download.sh)
-- 0b. download the specific samples by using accession #s to pull out specific wget lilnes from ena-file-downlad-read_run-*.sh (0b_Download_Specific_ENA_Files.sh)

1. Check raw metagenome sequence data with FastQC & MultiQC

2. Filter and trim metagenome reads with KneadData
-- 2a. Check quality of filtered, trimmed metagenome seq data with FastQC & MultiQC

3. Assign taxonomy to unassembled metagenomes with MetaPhlan4
-- 3a. Merge MetaPhlan4 results into one large txt file (tab delimited)

4. Functional annotation of unassembled metagenomes with Humann3

Run by Linton Freund 3/4-5/2025

Update 6/27/2025
- Ran StrainPhlan on ABS metagenomes (excluding the propensity matched HCs) to identify Ecoli and Kpneumoniae strains in the mgms, with particular interest in the FMT patient + Donor samples
- StrainPhlan created tree based on clade markers found in sample, but the tree doesn't include all the samples
- also ran SameStr on mgms to find Ecoli strains in the samples, and was able to get table describing 

Update 7/10/2025
- tried several methods to identify Ecoli and Kpneumoniae strains in ABS mgms (i.e., Strainphlan, SameStr)
- mapped ABS MGMs to Ecoli genome (K-12 strain), but only a single genome -- didn't really yield useful results
- created database of bacterial adhE genes from various taxa (taxa identified by Humann as having adhE KO ID or randomly selected taxa), mapped ABS mgm reads to these genes with bowtie2, used samtools to calculate coverage + stats
- ^^ will do the same thing with propsensity-matched HCs used for ABS manuscript (project ID PRJEB32762)
