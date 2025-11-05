#!/bin/bash

# set workplace variable as the main home directory
workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs/Kneaddata_Results//Bowtie2_adhE_Results/


# add sample ID info to last column in the coverage files from bowtie2

for i in *_adhE_genes_only_bowtie2_coverage.tsv;
do
    f=$(basename $i)
    SAMPLE=${f%_adhE_genes_only_*}
    #new=${SAMPLE//\_/\.}
    echo $SAMPLE
    #awk -i inplace '(gsub(/^c/, new"_c")) {print} new="${$SAMPLE//\_/\.}"' ${i}
    #awk -v A="$SAMPLE" '{OFS="\t"} {if (NR!=1) {print $0, A}}' ${i} > ${SAMPLE}_FMT_Donor_Ecoli_Markers_Coverage_bowtie2.txt
    awk -F '\t' -v A="$SAMPLE" '{OFS="\t"} NR==1{print $0 "\t" "SampleID";} NR>1{print $0, A} ' ${i} > ${SAMPLE}_adhE_genes_only_bowtie2_coverage_labeled.txt


done

awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {if (!a[$0]++) print}' *_adhE_genes_only_bowtie2_coverage_labeled.txt > ABS_adhE_bowtie2_Coverage_Results.tsv

for i in *_adhE_genes_only_bowtie2_depth.tsv;
do
    f=$(basename $i)
    SAMPLE=${f%_adhE_genes_only_*}
    #new=${SAMPLE//\_/\.}
    echo $SAMPLE
    #awk -i inplace '(gsub(/^c/, new"_c")) {print} new="${$SAMPLE//\_/\.}"' ${i}
    #awk -v A="$SAMPLE" '{OFS="\t"} {if (NR!=1) {print $0, A}}' ${i} > ${SAMPLE}_FMT_Donor_Ecoli_Markers_Coverage_bowtie2.txt
    awk -F '\t' -v A="$SAMPLE" '{OFS="\t"} NR==1{print "Chrom" "\t" "Pos" "\t" "ReadCount" "\t" "SampleID";} NR>1{print $0, A} ' ${i} > ${SAMPLE}_adhE_genes_only_bowtie2_depth_labeled.txt
    
done

awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {if (!a[$0]++) print}' *_adhE_genes_only_bowtie2_depth_labeled.txt > ABS_adhE_bowtie2_Depth_Results.tsv
