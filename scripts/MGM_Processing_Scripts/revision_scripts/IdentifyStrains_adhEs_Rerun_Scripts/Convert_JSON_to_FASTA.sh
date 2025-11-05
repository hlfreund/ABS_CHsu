#!/bin/bash

# set workplace variable as the main home directory
workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs/MetaPhlan_Results

# using json files containing consensus marker sequences from StrainPhlan
## need to run up to & including Step 3 in this tutorial first before running this script: https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1#usage
## must also have paste installed!

for FILE in ${workplace}/consensus_markers/FMT_Patient_Only/*.json;
do
    f=$(basename $FILE)
    SAMPLE=${f%_L00*.json}
    
    # grep "marker" from JSON file, cut the second field (aka containing the Uniref ID), remove the extra characters and create a new file
    ## new file will contain the sequence ID with a > in front so it's in fasta format (e.g. >UniRef90_A0A354L2I3|1__3|SGB15286)
    
    grep 'marker' ${f}  | cut -f2 -d : | sed 's/ "/>/g' | sed 's/",//g' | tr -d ' [ ' | sed '1{/^$/d}' > ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_marker_IDs.txt
    
    # grep "sequence" from JSON file, cut the second field which contains the DNA sequence, remove extra characters, and store that in a second txt file
    ## this will have a new sequence on each line in the file
    
    grep 'sequence' ${f} | cut -f2 -d : | tr -d '" ' > ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_marker_seqs.txt
    

    if [[ -f ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_marker_IDs.txt ]] && [[ -f ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_marker_seqs.txt ]]; then
        # if these two previous files exist, then use paste to interleave the two files, using the -d '\n' to indicate that the new line is the delimiter for merging
        ## paste will take line 1 from file 1, then add line 1 from file 2 on the next line, then interleave the files to create a FASTA file
        ## this step means you have to have paste installed!
    
        paste -d '\n' ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_marker_IDs.txt ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_marker_seqs.txt > ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_consensus_markers_FASTA.txt
            

    fi
    
    
done

# now to pull out just the Ecoli clade markers from new consensus_markers_FASTA.txt

for FILE in ${workplace}/consensus_markers/FMT_Patient_Only/*_consensus_markers_FASTA.txt;
do
    f=$(basename $FILE)
    SAMPLE=${f%_consensus_*.txt}
   
    # sed -n suppresses printing the pattern (i.e., normal sed output)
    ## '/SGB10068/{N;p;}' means we are looking for pattern SGB10068; when you find pattern, print that line & the next line
    sed -n '/SGB10068/{N;p;}' ${FILE} > ${workplace}/consensus_markers/FMT_Patient_Only/${SAMPLE}_Ecoli_consensus_markers_FASTA.txt
    
done
       
# test code below
#sed -n '/SGB10068/{N;p;}' PA_11_10_22_S74_consensus_markers_FASTA.txt > test.txt 
