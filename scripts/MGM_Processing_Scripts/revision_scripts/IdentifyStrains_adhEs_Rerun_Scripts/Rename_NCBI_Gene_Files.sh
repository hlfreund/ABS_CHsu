#!/bin/bash

# rename fna file by directory name we made when we downloaded from NCBI
for i in $(ls -d ./*_datasets); do
    echo ${i}
    f=$(basename $i)
    genome_name=${f%_datasets}
    mv ./${i}/ncbi_dataset/data/gene.fna ./${i}/ncbi_dataset/data/${genome_name}_gene.fna
done

# same code as above but in one line
for i in $(ls -d ./*_datasets); do f=$(basename $i); genome_name=${f%_datasets}; mv ./${i}/ncbi_dataset/data/gene.fna ./${i}/ncbi_dataset/data/${genome_name}_gene.fna; done

# code to move files from sub directories to current directory
for i in $(ls -d ./*_datasets/ncbi_dataset/data/*_gene.fna); do mv $i ./; done

# rename seq headers in fasta files
for i in *.fna; do sed -i 's/>[A-Z0-9]*[^[:space:]]*[[:space:]]/>/g' $i; done

# for a Mac: sed requires file extension included; see format below
for i in *.fna; do sed -i'.fna' -e 's/>[A-Z0-9]*[^[:space:]]*[[:space:]]/>/g' $i; done

# remove organism label from fna
sed -i'.fna' -e 's/organism\=//g' bacteria_adhE_gene_seqs.fna

# give each adhE header a unique identifier so it can be properly indexed by bowtie2
awk '/^>/{sub(/^>/,">Seq"++i"_");}1' bacteria_adhE_gene_seqs.fna > bacteria_adhE_gene_labeled_seqs.fna
