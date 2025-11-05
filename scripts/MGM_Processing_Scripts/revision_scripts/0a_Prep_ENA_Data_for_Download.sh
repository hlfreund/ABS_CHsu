#!/bin/sh

# first, pull out accession IDs from csv file
cut -d "," -f 14 seqIDmatched_clean.csv | tr -d '"' > accessions_list.txt

# there are hidden ^M at the end of each line now in accessions_list.txt, and we can get rid of them with vim:
# run vim, then enter each line
> vim accessions_list.txt
> :e ++ff=dos
> :set ff=unix
> :wq
# ^ this will remove ^M at the end of each line in accessions_list.txt

# then to pull out all the correct wget lines that match accessions in accessions_list.txt, run this loop:

for i in $(cat accessions_list.txt); do
    grep -e $i ena-file-download-read_run-PRJEB32762-fastq_ftp-20250304-0520.sh >> test.sh
done
