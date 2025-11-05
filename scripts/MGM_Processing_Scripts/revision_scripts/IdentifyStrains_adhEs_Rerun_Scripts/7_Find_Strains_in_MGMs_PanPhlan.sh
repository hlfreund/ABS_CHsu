#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G  #
#SBATCH --time=6-00:00:00     # 6 day
#SBATCH --output=Find_Ecoli_Strains_in_ABS_MGMs_PanPhlan.stdout
#SBATCH --mail-user=lfreund@health.ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Find_Ecoli_Strains_in_ABS_MGMs_PanPhlan"
#SBATCH -p hotel # partition
#SBATCH --qos=hotel
#SBATCH --account=htl119

# you can use any of the following: intel, batch, highmem, gpu

module load cpu/0.17.3
module load gcc/10.2.0-2ml3m2l
#module load anaconda3

source /tscc/nfs/home/lfreund/miniforge3/bin/activate BioBakerySet # must give full path to miniforge activate function to use source, otherwise must create alias for "activate"

# set the variables
workplace=/tscc/nfs/home/lfreund/scratch/ABS_MGMs
ecoli_panphlan=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/Reference_Genomes/Ecoli_PanPhlan
panphlan_dir=/tscc/nfs/home/lfreund/ps-hsulab/ProgramResources/PanPhlan_Git/panphlan

# MUST make the output directory or PanPhlan  mapping will fail
if [[ ! -d ${workplace}/PanPhlan_Map_Results ]]; then
    mkdir ${workplace}/PanPhlan_Map_Results
fi

# run PanPhlan on each sample
for i in ${workplace}/Kneaddata_Results/*_paired_1.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_paired_1.fastq}
    
    echo ${SAMPLE} "Concatenating paired_1.fastq and paired_2.fastq"
    
    if [[ ! -f ${workplace}/Kneaddata_Results/${SAMPLE}_paired_12.fastq ]]; then
        # concatenate R1 and R2 fastqs into single fastq file
        cat ${workplace}/Kneaddata_Results/${SAMPLE}_paired_1.fastq ${workplace}/Kneaddata_Results/${SAMPLE}_paired_2.fastq > ${workplace}/Kneaddata_Results/${SAMPLE}_paired_12.fastq
        
        echo ${SAMPLE} "concatenating fastqs is done"
    
    fi
    
    echo ${SAMPLE} "PanPhlan Mapping"
    
    if [[ ! -f ${workplace}/PanPhlan_Map_Results/${SAMPLE}_panphlan_ecoli.tsv ]]; then
        
        python ${panphlan_dir}/panphlan_map.py -i ${workplace}/Kneaddata_Results/${SAMPLE}_paired_12.fastq --indexes ${ecoli_panphlan}/Escherichia_coli/Escherichia_coli -p ${ecoli_panphlan}/Escherichia_coli/Escherichia_coli_pangenome.tsv -o ${workplace}/PanPhlan_Map_Results/${SAMPLE}_panphlan_ecoli.tsv
    fi
    
    echo ${SAMPLE} "PanPhlan Mapping Complete"
    
done

echo "PanPhlan Profiling and Compilation - All Samples"

${panphlan_dir}/panphlan_profiling.py -i ${workplace}/PanPhlan_Map_Results/ --o_matrix ${workplace}/PanPhlan_Map_Results/ABS_panphlan_result_profile_ecoli.tsv --o_covmat ${workplace}/PanPhlan_Map_Results/ABS_panphlan_result_profile_ecoli_coverage.tsv --o_covplot_normed ${workplace}/PanPhlan_Map_Results/ABS_panphlan_result_profile_ecoli_coverage_plot.png --func_annot  ${ecoli_panphlan}/Escherichia_coli/panphlan_Escherichia_coli_annot.tsv -p ${ecoli_panphlan}/Escherichia_coli/Escherichia_coli_pangenome.tsv --add_ref
