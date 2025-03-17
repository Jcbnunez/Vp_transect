#!/usr/bin/env bash  
#  
#SBATCH -J Index_ref  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x_%j.out    
#SBATCH -p general  

#--------------------------------------------------------------------------------
#index genome

genome=/netfiles/thermofly/Veromessor/genome/Veromessor_pergandei.UVM.Oct27.2024.fasta.softmasked

module load gcc/13.3.0-xp3epyt bwa/0.7.17-iqv3cxl picard/3.1.1-otrgwkh samtools/1.19.2-pfmpoam

bwa index $genome

picard CreateSequenceDictionary \
R=Veromessor_pergandei.UVM.Oct27.2024.fasta

samtools faidx $genome