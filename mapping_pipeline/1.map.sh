#!/usr/bin/env bash  
#  
#SBATCH -J Map_reads  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=X-Y

# ----------------------------------------------

genome=/netfiles/thermofly/Veromessor/genome/Veromessor_pergandei.UVM.Oct27.2024.fasta.softmasked
module load gcc/13.3.0-xp3epyt bwa/0.7.17-iqv3cxl samtools/1.19.2-pfmpoam

#--------------------------------------------------------------------------------

echo ${SLURM_ARRAY_TASK_ID}

# Set folders and file locations
working_folder=/netfiles/thermofly/Veromessor/mapping
meta=/netfiles/thermofly/Veromessor/metadata/metadata.final.Feb24.2025.txt
reads_folder=/netfiles/thermofly/Veromessor/Vp_clean_reads

# Use metadata file to extract sample names and forward and reverse reads
FIL=$(cat ${meta} | awk -F '\t' '{print $8}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
SA_PT1=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
SA_PT2=$(cat ${meta} | awk -F '\t' '{print $9}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

SAMP_NAME=$(echo ${SA_PT1}_${SA_PT2})

echo $FIL
echo $SAMP_NAME

#### Where to map

outfile=/netfiles/thermofly/Veromessor/mapping/jcbn_map
mkdir $outfile

#--------------------------------------------------------------------------------
# Parameters for software
CPU=6
QUAL=40 
#--------------------------------------------------------------------------------

# Map reads to reference
echo "Map" ${SAMP_NAME}

# Map with bwa mem2
bwa mem -M -t $CPU $genome \
$reads_folder/$FIL \
> $outfile/${SAMP_NAME}.sam

# Build bam files
samtools view -b -q $QUAL --threads $CPU  \
$outfile/${SAMP_NAME}.sam \
> $outfile/${SAMP_NAME}.bam

rm $outfile/${SAMP_NAME}.sam
