#!/usr/bin/env bash  
#  
#SBATCH -J clean_bams  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-96


module load gcc/13.3.0-xp3epyt samtools/1.19.2-pfmpoam picard/3.1.1-otrgwkh
#samtools=/netfiles/nunezlab/Shared_Resources/Software/samtools-1.19/samtools
#PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

QUAL=40
CPU=6
JAVAMEM=18G

meta=/netfiles/thermofly/Veromessor/metadata/metadata.final.Feb24.2025.txt
working_folder=/netfiles/thermofly/Veromessor/mapping_JCBN

# Use metadata file to extract sample names and forward and reverse reads
FIL=$(cat ${meta} | awk -F '\t' '{print $8}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
SA_PT1=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
SA_PT2=$(cat ${meta} | awk -F '\t' '{print $9}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

SAMP_NAME=$(echo ${SA_PT1}_${SA_PT2})
echo $SAMP_NAME

### Veromessor BAM cleaning
### 1. SAM flag cleaning
# Filter bam files and add flags
samtools view \
-b \
-q $QUAL \
-F 0x0004 \
--threads $CPU  \
$working_folder/bams/${SAMP_NAME}.bam \
> $working_folder/bams_clean/${SAMP_NAME}.bam


# Sort with picard
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$working_folder/bams_clean/${SAMP_NAME}.bam \
O=$working_folder/bams_clean/${SAMP_NAME}.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

# Remove duplicates with picard
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
I=$working_folder/bams_clean/${SAMP_NAME}.srt.bam \
O=$working_folder/bams_clean/${SAMP_NAME}.srt.rmdp.bam \
M=$working_folder/bams_clean/${SAMP_NAME}.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# Index with samtools
samtools index $working_folder/bams_clean/${SAMP_NAME}.srt.rmdp.bam

# Do QC on cleaned bams
$qualimap bamqc \
-bam $working_folder/bams_clean/${SAMP_NAME}.srt.rmdp.bam \
-outdir $working_folder/bams_qualimap/Qualimap_${SAMP_NAME} \
--java-mem-size=$JAVAMEM

