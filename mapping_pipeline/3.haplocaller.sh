#!/usr/bin/env bash  
#  
#SBATCH -J HaploCall  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 20:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-96

module load gatk/4.6.1.0
module load gcc/13.3.0-xp3epyt \
picard/3.1.1-otrgwkh \
htslib/1.19.1-6ivqauw \
picard/3.1.1-otrgwkh \
samtools/1.19.2-pfmpoam

#### Locations
BAMS_FOLDER=/netfiles/thermofly/Veromessor/mapping_JCBN/bams_clean
WORKING_FOLDER=/netfiles/thermofly/Veromessor/mapping_JCBN
REFERENCE=Veromessor_pergandei.UVM.Oct27.2024.softmasked.fasta

#### Metadata
meta=/netfiles/thermofly/Veromessor/metadata/metadata.final.Feb24.2025.txt
SUFFIX=srt.rmdp

#### Java info
JAVAMEM=18G
CPU=$SLURM_CPUS_ON_NODE
echo ${CPU}


#Read Information
Group_library="Vero_cline"
Library_Platform="G4"
Group_platform="Veromessor"

#HaploCaller
HET=0.005

#####
FIL=$(cat ${meta} | awk -F '\t' '{print $8}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
SA_PT1=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
SA_PT2=$(cat ${meta} | awk -F '\t' '{print $9}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

i=$(echo ${SA_PT1}_${SA_PT2})
echo $i

###########################################################################
###########################################################################
# Forcing a uniform read group to the joint bam file
###########################################################################
###########################################################################

mkdir $WORKING_FOLDER/RGSM_final_bams

picard AddOrReplaceReadGroups \
  I=$BAMS_FOLDER/${i}.$SUFFIX.bam \
  O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
  RGLB=$Group_library \
  RGPL=$Library_Platform \
  RGPU=$Group_platform \
  RGSM=${i}


###########################################################################
###########################################################################
# Index Bam files
###########################################################################
###########################################################################

picard BuildBamIndex \
      I=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
      O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bai

###########################################################################
###########################################################################
# Haplotype Calling
###########################################################################
###########################################################################

mkdir $WORKING_FOLDER/haplotype_calling

###NEED TO MAKE DICTIONARY
### CreateSequenceDictionary is finicky and needs the file to end on ".fasta"
#picard CreateSequenceDictionary \
#      R=Veromessor_pergandei.UVM.Oct27.2024.softmasked.fasta \
#      O=Veromessor_pergandei.UVM.Oct27.2024.softmasked.fasta.dict
#samtools faidx Veromessor_pergandei.UVM.Oct27.2024.softmasked.fasta

### GATK is also finicky and needs the file to end on ".fasta"

$GATK \
gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
HaplotypeCaller \
  -R Veromessor_pergandei.UVM.Oct27.2024.softmasked.fasta \
  -I $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
  -O $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf \
  --heterozygosity $HET \
  -ploidy 2 \
  -ERC GVCF
  
  
###########################################################################
###########################################################################
# Compress and index with Tabix
###########################################################################
###########################################################################

bgzip $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf
tabix $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf.gz

echo "done"