#!/usr/bin/env bash  
#  
#SBATCH -J BuildDB_GATK  
#SBATCH -c 5  
#SBATCH -N 1 # on one node  
#SBATCH -t 6:00:00   
#SBATCH --mem 50G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-487


###########################################################################
#Parameters
#Java
JAVAMEM=40G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

#Load Modules
module load gatk/4.6.1.0
WORKING_FOLDER=/netfiles/thermofly/Veromessor/mapping_JCBN
REFERENCE=Veromessor_pergandei.UVM.Oct27.2024.softmasked.fasta
intervals=Vper.intervals.txt

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

# Identify the Genome database to genotyoe
GenomeDB_path=`echo $WORKING_FOLDER/db_Vero/DB_${i}`

echo "now processing DB" ${i} $(date)


####
###########################################################################
###########################################################################
# Genotype call the samples in the DBI merged set
###########################################################################
###########################################################################

mkdir Genotyped_vcfs

apptainer -q exec $GATK gatk \
--java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
    GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GenomeDB_path \
   -O $WORKING_FOLDER/Genotyped_vcfs/${i}.genotyped.raw.vcf.gz

echo ${i} "done" $(date)
