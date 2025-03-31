#!/usr/bin/env bash  
#  
#SBATCH -J BuildDB_GATK  
#SBATCH -c 10  
#SBATCH -N 1 # on one node  
#SBATCH -t 20:00:00   
#SBATCH --mem 160G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-487

###########################################################################
#Parameters
#Java
JAVAMEM=150G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

#Load Modules
module load gatk/4.6.1.0
WORKING_FOLDER=/users/j/c/jcnunez/scratch/thermofly/basisetae/mapping

# User defined inputs -- this represents the name of the samples
intervals=Vper.intervals.txt
sample_map=files_to_haplo.txt
DB_location=db_Vero
mkdir $DB_location

#This file looks like this
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gzs

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`
echo ${i} "is being processed" $(date)

#Working folder is core folder where this pipeline is being run.
mkdir $WORKING_FOLDER/TEMP_MERGEVCF_${i}

# Linearized fasta
# awk '{ if(/^>/){ print NR==1 ? $0"\r" : "\r\n"$0"\r"}else{ printf "%s",$0}} END{print "\r"}' D.basisetae_nanopore.fasta.masked

###########################################################################
###########################################################################
# Merge VCFs using GenomicsDBImport
###########################################################################
###########################################################################

### UPDATE TO GATK4
apptainer -q exec $GATK gatk \
--java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
        GenomicsDBImport \
       --genomicsdb-workspace-path $DB_location/DB_${i} \
       --batch-size 50 \
       --sample-name-map $sample_map \
       --tmp-dir $WORKING_FOLDER/TEMP_MERGEVCF_${i} \
       --reader-threads $CPU \
       -L ${i}

rm -r $WORKING_FOLDER/TEMP_MERGEVCF_${i}

echo ${i} "done" $(date)