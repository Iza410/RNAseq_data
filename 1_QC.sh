#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=02:00:00
#SBATCH --job-name=qc
#SBATCH --mail-user=izabela.biedron@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=logs/output_qc_%j.o
#SBATCH --error=logs/error_qc_%j.e

#add module
module add UHTS/Quality_control/fastqc/0.11.7

#path to the reads
reads=/data/courses/rnaseq_course/breastcancer_de/reads
output=/data/users/ibiedron/rnaseq/results/QC
#create a folder for results
mkdir /data/users/ibiedron/rnaseq/results/QC
cd /data/users/ibiedron/rnaseq/results/QC
#run QC
fastqc ${reads}/*.fastq.gz --threads 1 -o ${output}
