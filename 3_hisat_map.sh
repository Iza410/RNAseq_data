#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=9G
#SBATCH --time=14:00:00
#SBATCH --job-name=hs
#SBATCH --mail-user=izabela.biedron@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=logs/output_hs_%j.o
#SBATCH --error=logs/error_hs_%j.e
#SBATCH --array=0-11%6

#load module
module load UHTS/Aligner/hisat/2.2.1


#samples
samples=('HER21' 'HER22' 'HER23' 'NonTNBC1' 'NonTNBC2' 'NonTNBC3' 'Normal1' 'Normal2' 'Normal3' 'TNBC1' 'TNBC2' 'TNBC3')

#paths to folders

read_1=/data/courses/rnaseq_course/breastcancer_de/reads/${samples[$SLURM_ARRAY_TASK_ID]}_R1.fastq.gz
read_2=/data/courses/rnaseq_course/breastcancer_de/reads/${samples[$SLURM_ARRAY_TASK_ID]}_R2.fastq.gz

results=/data/users/ibiedron/rnaseq/results



#go to working directory
cd ${results}/hisat_ind

#run mapping
hisat2 -x genome_tran -1 ${read_1} -2 ${read_2} -S ${results}/mapped/${samples[$SLURM_ARRAY_TASK_ID]}.sam -p 4