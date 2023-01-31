#!/usr/bin/env bash

#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=30G
#SBATCH --time=07:00:00
#SBATCH --job-name=sam
#SBATCH --mail-user=izabela.biedron@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=logs/output_sam_%j.o
#SBATCH --error=logs/error_sam_%j.e
#SBATCH --array=0-11%4

#load module
module add UHTS/Analysis/samtools/1.10

#for array
samples=('HER21' 'HER22' 'HER23' 'NonTNBC1' 'NonTNBC2' 'NonTNBC3' 'Normal1' 'Normal2' 'Normal3' 'TNBC1' 'TNBC2' 'TNBC3')

#paths
mapped_s=/data/users/ibiedron/rnaseq/results/mapped/${samples[$SLURM_ARRAY_TASK_ID]}.sam
bam_s=/data/users/ibiedron/rnaseq/results/mapped/${samples[$SLURM_ARRAY_TASK_ID]}.bam
sort_bam=/data/users/ibiedron/rnaseq/results/sorted/${samples[$SLURM_ARRAY_TASK_ID]}_sort.bam

cd /data/users/ibiedron/rnaseq/results/mapped

#run samtools
samtools view -hbS ${mapped_s} > ${bam_s}
samtools sort -m 30000M -@ 6 -o ${sort_bam} -T temp ${bam_s}
samtools index ${sort_bam}