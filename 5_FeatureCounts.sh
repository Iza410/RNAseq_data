#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=15:00:00
#SBATCH --job-name=FC
#SBATCH --mail-user=izabela.biedron@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=logs/output_FC_%j.o
#SBATCH --error=logs/error_FC_%j.e


#load the module
module load UHTS/Analysis/subread/2.0.1

anno=/data/users/ibiedron/rnaseq/ref/Homo_sapiens.GRCh38.108.gtf
out=/data/users/ibiedron/rnaseq/results/featureCounts/featureCounts.txt
sorted_bam=/data/users/ibiedron/rnaseq/results/sorted/*_sort.bam


module add UHTS/Analysis/subread/2.0.1

featureCounts -p -C -s 0 -a $anno -o $out $sorted_bam -T 8 -Q 10 --tmpDir $SCRATCH -g gene_id