#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=160G
#SBATCH --time=04:00:00
#SBATCH --job-name=hs
#SBATCH --mail-user=izabela.biedron@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=logs/output_hs_%j.o
#SBATCH --error=logs/error_hs_%j.e

#load module
module load UHTS/Aligner/hisat/2.2.1

#paths to folders
reads=/data/courses/rnaseq_course/breastcancer_de/reads
ref=/data/users/ibiedron/rnaseq/ref
results=/data/users/ibiedron/rnaseq/results



#create a folder for results
cd ${results}
mkdir hisat_ind
cd hisat_ind

hisat2_extract_splice_sites.py ${ref}/Homo_sapiens.GRCh38.108.gtf > ${results}/hisat_ind/genome.ss
hisat2_extract_exons.py ${ref}/Homo_sapiens.GRCh38.108.gtf > ${results}/hisat_ind/genome.exon

hisat2-build -p 16 --exon ${results}/hisat_ind/genome.exon \
--ss ${results}/hisat_ind/genome.ss ${ref}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
${results}/hisat_ind/genome_tran
