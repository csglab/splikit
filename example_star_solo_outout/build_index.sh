#!/bin/bash
#SBATCH --job-name=build_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arshammikaeili@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=03:00:00
#SBATCH --output=/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/logs/build_index_%j.log
#SBATCH --account=rrg-hsn


GTF="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/data/mouse_ref/Mus_musculus.GRCm39.110.gtf"
FA="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/data/mouse_ref/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
OUTDIR="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/results/index/"


module load StdEnv/2020 star/2.7.9a
STAR --runThreadN 32 \
     --runMode genomeGenerate \
     --genomeDir ${OUTDIR} \
     --genomeFastaFiles ${FA} \
     --sjdbGTFfile ${GTF} \
     --sjdbOverhang 99