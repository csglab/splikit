#!/bin/bash
#SBATCH --job-name=STAR_solo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arshammikaeili@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=128G
#SBATCH --time=40:00:00
#SBATCH --output=/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/logs/%A_STAR_solo.log
#SBATCH --account=rrg-hsn

FILE="L8TX_181211_01_C01_S01_L003"



GENOMEDIR="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/results/index"
FASTQDIR="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/data/10x_v3"
R1=${FASTQDIR}/${FILE}/${FILE}_R1_001.fastq.gz
R2=${FASTQDIR}/${FILE}/${FILE}_R2_001.fastq.gz
CBWL="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/data/barcode_whitelist/3M-february-2018.txt"
OUTDIR="/home/arsham79/projects/rrg-hsn/arsham79/alt_splicing/results/star_solo_out"
cd ${OUTDIR}
mkdir ${FILE}

# Loging
echo "The mapping starts..."
echo ${FILE}
echo "=========================================="

# Load the required modules
module load star

# Run STARsolo
STAR --runThreadN 4 \
     --genomeDir ${GENOMEDIR} \
     --readFilesCommand zcat \
     --readFilesIn ${R2} ${R1} \
     --outFileNamePrefix ${OUTDIR}/${FILE}/ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM CY UY \
     --soloType Droplet \
     --soloCBwhitelist ${CBWL} \
     --soloUMIstart 17 \
     --soloUMIlen 12 \
     --soloBarcodeReadLength 0 \
     --clipAdapterType CellRanger4 \
     --outFilterScoreMin 30 \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
     --soloUMIfiltering MultiGeneUMI_CR \
     --soloUMIdedup 1MM_CR \
     --soloFeatures Gene GeneFull SJ Velocyto \
     --soloMultiMappers EM \
     --limitBAMsortRAM 90000000000

echo "==================== done ===================="