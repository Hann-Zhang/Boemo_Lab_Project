#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL3-CPU
#SBATCH --job-name=align
#SBATCH --nodes=1
#! (<= nodes*32)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=76
#SBATCH --output="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/2022_05_31.minimap2.barcode0%a.stdout"
#SBATCH --error="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/2022_05_31.minimap2.barcode0%a.stderr"
#SBATCH --time=6:00:00
#SBATCH --array=10-12

module load ceuadmin/samtools/1.2

GENOME="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/ONT_data/2021_10_04_MG_ONT_TAnalogue_Amplicons/referenceSequence.fasta"
TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/barcode${SLURM_ARRAY_TASK_ID}" #! ${SLURM_ARRAY_TASK_ID}"
QUERY="${TARGETDIR}/reads_barcode${SLURM_ARRAY_TASK_ID}.fastq"
MINIMAP2="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/minimap2-2.17_x64-linux/minimap2"
OUTSAM="${TARGETDIR}/alignments.sam"
OUTBAM="${TARGETDIR}/alignments.bam"
OUTPREF="${TARGETDIR}/alignments.sorted"

srun $MINIMAP2 -L -ax map-ont -t 152 -a -o $OUTSAM $GENOME $QUERY
srun samtools view -Sb -o $OUTBAM $OUTSAM
srun samtools sort -@ 152 $OUTBAM $OUTPREF
srun samtools index ${OUTPREF}.bam