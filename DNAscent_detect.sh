#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL3-CPU
#SBATCH --job-name=DNAscent
#SBATCH --nodes=1
#! (<= nodes*32)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_01062022_DNAscent/2022_06_01.DNAscent.detect.barcode0%a.stdout"
#SBATCH --error="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_01062022_DNAscent/2022_06_01.DNAscnet.detect.barcode0%a.stderr"
#SBATCH --time=6:00:00
#SBATCH --array=5-7

GENOME="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/ONT_data/2021_10_04_MG_ONT_TAnalogue_Amplicons/referenceSequence.fasta"
TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/barcode0${SLURM_ARRAY_TASK_ID}" #! ${SLURM_ARRAY_TASK_ID}"
QUERY="${TARGETDIR}/alignments.sorted.bam "
DNAscent="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/DNAscent/bin/DNAscent"
OUTPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_01062022_DNAscent/output_barcode0${SLURM_ARRAY_TASK_ID}.detect" #! ${SLURM_ARRAY_TASK_ID}"


srun $DNAscent detect -b $QUERY -r $GENOME -i /home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_01062022_DNAscent/index.dnascent -o $OUTPUT -l 500