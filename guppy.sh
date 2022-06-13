#!/bin/bash
#
#SBATCH -p ampere
#SBATCH -A BOEMO-SL3-GPU
#SBATCH --job-name=guppy
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/2022_05_31.guppy.stdout"
#SBATCH --error="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/2022_05_31.guppy.stderr"
#SBATCH --time=12:00:00
#! SBATCH --array=10-12

GUPPY="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/guppy_v5.0.16/ont-guppy/bin/guppy_basecaller"
INPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/ONT_data/2022_02_16_MG_ONT_AnalogueAmplicons/fast5"
OUTPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_31052022/fast5_basecalled"
CONFIG="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/guppy_v5.0.16/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg"

srun $GUPPY -i $INPUT -s $OUTPUT -x 'auto' -r -c $CONFIG --cpu_threads_per_caller 32 --barcode_kits "EXP-NBD104"
# --barcode_kits "EXP-NBD104"
