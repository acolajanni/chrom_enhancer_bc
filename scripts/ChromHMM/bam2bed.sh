#!/bin/bash

################################ Slurm options #################################
### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=10:00:00

### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=12800MB

################################################################################

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
# Detail Information:
#scontrol show job $SLURM_JOB_ID
echo '########################################'

# modules loading
module load bedtools/2.26.0
module add samtools/1.9



MAIN_DIR="/shared/projects/chic_tcell_activ/data/ATAC_seq/BAM"
#BAM_DIR=$MAIN_DIR/bam
BED_DIR=/shared/projects/chic_tcell_activ/results/CHiP_seq/bed
prefix=${cond}-${donor}


#echo "samtools index"
#samtools index $BAM_DIR/$file
echo "$prefix samtools bam to bed"
bedtools bamtobed -i $MAIN_DIR/${prefix}_sorted_BL.bam > $BED_DIR/${prefix}_ATAC.bed

# What you actually want to launch
echo 'Waooouhh. Awesome.'


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
