#!/bin/bash

################################ Slurm options #################################
### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=10:00:00

### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=12500MB

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
#module load python/2.7
source ~/.bashrc

module load conda
conda activate env_chromhmm

MAIN_DIR="/shared/projects/chic_tcell_activ/results"
binar_DIR=$MAIN_DIR/chromHMM/binar/binarized_histones
BED_DIR=$MAIN_DIR/CHiP_seq/bed

#echo "BinarizeBed histones"
#ChromHMM.sh BinarizeBed -b 200 /shared/home/edarbo/chic_tcell_activ/data/hg19.chrom.sizes $BED_DIR $MAIN_DIR/CHiP_seq/cellmarkfiletable_noATAC.txt $binar_DIR

#rm $binar_DIR/*gl*
#rm $binar_DIR/*hap*

echo "Learn model"
ChromHMM.sh LearnModel -p 1 $binar_DIR $MAIN_DIR/chromHMM/model_noATAC_${nS}states ${nS} hg19  


conda deactivate
# What you actually want to launch
echo 'Waooouhh. Awesome.'


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
