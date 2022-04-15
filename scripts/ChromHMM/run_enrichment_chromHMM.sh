#!/bin/bash

################################ Slurm options #################################
### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=1:00:00

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

#echo "Enrichment G0"
#ChromHMM.sh OverlapEnrichment $MAIN_DIR/chromHMM/model${mod}${nS}states/G0_${nS}_segments.bed $MAIN_DIR/chromHMM/REGIONS $MAIN_DIR/chromHMM/model${mod}${nS}states/ATAC_enrichment_G0_${nS}states 

echo "Anchors G0 and G1"

for cl in $(ls $MAIN_DIR/chromHMM/ANCHORS | grep ATAC); do
    j=$(echo $cl | sed s/.txt//g) 
    echo $j
    ChromHMM.sh NeighborhoodEnrichment $MAIN_DIR/chromHMM/model${mod}${nS}states/G0_${nS}_segments.bed $MAIN_DIR/chromHMM/ANCHORS/$cl $MAIN_DIR/chromHMM/model${mod}${nS}states/ATAC_neibourhood_G0_${nS}states_${j}

    ChromHMM.sh NeighborhoodEnrichment $MAIN_DIR/chromHMM/model${mod}${nS}states/G1_${nS}_segments.bed $MAIN_DIR/chromHMM/ANCHORS/$cl $MAIN_DIR/chromHMM/model${mod}${nS}states/ATAC_neibourhood_G1_${nS}states_${j}
done
#echo "Enrichment G1"

#ChromHMM.sh OverlapEnrichment $MAIN_DIR/chromHMM/model${mod}${nS}states/G1_${nS}_segments.bed $MAIN_DIR/chromHMM/REGIONS $MAIN_DIR/chromHMM/model${mod}${nS}states/ATAC_enrichment_G1_${nS}states


conda deactivate
# What you actually want to launch
echo 'Waooouhh. Awesome.'


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
