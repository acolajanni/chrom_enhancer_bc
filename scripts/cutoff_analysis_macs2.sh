#!/usr/bin/bash


#SBATCH --mail-user antonin.colajanni@etu.u-bordeaux.fr
#SBATCH --cpus-per-task 6
#SBATCH --mem 50GB

module load macs2

save_file=/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/cutoff

cat /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/bedGraph_location.txt | while read path name; do
    echo "Path : $path"
    echo "Filename : $name"
    
    macs2 bdgpeakcall -i $path  -o "$name.narrowPeak" --outdir $save_file --cutoff-analysis
    
    echo "__________________"

done
