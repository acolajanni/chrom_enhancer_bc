#!/usr/bin/bash


#SBATCH --mail-user antonin.colajanni@etu.u-bordeaux.fr
#SBATCH --cpus-per-task 1
#SBATCH --mem 4GB

module load macs2

save_file=/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/cutoff_diff

cat /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/bedGraph_location_score8.txt | while read path name score ; do
    echo "Path : $path"
    echo "Filename : $name"
    echo "Score cutoff : $score"
    
    sbatch --cpus-per-task 2 --mem 4GB macs2 bdgpeakcall -i $path  -o "$name.narrowPeak" --outdir $save_file -c $score 

    
    echo "__________________"

done
