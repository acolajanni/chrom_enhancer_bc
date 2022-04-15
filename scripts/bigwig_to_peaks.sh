#!/usr/bin/bash

#SBATCH --mail-user antonin.colajanni@etu.u-bordeaux.fr
#SBATCH --cpus-per-task 4
#SBATCH --mem 50GB

module load macs2


path=/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/bed_file/

cat /shared/projects/chrom_enhancer_bc/data/interaction/ATAC/BigWig_paths.txt | while read path name; do
    echo "Path : $path"
    echo "Filename : $name"
    
    # default options (chr5)
    #bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 
    #macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/chr5_default/
    #rm tmp.bedgraph
    
    # 10 Mb option (chr5)
    #bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 -start=84000000 -end=94000000 
    #macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/10mb_default/ 
    #rm tmp.bedgraph
    
    # 10 Mb option (chr5) + cutoff 1e-2
    #bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 -start=84000000 -end=94000000 
    #macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" -c 2 --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/10mb_cutoff/
    #rm tmp.bedgraph
       
    # default options (chr5) + cutoff 1e-2
    #bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5  
    #macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" -c 2 --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/chr5_cutoff/
    #rm tmp.bedgraph    
    
    # default options whole genome + cutoff 1e-2
    bigWigToBedGraph "$path" "/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/bed_file/$name.bigwig.bedgraph" 
    macs2 bdgpeakcall -i "/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/bed_file/$name.bigwig.bedgraph" -o "$name.bedGraph" -c 2 --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/whole_genome_cutoff/
    #rm tmp.bedgraph  
    
    # default options whole genome + cutoff 1e-5 (default)
    #bigWigToBedGraph "$path" tmp.bedgraph  
    macs2 bdgpeakcall -i "/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/bed_file/$name.bigwig.bedgraph" -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/whole_genome/
    #rm tmp.bedgraph  
    
    echo " "
done


