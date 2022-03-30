#!/usr/bin/bash



cat /shared/projects/chrom_enhancer_bc/data/interaction/ATAC/BigWig_paths.txt | while read path name; do
    echo "Path : $path"
    echo "Filename : $name"
    
    # default options (chr5)
    bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 
    macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/chr5_default/
    rm tmp.bedgraph
    
    # 10 Mb option (chr5)
    bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 -start=84000000 -end=94000000 
    macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/10mb_default/ 
    rm tmp.bedgraph
    
    # 10 Mb option (chr5) + cutoff 1e-2
    bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 -start=84000000 -end=94000000 
    macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" -c 2 --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/10mb_cutoff/
    rm tmp.bedgraph
       
    # default options (chr5) + cutoff 1e-2
    bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5  
    macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" -c 2 --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/chr5_cutoff/
    rm tmp.bedgraph    
    
    
    echo " "
done


