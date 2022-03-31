#!/usr/bin/bash



cat /shared/projects/chrom_enhancer_bc/data/interaction/GROseq/BigWig_paths.txt | while read path name; do
    echo "Path : $path"
    echo "Filename : $name"
    
    # default options (chr5)
    bigWigToBedGraph "$path" tmp.bedgraph -chrom=chr5 
    macs2 bdgpeakcall -i tmp.bedgraph -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/GROseq/
    rm tmp.bedgraph
    echo " "
done

echo " "
echo " "
echo " "
echo "___________________________________________________ "
echo " "
echo " "


cat /shared/projects/chrom_enhancer_bc/data/interaction/GROseq/bedgraph_paths.txt | while read path name; do
    echo "Path : $path"
    echo "Filename : $name"
    # default options (chr5)
    macs2 bdgpeakcall -i "$path" -o "$name.bedGraph" --outdir /shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/GROseq/
    echo " "
done    

