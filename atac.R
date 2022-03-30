################################################################################                                                                            #
# > March 2022                                                                                                                 
# > Script : atac.R                                                                                                         
# > Function : Analyse ATACseq bigwig files                                        
# @ COLAJANNI Antonin                                                          
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/ATAC")


# I M P O R T
library(rtracklayer)
source("./scripts/fun_interactions.r")
data.dir = file.path(main.dir,"data/interaction/ATAC")
save_file = file.path(main.dir,"data/interaction/bed_file/bigwig" )

setwd(data.dir)
bigwig_atac = import.bw("./MDAMB231/GSE92898/GSM2439558_MDA-MB-231-DMSO-4w.bigwig")
bedgraph_atac = BigWig_to_bed("./MDAMB231/GSE92898/GSM2439558_MDA-MB-231-DMSO-4w.bigwig")
setwd(main.dir)


###
# R : MACSr : 3.14/3.13 : on a 3.12
#source("./docs/pkg/MACSr_1.2.0.tar.gz")
#install.packages("./docs/pkg/MACSr_1.1.0.tar.gz", repos = NULL, type = "source")
library(MACSr)

cp2 <- callpeak(bigwig_atac, store_bdg = TRUE,
                name = "run_callpeak_broad", outdir = tempdir(),
                broad = TRUE)





