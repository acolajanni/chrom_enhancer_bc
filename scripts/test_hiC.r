################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : test_hiC.r                                                        #                                                        
# > Fonction : Explore how to use hiC data                                     #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")


# I M P O R T
source("./scripts/fun_interactions.r")
load(file.path(data.dir,"3div_0.05.Rdata"))


require(HiTC)
library(GenomicRanges)
library(stringr)


H1 = File_list_filtered0.05$H1

# POLR3G chromosomic location : https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000113356;r=5:89767565-89810370
#Chromosome 5: 89767565-89810370

## Modification donnees : 


#Selection du chromosome 5 : 
chr5 = H1[grepl("chr5",H1$frag1),]
chr5$frag1 = str_remove_all(chr5$frag1, "chr5:")
chr5$frag2 = str_remove_all(chr5$frag2, "chr5:")

range_start = as.data.frame(str_split_fixed(chr5$frag1, pattern = "-", n=2))

colnames(range_start) = c("start","end")




# Autre Test :
chr5 = H1[grepl("chr5",H1$frag1),]
ranges = c("frag1","frag2")
chr5_range = chr5[colnames(chr5)%in%ranges]
test = GRangesList(chr5_range)


# Visualisation 

library(Gviz)

GeneRegionTrack(test)


#-log10(.05)
#-log10(0.001)

#
# sapply(strsplit(monTruc,":"),"[[",1) 
