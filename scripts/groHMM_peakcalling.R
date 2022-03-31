################################################################################                                                                            #
# > March 2022                                                                                                                 
# > Script : groHMM peakcalling                                                                                                         
# > Function : extract peaks region of GROseq data                
# @ COLAJANNI Antonin                                                          
################################################################################

## GRO_seq
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

# I M P O R T
library(groHMM)
source("./scripts/fun_interactions.r")
data.dir = file.path(main.dir,"data/interaction/bed_file/bigwig")

GROseq_files = list.files(data.dir, pattern = "GRO")
GROseq_files = GROseq_files[grepl(pattern = ".bedGraph", GROseq_files)]

columns = c("chrom","start","end","score")

GROseq_list = list()
for (file in GROseq_files){
  path = file.path(data.dir, file)
  tmp=read.table(path, skip=1, col.names = columns)
  GROseq_list[[file]]=tmp
}

test = GROseq_list$`H1_GSE64758_GRO-seq.bedGraph`

grange = makeGRangesFromDataFrame(test)


                                  