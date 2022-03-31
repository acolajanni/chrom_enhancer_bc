################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : get_BENGI                                                         #                                                        
# > Function : get pcHiC files of BENGI datasets                               #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/pcHiC/BENGI/Benchmark/")


# I M P O R T
files = list.files(file.path(data.dir,"All-Pairs.Natural-Ratio/"), full.names=TRUE)
pchic_files = files[grepl(pattern = "CHiC", x = files)]

cCRE_keys = read.table(file.path(data.dir,"Annotations/hg19-cCREs.bed.gz"),
                       col.names = c("chromosome","start","end","rDHS","cCRE_accession","cCRE_group"))

ENS_keys = read.table(file.path(data.dir,"Annotations/GENCODEv19-TSSs.bed.gz"),
                      col.names = c("chromosome","start","end","transcript_id",".","strand","gene_id"))

#prendre 1 dans la colonne 3
pcHIC = lapply(pchic_files, function(x) read.table(x, 
                col.names = c("cCRE","locus","interaction","cluster")))

names(pcHIC) = c("CD34","GM12878")
pcHIC = lapply(pcHIC, function(x) filter(x, interaction == 1))



save(cCRE_keys,ENS_keys,pcHIC, file = "./data/interaction/rdata/BENGI_benchmark_pchic.RData")

################################################################################
## B E N G I   R A W
data.dir = file.path(main.dir,"data/interaction/pcHiC/BENGI_raw/")
files = list.files(paste0(main.dir,"/","data/interaction/pcHiC/BENGI_raw/"), full.names=TRUE)

raw_bengi_benchmark = list("GM12878" = read.delim(files[1]), 
                           "CD34" = read.delim(files[2]))

raw_bengi_benchmark = lapply(raw_bengi_benchmark, function(data) 
  data[grepl(data$Symbol, pattern = "POLR3G-") & data$log.observed.expected. > 10, ])

curated_benchmark = lapply(raw_bengi_benchmark, function(data)
  select(data, chr, start, end ))

raw_bengi_polr3g = lapply(curated_benchmark, function(data) 
  data %>% mutate(polr3g_interaction = paste0(start,'-',end)))


save(raw_bengi_polr3g, file = "./data/interaction/rdata/raw_BENGI_interactions_POLR3G.RData")

