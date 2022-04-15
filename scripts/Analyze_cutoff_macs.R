################################################################################                                                                            #
# > April 2022                                                                                                                 
# > Script : Analyze cutoff macs                                                                                          
# > Function : Analysis of macs2 cuotff-analysis output
# @ COLAJANNI Antonin                                                          
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir ="/shared/projects/chrom_enhancer_bc/data/interaction/bedgraph_analysis/cutoff/"

# I M P O R T
source("./scripts/fun_interactions.r")

# retrieve cutoff_analysis files
cutoff_files = list.files(data.dir) 
# removing merged file (LM / BrM / Par(mb231))
cutoff_files = cutoff_files[! grepl(pattern = "nomt.bs5.RPKM", cutoff_files)] 


cutoff_list = list()
GSE = unique(substr(cutoff_files, 1 , 9))
for (expr in GSE){
  cutoff_list[[ expr ]] = list()
} 

for (f in cutoff_files) {
  GSE = substr(f, 1, 9)
  cutoff_list[[GSE]] [[ f ]] = read.table(file.path(data.dir,f), header = TRUE)
}

# Filtering to keep pscore of interest
cutoff_list_filtered = lapply(cutoff_list, function(liste)
  lapply(liste, function(x) filter(x, ((npeaks < 80000 & npeaks > 30000 ) | pscore == min(x$pscore) )) ))

# merging dataframe for each sublist
cutoff_score_df = lapply(cutoff_list_filtered, function(liste) liste = ldply(liste))


compute_cutoff_score = function(df){
  x=NULL
  # if not enough peaks and pscore sufficently high :  
  if(max(df$npeaks) < 30000 & between(min(df$pscore), 2, 10 )){
    x=1 }
  
  # If score very low : no peaks, sample should be removed
  else if (mean(df$pscore) < 2){
    x = 0 }
  
  # if the cutoff is good enough
  else if (mean(df$pscore) > 2 & between(mean(df$npeaks), 10000,80000 )){
    df = filter(df, pscore > 2)
    x = mean(df$pscore) }
  
  # If the cutoff analysis could not go low enough
  else { x=10 }
  
  return(x)
}

cutoff_score = lapply(cutoff_score_df, function(x) compute_cutoff_score(x))
# Keep the only good cutoff so far
cutoff_score[["GSE107176"]] = 2
cutoff_score = ldply(cutoff_score)
colnames(cutoff_score) = c("GSE","pscore")

# Import known table location of bedgraph - name of narrowpeak file to create
bedGraph = read.table("./data/interaction/bedgraph_analysis/bedGraph_location.txt")
bedGraph$GSE = substr(bedGraph$V2, 1, 9)

# merge dataframe
bedGraph_score = merge(bedGraph, cutoff_score, by = "GSE")

# keep only interesting columns and rows of exploitable dataset
bedGraph_score = bedGraph_score %>%
  select(V1,V2,pscore) %>%
  filter(pscore != 0)

write.table(bedGraph_score, file = "./data/interaction/bedgraph_analysis/bedGraph_location_score.txt",
            quote = F, row.names = F, col.names = F )



