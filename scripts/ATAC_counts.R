################################################################################                                                                            #
# > April 2022                                                                                                                 
# > Script : ATAC_counts.R                                                                                             
# > Function : Analyse ATAC.counts file
# @ COLAJANNI Antonin                                                          
################################################################################

## ATAC_seq
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = "./data/interaction/ATAC/"

# I M P O R T
library(edgeR)
source("./scripts/fun_interactions.r")

# Read table
count.dir = list.files("./", pattern = "atac.counts.txt$", recursive = TRUE)
atac_count = as.data.frame(read_table(count.dir, skip = 1))

# change format of the "Chr" column
atac_count = mutate(atac_count, Chr = paste0("chr",Chr) )

# Select columns
experiment = colnames(atac_count)[7:length(colnames(atac_count))]
count = select(atac_count, all_of(c("Geneid",experiment)))

# change format
rownames(count) = count$Geneid
count$Geneid = NULL

# EdgeR treatment
design = c(
  rep("Par",2),
  rep("BrM",2),
  rep("LM",2))

DGE =  DGEList(counts = count, group = factor(design))

# Filtering lowly present region
keep = filterByExpr(DGE, group = design)
DGE_filtered = DGE[keep,]

# Extract count data
DGE_filtered = calcNormFactors(DGE, method="TMM")
count_filtered <- as.data.frame(cpm(DGE_filtered, normalized.lib.sizes=TRUE, log = F))

# Retrive rows that pass filters
atac_count = atac_count[atac_count$Geneid %in% rownames(count_filtered) , ]

# change row values with normalized ones
count_filtered$Geneid = rownames(count_filtered)
count_filtered = merge(select(atac_count,Geneid,Chr,Start,End), count_filtered, by = "Geneid" )
count_filtered$Geneid = NULL

save_path = "./data/interaction/bedgraph_analysis/bed_file/"
file_name = "GSE129646_atac_counts_"
# Change to 3x2(replicates) bedgraph format and save them
for (exp in experiment) {
  tmp_filename = paste0(file_name,exp,".bedgraph")  
  tmp = select(count_filtered, all_of(c("Chr","Start","End",exp)))  
  write.table(tmp, file = file.path(save_path,tmp_filename), quote = F, row.names = F, col.names = F)
}





