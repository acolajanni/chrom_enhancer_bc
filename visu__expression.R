# > February 2022                                                              #                                                
# > Script : visu expresion                                                    #                                                        
# > Fonction : visualisation data depmap                                       # 
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/expression/rdata")


# I M P O R T
library(stringr)
library(dplyr)
library(reshape2)
source("./scripts/Normalization.R")
load(file.path(data.dir,"cancer_healthy_expression.RData"))

# change format
cell_lines = names(expr)[! names(expr) %in% c("gene_name","ENSEMBL_id")]
expr = select(expr, -gene_name)
# as numeric
expr <- cbind("ENSEMBL_id"=expr$ENSEMBL_id, 
              mutate_all(select(expr,all_of(cell_lines)), function(x) as.numeric(as.character(x))))

# Remove lowly expressed genes / 1 duplicated id
expr <- expr %>%
  filter(! duplicated(ENSEMBL_id))

row.names(expr) = expr$ENSEMBL_id
expr$ENSEMBL_id = NULL

Row_sums = rowSums(expr)
filtered_expr = expr[Row_sums > 3*ncol(expr),]
  






##################################
# Normalisation 
source("./scripts/Normalization.R")
library(edgeR)
library(DESeq2)

# Parameters
design = c(rep("CCLE", 20 ), rep("RoadMap",length(cell_lines)-20))
count.matrix = data.matrix(filtered_expr)
tools=c("TMM", "TMMwsp", "RLE", "upperquartile")
##     ##     ##     ##     ##     ##     ##     ##     ##     ## 
    # Testing normalizations : edgeR
norm_list = list()

tool = 'TMM'
for (tool in tools){
  edgeR.dgelist = DGEList(counts = count.matrix, group = factor(design))
  nf = calcNormFactors(edgeR.dgelist, method=tool)
  edgeR_norm <- cpm(nf, normalized.lib.sizes=TRUE, log = TRUE)
  #melted_expr_edgeR = reshape2::melt(edgeR_norm)
  norm_list[[tool]] = edgeR_norm  
}
##     ##     ##     ##     ##     ##     ##     ##     ##     ## 
    ## Testing normalizations : DESeq2
deseq_norm = tools.norm.RNAseq(count.matrix, tool = "vst2", design = design )
#melted_expr = reshape2::melt(deseq_norm)
norm_list[["DESeq2"]] = deseq_norm

## Applying limma batch correction
norm_list_limma = lapply(norm_list, function(x) removeBatchEffect(x, design) )
norm_list_limma = lapply(norm_list, function(x) reshape2::melt(x) )


##     ##     ##     ##     ##     ##     ##     ##     ##     ## 
    ## Batch correction : ComBat_seq
adjusted = ComBat_seq(count.matrix, batch=design, group=NULL)

norm_list_combat = list()

for (tool in tools){
  edgeR.dgelist = DGEList(counts = adjusted, group = factor(design))
  nf = calcNormFactors(edgeR.dgelist, method=tool)
  edgeR_norm <- cpm(nf, normalized.lib.sizes=TRUE, log = TRUE)
  #melted_expr = reshape2::melt(edgeR_norm)
  norm_list_combat[[tool]] = edgeR_norm  
}
# Testing normalizations : DESeq2
deseq_norm = tools.norm.RNAseq(adjusted, tool = "vst2", design = design )
#melted_expr = reshape2::melt(deseq_norm)
norm_list[["DESeq2"]] = deseq_norm

norm_list_combat = lapply(norm_list, function(x) reshape2::melt(x) )


##     ##     ##     ##     ##     ##     ##     ##     ##     ## 
  ## Batch correction : RUVseq
library(RUVSeq)




##     ##     ##     ##     ##     ##     ##     ##     ##     ## 



################################################################################
# Une seule normalisation, plusieurs corrections de batch
# Parameters
design = c(rep("CCLE", 20 ), rep("RoadMap",length(cell_lines)-20))
count.matrix = data.matrix(filtered_expr)

 # LIMMA
edgeR.dgelist = DGEList(counts = count.matrix, group = factor(design))
nf = calcNormFactors(edgeR.dgelist, method=tool)
edgeR_norm <- cpm(nf, normalized.lib.sizes=TRUE, log = TRUE)
edgeR_norm = removeBatchEffect(edgeR_norm, design)
edgeR_norm = reshape2::melt(edgeR_norm)

 # ComBat-seq
library(sva)
adjusted = ComBat_seq(count.matrix, batch=design, group=NULL)
edgeR.dgelist = DGEList(counts = adjusted, group = factor(design))
nf = calcNormFactors(edgeR.dgelist, method=tool)
ComBat_norm <- cpm(nf, normalized.lib.sizes=TRUE, log = TRUE)
ComBat_norm = reshape2::melt(ComBat_norm)

  # RUVseq
edgeR.dgelist = DGEList(counts = count.matrix, group = factor(design))
nf = calcNormFactors(edgeR.dgelist, method=tool)
nf = estimateGLMCommonDisp(nf)
nf = estimateGLMTagwiseDisp(nf)
fit <- glmFit(nf)
res <- residuals(fit, type="deviance")
colnames(res) = cell_lines
## Fonctionne pas
################################################################################


Roadmap_ensembl = read.table(file.path(data.dir,"./RoadMap/Ensembl_v65.Gencode_v10.ENSG.gene_info"), 
                             header = FALSE, col.names = c("ID","chr","start","end","transcription","type","gene_name","hgnc") )



ctcfl = "ENSG00000124092"
ctcf =  "ENSG00000102974"
pou5f1 ="ENSG00000204531"
myc =   "ENSG00000136997"
polr3gl="ENSG00000121851"
polr3g ="ENSG00000113356"
#melted_expr[grepl(pattern = polr3g, melted_expr$Var1),]

CCLE = cell_lines[1:20]
ROADMAP = cell_lines[21:length(names(expr))]
#melted_expr = reshape2::melt(edgeR_norm)

condition = ifelse(ComBat_norm$Var2 %in% CCLE, yes = "CCLE", no = "RoadMap")

# cancer lines of interest : 
of_interest = c("MDAMB231","MDAMB436","HCC1937", "HCC1806",
                "HCC1954","MDAMB468","K562")
color_cell = ifelse(grepl(paste(of_interest,collapse = "|"),  cell_lines ), yes = "red", no ="black")


expression_plot = function(melted_df){
plot = ggplot(data = melted_df, aes(x=Var2, y=value, color = condition)) + 
  geom_boxplot()   +
  geom_text(aes(x=Var2,
                label=ifelse(grepl(polr3g,Var1),".",''))
            , size= 10, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = color_cell),
          axis.text=element_text(size=6))  +
  xlab("cell_lines") + ylab("normalized_expression")
  return(plot)
}
  
# Raw + log2
melt_raw = reshape2::melt(log2(count.matrix+1))
expression_plot(melt_raw) + ggtitle("Expression boxplot : Raw expression (log2+1)")

# Norm
deseq_norm = tools.norm.RNAseq(count.matrix, tool = "vst2", design = design )
melt_norm = reshape2::melt(deseq_norm)
expression_plot(melt_norm) + ggtitle("Expression boxplot : Normalized expression")

# Norm + batch(limma)
batch_norm = removeBatchEffect(deseq_norm, design)
batch_norm_melt = reshape2::melt(batch_norm)
expression_plot(batch_norm_melt) + ggtitle("Expression boxplot : Normalized expression + batch correction (limma)")

# Norm + batch(ComBat)
deseq_norm_combat = tools.norm.RNAseq(adjusted, tool = "vst2", design = design )
batch_norm_combat_melt = reshape2::melt(deseq_norm_combat)
expression_plot(batch_norm_combat_melt) + ggtitle("Expression boxplot : Normalized expression + batch correction (Combat)")

# Raw + batch(ComBat)
batch_raw_combat_melt = reshape2::melt(log2(adjusted+1))
expression_plot(batch_raw_combat_melt) + ggtitle("Expression boxplot : Raw expression (log2+1) + batch correction (Combat)")


# raw + batch(limma)
batch_raw = removeBatchEffect(log2(count.matrix+1), design)
batch_raw_melt = reshape2::melt(batch_raw)
expression_plot(batch_raw_melt) + ggtitle("Expression boxplot : Raw expression (log2+1) + batch correction (limma)")





combat_plot = expression_plot(ComBat_norm) + ggtitle("Expression boxplot : ComBat_seq")
combat_plot

edgeR_plot = expression_plot(edgeR_norm) + ggtitle("Expression boxplot : limma batch correction")
edgeR_plot

###############################################################################"
# Heatmap

of_interest = c(
  "CTCFL" = "ENSG00000124092",
  "CTCF" =  "ENSG00000102974",
  "POU5F1" ="ENSG00000204531",
  "MYC" =   "ENSG00000136997",
  "POLR3GL"="ENSG00000121851",
  "POLR3G" ="ENSG00000113356" )
of_interest = data_frame("ENS_id" = of_interest, "gene_name" = names(of_interest))




tmp = filter(ComBat_norm, grepl(paste(of_interest$ENS_id,collapse="|"),edgeR_norm$Var1))
tmp$copie = substr(tmp$Var1,1,15)



to_heatmap = select(merge(tmp,of_interest, by.x = "copie", by.y = "ENS_id"),
                    gene_name,Var2 ,value)
colnames(to_heatmap) = c("gene_name","cell_line","expression")


heat = reshape2::dcast(to_heatmap, gene_name~cell_line)
row.names(heat) = heat$gene_name
heat$gene_name = NULL
annotation_heatmap = data.frame("data_source" = design, row.names = cell_lines)

color = c("red","blue")
names(color) = c("CCLE","RoadMap")



TNBC = c("MDAMB231","MDAMB436","HCC1937", "HCC1806",
         "HCC1954","MDAMB468")

annotation_heatmap$state = ifelse(
  grepl(paste(TNBC,collapse="|") ,row.names(annotation_heatmap)), 
  yes = "TNBC", no = "other")

colors = list(data_source = color)

map = pheatmap(t(heat),
         #annotation_col = annotation_heatmap,
         #annotation_row = colors,
         cutree_rows = 4,
         angle_col = 315,
         fontsize_row = 5)
         

#map


tmp = data.frame(t(heat))
#tmp["colors",] <-ifelse(grepl(paste(TNBC,collapse="|") ,row.names(annotation_heatmap)), 
#                       yes = "red", no = "white")

tmp$colors = ifelse(grepl(paste(TNBC,collapse="|") ,row.names(annotation_heatmap)), 
                                                  yes = "red", no = "darkgrey")

col = tmp[order(match(colnames(heat), map$gtable$grobs[[5]]$label)),]$colors
#col = select(tmp, all_of(col))["colors",]

map$gtable$grobs[[5]]$gp=gpar(col= col, fontsize = 5)
map



######################################################
library(pheatmap)

myColor <- colorRampPalette(c("blue", "white", "red"))(70)


#color = brewer.pal(11, name = "RdBu")
#color_scale = rev(c(color[1:5],"#FFFFFF","#FFFFFF",color[7:11]))
#show_col(myColor)

test = heat_norm


get_color_scale_pheatmap = function(matrix,ncolors = 50,threshold = 1){
  
  range <- max(abs(matrix))
  myBreaks = seq(-range, range, length.out = ncolors)
  myBreaks[myBreaks > -threshold  &  myBreaks < threshold ] = 0
  myBreaks = unique(myBreaks)
  
  paletteLength <- length(myBreaks)
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-1)
  return(list("colorscale" = myColor, "breaks" = myBreaks)) }

range <- max(abs(heat_norm))
myBreaks = seq(-range, range, length.out = 50)
myBreaks[myBreaks > -1  & myBreaks < 1 ] = 0
myBreaks = unique(myBreaks)

paletteLength <- length(myBreaks)
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-1)


scale = get_color_scale_pheatmap(heat_norm,25)

test = data.frame(t(heat_norm))
test = test[order(-test$POLR3G),]


pheatmap(test,
         color = scale$colorscale,
         breaks = scale$breaks,
         cluster_rows = F, cluster_cols = F,
         cutree_rows = 2,
         angle_col = 0,
         fontsize_row = 5)
