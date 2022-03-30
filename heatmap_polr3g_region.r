################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : heatmap around polr3g + GVIZ                                      #                                                        
# @ COLAJANNI Antonin                                                          #
################################################################################


# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")


# I M P O R T
source("./scripts/fun_interactions.r")
load(file.path(data.dir,"3div_0.05.Rdata"))


data = File_list_filtered0.05
rm(File_list_filtered0.05)
data_05 = lapply(data, hic_format )


range_heatmap = 5e+6
binsize = 10000

heatmap_range = GRanges(
  seqnames = Rle("chr5"),
  IRanges(89768410-range_heatmap, 89777404+range_heatmap)
)


# bins of 5kb for heatmap visualisations center around polr3g bait
borne_inf = seq(89768410-range_heatmap , 89777404+range_heatmap, by = binsize)
borne_sup = seq(89768410-range_heatmap+binsize-1 , 89777404+range_heatmap+binsize-1, by = binsize)

# Subset to keep chr5 (faster computation)
chr5 = lapply(data_05, function(x) subset_chr(x, "chr5"))

# Creating the list of subsetted dataframe around the interval 
heatmap = lapply(chr5, function(x) subset_by_interval(x,89768410-range_heatmap, 89777404+range_heatmap) )

# add empty interactions (diagonal) for heatmap visualisation
diagonal = data.frame("frag1"=paste0(borne_inf,"-",borne_sup), 
                      "frag2"=paste0(borne_inf,"-",borne_sup))

heatmap = lapply(heatmap, function(x) rbind(select(x,frag1,frag2), diagonal))


# transform it into a single dataframe
heatmap = ldply(heatmap)

# bins into dataframe for easier steps
bornes = data.frame("start"=borne_inf, "end"=borne_sup)

# correspondence between bins and frag1 and 2 chromosomal positions
include_in_bornes = get_position_inside_bins(heatmap,bornes, binsize , "frag1")
include_in_bornes2 = get_position_inside_bins(heatmap,bornes, binsize, "frag2")


# Merge all the bins interactions
tmp = merge(select(heatmap,frag1,frag2),include_in_bornes, by.x = "frag1", by.y = ".id") 
interactions = merge(tmp, include_in_bornes2, by.x = "frag2", by.y = ".id")


# Middle position of the bin
interactions = interactions %>% mutate(
  frag1 = paste0(start.x,"-",end.x),
  frag2 = paste0(start.y,"-",end.y)
  #frag1 = round((start.x + end.x )/ 2),
  #frag2 = round((start.y + end.y )/ 2)
) %>% select(frag1, frag2)


# Symetric matrix :
interactions_reverse = filter(interactions, frag1 != frag2)
interactions_reverse = interactions_reverse %>% mutate(
  frag1_r = frag2,
  frag2_r = frag1
  ) %>% select(frag1_r, frag2_r)

colnames(interactions_reverse) = c("frag1","frag2")
interactions = rbind(interactions, interactions_reverse)

# To create square matrix
mat = as.matrix(xtabs(~frag1+frag2,interactions, sparse = TRUE))

# row/col annotation
labs.row = ""
labs.row[1:nrow(mat)] = "-"
labs.row[c(nrow(mat)/2,(nrow(mat)/2)-1)] <- "--------- POLR3G"

hm = pheatmap(log2(mat+1),
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         show_rownames = TRUE, 
         show_colnames = TRUE,
         scale = "none",
         #breaks = breaks,
         #inferno(8, direction = -1)
         brewer.pal(8, name = "YlOrRd"),
         fontsize = 12,
         labels_row = labs.row,
         labels_col = labs.row)  





################################################################################
# G V I Z

polr3g_exact = GRanges(
  seqnames = Rle("chr5"),
  IRanges(89768410, 89777404)
)

interaction_polr3g = 
  lapply(data_05, function(x) subset_by_overlap(ROI = polr3g_exact, df = x))

interaction_polr3g = ldply(interaction_polr3g)

ref = GRanges(
  unique(interaction_polr3g$interaction_polr3g),
  seqnames = Rle("chr5"))  

atrack <- AnnotationTrack(ref)
gtrack <- GenomeAxisTrack(IRanges(start = 89768410 - 2000000, end = 89777404 + 2000000), chromosome = "chr5")


library(rtracklayer)
library(GenomicFeatures)
itrack <- IdeogramTrack(genome="hg19", 
                        chromosome="chr5", # specify chromosome in ucsc naming
                        from =89768410 - 2000000,
                        to=89777404 + 2000000)

ht = HighlightTrack(trackList = list(itrack, gtrack,atrack),
                    start = 89768410, end = 89777404, name = "POLR3G_Bait")

plotTracks(list(ht))

