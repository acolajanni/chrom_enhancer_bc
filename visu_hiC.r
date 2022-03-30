################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : visu_hiC                                                          #                                                        
# > Fonction : visualize features of hiC datasets                              #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")


# I M P O R T
source("./scripts/fun_interactions.r")
source("./scripts/functions_vizu.r")
load(file.path(data.dir,"3div_0.001.Rdata"))



#' Change the format of the original 3div dataframes to something readable for 
#' Hi-C analysis
#'
#' @param df dataframe
#'
#' @return
#' @examples
#' 
hic_format <- function(df){
  hic = df %>% 
    mutate(seq1 = 
      #paste0("chr", formatC(str_remove(substring(frag1, 4 , 5), ":"), width=2, flag="0"))
      #paste0("0",str_remove(substring(frag1, 4 , 5), ":"))
      str_remove(substring(frag1, 1 , 5), ":")
      ) %>% 
    mutate(seq2 = 
      #paste0("chr", formatC(str_remove(substring(frag2, 4 , 5), ":"), width=2, flag="0"))
      str_remove(substring(frag1, 1 , 5), ":")
      ) %>% 
    mutate(
      frag1 = str_remove(frag1, pattern = "chr(..?):"),
      frag2 = str_remove(frag2, pattern = "chr(..?):")
    ) 
  return(hic)
}

#' Count the number of contact for each chromosome (in column "seq")
#'
#' @param df dataframe
#'
#' @return list of dataframe with 2 column (chromosome, count value)
count_chrom = function(df){
  return(data.frame(table(df$seq1)))}

################################################################################

result = extract_features_from_list(File_list_filtered0.001,"dist")
data = File_list_filtered0.001
rm(File_list_filtered0.001)

# distance dataframe
distances = list_to_2colDF(result)

# Number of contact dataframe 
interact = lapply(result,length)
interact = list_to_2colDF(interact)

# contact number per chromosomes per cell line
data = lapply(data, hic_format )
inter_by_chrom = lapply(data, count_chrom)
df <- ldply (inter_by_chrom, data.frame)
colnames(df) = c("cell_line","chr","contact")




ad2 = data$AD2
ft2 = data$FT2

A = data.frame(table(ft2$seq1))


sum(count_chrom(ad2)$x)
sum(count_chrom(ft2)$x)


sum(count_chrom2(ad2))
sum(count_chrom2(ft2))





# Vizu : 
hist_contact(distances, 
  title = "Distribution of contact distances of the cell lines from 3div database (PCHi-C filtered at p < 0.05)")
contact_distrib(interact, 
  title = "Number of contact of the cell lines from 3div database (PCHi-C filtered at p < 0.05)")


cells = unique(df$cell_line)

df1 = df[df$cell_line %in% cells[1:8], ]
df2 = df[df$cell_line %in% cells[9:17],]
df3 = df[df$cell_line %in% cells[18:28],]

ggplot(df3, aes(x=chr, y=contact)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ cell_line, nrow = 3 ) +
  theme(axis.text.x = element_text(angle = 90))




################################################################################
# Chercher les intéractions :
# POLR3G chromosomic location : https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000113356;r=5:89767565-89810370
#Chromosome 5: 89767565-89810370
load("./data/phic/rdata/3div_0.05.Rdata")
data = File_list_filtered0.05
rm(File_list_filtered0.05)
data = lapply(data, hic_format )

polr3g = GRanges(
  seqnames = Rle("chr 5"),
  ranges = IRanges(start = 89,767,565, end = 89,810,370)
)
polr3g = GRanges(
  seqnames = Rle("chr 5"),
  ranges = IRanges(start = 89770680, end = 89777400)
) 

A = as.vector(seqnames(polr3g))


range_frag1 = GRanges(
  seqnames = Rle(df$seq1),
  IRanges(df$frag1)
)
range_frag2 = GRanges(
  seqnames = Rle(df$seq2),
  IRanges(df$frag2)
)



#' Subset a dataframe by its chromosome number in the column "seq1"
#'
#' @param df dataframe
#' @param chr string
#'
#' @return subset of the original dataframe
#' @export
#'
#' @examples
subset_chr<-function(df,chr){
  return(filter(ad2, seq1 == chr))}

#' Get the regions that overlap with the Region Of Interest (ROI)
#'
#' @param GRange GenomicRange object, region to compare to the ROI
#' @param ROI GenomicRange object, Region of interest
#'
#' @return Character string vector. 
#' String format : "position1-position2"
#' @export
#'
#' @examples
get_overlappped_regions<-function(GRange,ROI) {
  query = subsetByOverlaps(GRange,ROI)
  if (length(query) == 0){return(NA)}
  range = ranges(query)
  start = range@start
  width = range@width
  return(paste0(start,"-",start+width-1))
  
}

#' Get the regions that overlap with the Region Of Interest (ROI) for both
#' fragments that are in contact in HiC experiment
#'
#' @param df dataframe
#' @param ROI GenomicRange object, Region of interest
#'
#' @return dataframe with two columns : 
#' "frag1" contains the regions that overlap the ROI for the 1st fragment
#' "frag2" contains the regions that overlap the ROI for the 2nd fragment
#' @export
#'
#' @examples
get_all_overlapped_regions<-function(df,ROI){
  df = subset_chr(df , chr = as.vector(seqnames(ROI)))
  
  f1 = chr_location_to_ranges(df, "frag1")
  ROI_reg = unique(get_overlappped_regions(f1,ROI))
  
  f2 = chr_location_to_ranges(df, "frag2")
  ROI_reg2 = unique(get_overlappped_regions(f2,ROI))
  return(data.frame("frag1" = ROI_reg, "frag2" = ROI_reg2))
}

#' Subset the dataframe to only keep the regions that are supposed to be in 
#' contact with the ROI
#'
#' @param df dataframe
#' @param ROI GenomicRange object, Region of interest
#'
#' @return dataframe
#' @export
#'
#' @examples
subset_by_overlap<-function(df,ROI){
  interactions = get_all_overlapped_regions(df,ROI)
  
  return( df %>% 
    mutate(
      is_overlapped_f1 = ifelse(frag1 %in% interactions$frag1, TRUE, FALSE),
      is_overlapped_f2 = ifelse(frag2 %in% interactions$frag2, TRUE, FALSE)
    ) %>% filter(
      xor(is_overlapped_f1,is_overlapped_f2)
    ) %>% mutate(
      interaction_polr3g = 
        ifelse(!is_overlapped_f1 & is_overlapped_f2, yes = frag1, no = frag2)
    ) %>% select(
      c(interaction_polr3g,freq,dist,X.log10.result.) )
  ) }






inters = lapply(data, function(x) subset_by_overlap(ROI = polr3g, df = x))












################################################################################
load("./data/phic/rdata/3div_0.05.Rdata")
data = File_list_filtered0.05
rm(File_list_filtered0.05)
data = lapply(data, hic_format )

polr3g = GRanges(
  seqnames = Rle("chr 5"),
  ranges = IRanges(start = 89770680, end = 89777400)
  ) 

inters = lapply(data, function(x) get_overlap_region(ROI = polr3g, df = x))


ft = inters$HCmerge
bornes = get_individual_position(ft)

frag1 = GRanges(
  seqnames = Rle(unique(ft$seq1)),
  ranges = IRanges(unique(ft$frag1))
)
frag2 = GRanges(
  seqnames = Rle(unique(ft$seq2)),
  ranges = IRanges(unique(ft$frag2))
)

#query_f2 = subsetByOverlaps(frag2,polr3g)

query_f1 = subsetByOverlaps(frag1,polr3g)

A = ranges(query_f1)
start = A@start
width = A@width

A = paste0(start,"-",start+width-1)



get_overlappped_regions<-function(GRange,ROI) {
  query = subsetByOverlaps(GRange,ROI)
  if (length(query) == 0){return(NA)}
  range = ranges(query)
  start = range@start
  width = range@width
  return(paste0(start,"-",start+width-1))

}


f1 = chr_location_to_ranges(ft, "frag1", unique = TRUE)
ROI_reg = get_overlappped_regions(f1,polr3g)

f2 = chr_location_to_ranges(ft, "frag2", unique = TRUE)
ROI_reg2 = get_overlappped_regions(f2,polr3g)






################################################################################
query_f2
query_f1


length(unique(bornes$start1))
unique(bornes$end1)

length(unique(bornes$start2))
length(unique(bornes$end2))


#in_POLR3G_frag1 = unique(as.matrix(ranges(query_f1))[,1])
#in_POLR3G_frag2 =c(in_POLR3G, unique(as.matrix(ranges(query_f2))[,1]))


#test = bornes[ !bornes$start1 %in% in_POLR3G_frag1 ,]
#test2 = bornes[ !bornes$start2 %in% in_POLR3G_frag2 ,]





















