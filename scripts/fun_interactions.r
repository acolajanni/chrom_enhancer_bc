################################################################################                                                                            
# > February 2022                                                              #                                                
# > Script : visu_hiC                                                          #                                                        
# > Fonction : visualize features of hiC datasets                              #        
# @ COLAJANNI Antonin                                                          #
################################################################################




library(InteractionSet)
library(tidyverse)
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(Gviz)




#' Title
#'
#' @param df 
#' @param feature 
#'
#' @return
#' @export
#'
#' @examples
chr_location_to_ranges <- function(df,feature){
  if (feature == "frag1"){
    seqnames = df$seq1
    reg = df$frag1
  }
  else if(feature == "frag2"){
    seqnames = df$seq2
    reg = df$frag2
  }
  else{ stop("enter a valid feature name") }
  
  ranges = IRanges(reg)
  return(GRanges(Rle(seqnames),ranges)) }


get_interaction_set <- function(df){
  source = chr_location_to_ranges(df, "frag1")
  target = chr_location_to_ranges(df, "frag2")
  gi = GInteractions(source, target)
  gi$freq = df$freq
  return(list("interaction_set" = gi, "pos1" = source, "pos2" = target) )
}


get_individual_position <- function(df){
  
  range_start = as.data.frame(str_split_fixed(df$frag1, pattern = "-", n=2))
  colnames(range_start) = c("start1","end1")
  
  range_start2 = as.data.frame(str_split_fixed(df$frag2, pattern = "-", n=2))
  colnames(range_start2) = c("start2","end2")
  
  df = cbind(range_start,range_start2)
  df = apply(df, 2 , function(x) as.numeric(as.character(x)))
  return(data.frame(df))
}




#' Old version of "subset_by_overlap" : fonctionne mal
#'
#' @param ROI region of interest
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
get_overlap_region <-function(ROI, df){
  
  #if (missing(ROI)){
  #  ROI = GRanges(
  #    seqnames = Rle("chr5"),
  #    ranges = IRanges(start = 89767565, end = 89810370)
  #  )
  #}
  interactions = get_interaction_set(df)
  
  Olap = findOverlaps(interactions$interaction_set,ROI,type ="any",use.region="both")
  R1 = unique(as.matrix(ranges(interactions$pos1[queryHits(Olap)]))[,1])
  R2 = unique(as.matrix(ranges(interactions$pos2[queryHits(Olap)]))[,1])
  
  ranges = get_individual_position(df)
  interactions = df[ranges$start1 %in% R1 | ranges$start2 %in% R2, ]
  return( interactions )
}




################################################################################
# Fonctions V2


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
  return(filter(df, seq1 == chr))}

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
  
  #if (length(ROI_reg) == length(ROI_reg2)){
  #  return(data.frame("frag1" = ROI_reg, "frag2" = ROI_reg2)) 
  #}
  return(list("frag1" = ROI_reg, "frag2" = ROI_reg2))  
}

#' Subset the dataframe to only keep the regions that are supposed to be in 
#' contact with the ROI
#'
#' @param df dataframe
#' @param ROI GenomicRange object, Region of interest
#' @param remove_col logical TRUE to remove unnecessary columns
#'
#' @return dataframe
#' @export
#'
#' @examples
subset_by_overlap<-function(df,ROI,remove_col = TRUE){
  interactions = get_all_overlapped_regions(df,ROI)
  
  df = df %>% dplyr::mutate(
      is_overlapped_f1 = ifelse(frag1 %in% interactions$frag1, TRUE, FALSE),
      is_overlapped_f2 = ifelse(frag2 %in% interactions$frag2, TRUE, FALSE)
    ) %>% dplyr::filter(
      xor(is_overlapped_f1,is_overlapped_f2)
    ) %>% dplyr::mutate(
      interaction_polr3g = 
        ifelse(!is_overlapped_f1 & is_overlapped_f2, yes = frag1, no = frag2)) 
  
  if (remove_col){
    return (dplyr::select(df,c(frag1,frag2,interaction_polr3g
                               ,frag1_cov,frga2_cov,freq,dist,X.log10.result.)))
  }
  return(df)
}



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
      str_remove(substring(frag1, 1 , 5), ":")
    ) %>% 
    mutate(seq2 = 
      #paste0("chr", formatC(str_remove(substring(frag2, 4 , 5), ":"), width=2, flag="0"))
      str_remove(substring(frag2, 1 , 5), ":")
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


count_row = function(df){
  return(nrow(df))}

list_to_2colDF <- function(data){
  df = data.frame(do.call(cbind, data))
  return (melt(df))
}




#' Dataframe with chromosomic location "XXXX-XXXX" format
#'
#' @param df dataframe
#' @param start int
#' @param end int 
#'
#' @return dataframe
#' @export
#'
#' @examples
subset_by_interval<-function(df, start, end ){
  positions = get_individual_position(df)
  return(
    df[positions$start1 >= start & positions$start1 <= end &
         positions$start2 >= start & positions$start2 <= end, ]
  )}

character_to_nums = function(char){
  char = str_split(char, "-", simplify = TRUE)
  return(as.numeric(char))
}

#' Returns all the chromosomal positions and the bins they are associated with 
#'
#' @param df dataframe with chromosomal positions
#' @param bornes dataframe with start & end columns
#' @param binsize 
#' @param feature 
#'
#' @return dataframe
#' @export
#'
#' @examples
get_position_inside_bins = function(df,bornes,binsize,feature){
  if (feature == "frag1"){
    start = "frag1_start"
    end = "frag1_end" }
  else if(feature == "frag2"){
    start = "frag2_start"
    end = "frag2_end" }
  else{ stop("feature name not recognized")}
  
  binsize = binsize-1
  
  include_in_bornes = lapply(unique(df[[feature]]), function(x){
    pos = character_to_nums(x)
    sup = ifelse(bornes$start > pos[1]-binsize, TRUE, FALSE )
    inf = ifelse(bornes$end < pos[2]+binsize, TRUE, FALSE )
    return(bornes[sup & inf,]) })

  names(include_in_bornes) = unique(df[[feature]])
  return(ldply(include_in_bornes))  }

filter_pvalue <-function(df,seuil = 2){
  return(filter(df,X.log10.result. >= seuil))}



#### BED creation 

#' Transform a vector of numerical value into RGB value
#'
#' @param score vector of numerical value
#' @param color_scale character string (viridis / plasma)
#'
#' @return dataframe with 3 columns : red,green,blue
#' @export
#'
#' @examples
num_to_colors = function(score, color_scale = "viridis"){
  
  if (color_scale == "viridis"){
    return(colors = data.frame(t(col2rgb(viridis(max(score),direction = -1)
                                         [cut(score, max(score))] ))) ) } 
  
  else if (color_scale == "plasma"){
    return(colors = data.frame(t(col2rgb(plasma(max(score))
                                         [cut(score, max(score))] ))) )}
  
  else if (color_scale == "cividis"){
    return(colors = data.frame(t(col2rgb(cividis(max(score),direction = -1)
                                       [cut(score, max(score))] ))) ) }
  
  else if (color_scale == "mako"){
    return(colors = data.frame(t(col2rgb(mako(max(score),direction = -1)
                                         [cut(score, max(score))] ))) ) }
    
    
}

#' Transform dataframe with RGB values (3 columns) into one column
#'
#' @param df_rgb dataframe with 3 columns : red,green,blue
#'
#' @return dataframe with one column "color"
#' @export
#'
#' @examples
rgb_to_1col = function(df_rgb){
  return(df_rgb%>% mutate(
    color = paste0(red,",",green,",",blue)
  ) %>% select(color)) }


#' modify a dataframe with bed format to make a "full bed" format
#'
#' @param df Datframe with 5 columns "chr","start","end","score","name"
#'
#' @return Dataframe with 4 other columns strand, thick1, thick2, color
#' @export
#'
#' @examples
full_bed_format = function(df,colorscale = "viridis"){
  
  colors = num_to_colors(df$score, color_scale = colorscale)
  colors = rgb_to_1col(colors)
  
  df_col = cbind(df,colors) %>%  mutate(
    strand = ".",
    thick1 = start,
    thick2 = end
  ) %>% select(chr,start,end,name,score,strand,thick1,thick2,color)
  
  return(df_col) }



#' Change the color-bed format to bigInteract format 
#'
#' @param bed dataframe with bed 9columns 
#' chr,start,end,name,score,strand,thick1,thick2,color
#'
#' @return dataframe with 18 columns (bigInteract format)
#' @export
#'
#' @examples
bed_to_bigInteract = function(bed){
  bed$value = bed$score
  bed$exp = "."
  bed = select(bed, chr,start,end,name,score,value,exp,color)
  
  bed$sourceChrom = bed$chr
  bed$sourceStart = bed$start
  bed$sourceEnd = bed$start
  bed$sourceName = NA
  bed$sourceStrand = "."
  
  bed$targetChrom = bed$chr
  bed$targetStart = bed$end
  bed$targetEnd = bed$end
  bed$targetName = NA
  bed$targetStrand = "."
  return(bed)
}

#' Convert bigwig file to a bed-type dataframe
#'
#' @param file_path path to a bigwig file
#'
#' @return dataframe with 4 columns (chr, strat, end, score)
#' @export
#'
#' @examples
BigWig_to_bed = function(file_path, liftover=FALSE, lift_path = NULL){
  gr = rtracklayer::import.bw(file_path)
  if (liftover){
    chain <- rtracklayer::import.chain(lift_path)
    gr <- unlist(rtracklayer::liftOver(gr, chain))
  }
  
  bed <- data.frame(chr=seqnames(gr),
                    start=start(gr)-1,
                    end=end(gr),
                    score= gr$score)
  
  if (!str_detect(bed$chr[5], "chr")){
    bed$chr = paste0("chr",bed$chr)
  }
  
  return(bed) }



#' Convert bed dataframe to bedGraph format
#'
#' @param bed dataframe with 4 columns (chr, strat, end, score) 
#' @param name name of the bed file
#' @param color color to display the signal
#'
#' @return
#' @export
#'
#' @examples
bed_to_bedGraph = function(bed, name = "test",color = "0,0,255"){
  track_type = "bedGraph"
  options = 'visibility=full yLineOnOff=on autoScale=on yLineMark="0.0" alwaysZero=on graphType=bar maxHeightPixels=128:75:11 windowingFunction=maximum smoothingWindow=off'
  first_line = paste0("track type=",track_type, " name=", name, " description=", name, " color=",color, " ",options )
  headline = data.frame(chr = first_line,
                        start = "",
                        end = "",
                        score = "")
  
  
  bedGraph = rbind(headline, bed)
  return(bedGraph)
}

#' Convert a bigwig file (from path) into a bedgraph file
#' 
#' bedGraph are used to plot an histogram of signal in UCSC genome browser
#'
#' @param bigwig_path path to a bigwig file
#' @param file_path path to save the .bedGraph file (FALSE if you don't want to save a bedgraph file) 
#' @param chr_filter logical value to apply a filter
#' @param chr_to_keep character string "chr8" for example
#' @param name character string. name of the bedgraph file
#' @param color character string RGB format : 0,0,255 for example 
#'
#' @return
#' @export
#'
#' @examples
BigWig_to_bedGraph = function(bigwig_path, file_path, chr_filter = TRUE, chr_to_keep = "chr5", 
                              name = "test",color = "0,0,255", apply.log = TRUE, region_to_keep = c(87500000,91000000) ){
  bed = BigWig_to_bed(bigwig_path)
  
  if (missing(region_to_keep)){
    region_to_keep = c(min(bed$start),max(bed$end))
  }
  
  if (chr_filter){
    bed = filter(bed, chr == chr_to_keep)
    bed = filter(bed, between(start,region_to_keep[1],region_to_keep[2]))
  }
  
  #if (max(abs(bed$score)) > 20 ){ apply.log = TRUE }
  #else {apply.log = FALSE}
  
  if (apply.log) {
    bed$score = log2(bed$score+1)
  }
  
  
  bedGraph= bed_to_bedGraph(bed, name, color)
  if (!is.logical(file_path)) {
    write.table( bedGraph, file_path, sep="\t", quote = F, row.names = F, col.names = F ) }
  
  return(bedGraph)
}

#' Merge and log2 normalize plus and minus strand score
#'
#' @param plus dataframe : bed format (chr,start,end,score)
#' @param minus dataframe : bed format (chr,start,end,score) 
#' @param name character string, name of the bedgraph file
#' @param color character string, RGB format (0,0,255)
#' @param chr_filter logical value to apply a filter
#' @param chr_to_keep character string "chr8" for example
#'
#' @return
#' @export
#'
#' @examples
merge_plus_minus_GROseq = function(plus_strand,minus_strand,name,color, chr_filter = TRUE, chr_to_keep = "chr5"){
  
  if (chr_filter){
    plus_strand = filter(plus_strand, chr == chr_to_keep)
    minus_strand = filter(minus_strand, chr == chr_to_keep) 
  }
  
  minus_strand$score = log2(abs(minus_strand$score)+1)
  plus_strand$score = log2(abs(plus_strand$score)+1)
  minus_strand$score = minus_strand$score * -1
  merged = rbind(minus_strand,plus_strand)
  merged_bedGraph = bed_to_bedGraph(merged, name = name, color = color)
  return(merged_bedGraph)
}


#' Convert multiple bigwig files (in path format) in bedgraph 
#'
#' @param paths vector of character string. Paths to bigwig files
#' @param filenames vector of the name of the file to be created
#' @param save_file character string. Path to save the bedfiles
#' @param apply.log logical. 
#'
#' @return
#' @export
#'
#' @examples
multiple_bigWig_to_bedGraph = function(paths,filenames,save_file, apply.log = TRUE){
  count = 1  
  ATAC_dataset = list()
  for(f in paths ){
    print(f)
    name = filenames[count]
    save.path = file.path(save_file, paste0(name,".bedGraph")   )
    
    
    ATAC_dataset[[ filenames[count] ]] = BigWig_to_bedGraph(bigwig_path = f, 
                                                            file_path = save.path,
                                                            name = name,
                                                            apply.log = apply.log)

    count = count+1
  }
  return(ATAC_dataset) }
