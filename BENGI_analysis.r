################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : BENGI_analysis                                                    #                                                        
# > Function : Analyse of BENGI pcHIC dataset                                  #        
# @ COLAJANNI Antonin                                                          #
################################################################################


# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")


# I M P O R T
source("./scripts/fun_interactions.r")
load(file.path(data.dir,"/BENGI_benchmark_pchic.RData"))


#' Add columns with genomic location of cCRE ID
#'
#' @param df dataframe with cCRE and ENSEMBL ID
#' @param keys dataframe with cCRE ID and their location
#'
#' @return
#' @export
#'
#' @examples
add_cCRE_location <- function(df,keys){
  merged_df = merge(select(df,cCRE,locus), keys, by.x = "cCRE", by.y = "cCRE_accession")
  if (nrow(merged_df) == 0){
    return(data.frame(NA))
  }
  curated_df = merged_df %>% mutate(
    cCRE = paste0(start,"-",end)
  ) %>% select(cCRE, locus, cCRE_group, start,end,chromosome)

  return(curated_df)
}

#' Add columns with genomic location of ENSEMBL ID
#'
#' @param df dataframe with cCRE and ENSEMBL ID
#' @param keys dataframe with ENSEMBL ID and their location
#'
#' @return
#' @export
#'
#' @examples
add_ENSEMBL_location <- function(df, keys){
  merged_df = merge(df, select(keys, chromosome, start,end,gene_id, transcript_id), by.x = "locus", by.y = "gene_id")
  curated_df = merged_df %>% mutate(
    locus_location = paste0(chromosome,':',start,"-",end)
  ) %>% select(cCRE, locus, locus_location , transcript_id,start,end)
  names(curated_df)[names(curated_df) == "start"] = "locus_start"
  names(curated_df)[names(curated_df) == "end"] = "locus_end"
  return(curated_df)
}


#' Add columns with the genomic location of cCRE and ENSEMBL ID
#'
#' @param df dataframe with cCRE and ENSEMBL ID
#' @param cCRE dataframe with cCRE ID and their location
#' @param ENSEMBL dataframe with ENSEMBL ID and their location 
#'
#' @return dataframe
#' @export
#'
#' @examples
add_ID_location <- function(df, cCRE, ENSEMBL){
  cCRE_df = add_cCRE_location(df, cCRE)
  ENS_df = add_ENSEMBL_location(df, ENSEMBL)
  
  cols_to_keep = c("cCRE_location","locus_location","cCRE_group",
                   "cCRE_start","cCRE_end","locus_start","locus_end",
                   "transcript_id")
  
  merged_df = merge(cCRE_df,ENS_df)
  curated_df = merged_df[, colnames(merged_df)%in%cols_to_keep]
  return(arrange(curated_df, cCRE_start, locus_start))
  
}
################################################################################

# Filter only the part that interest us (chr5)
cCRE_keys = filter(cCRE_keys, chromosome == "chr5" )
ENS_keys = ENS_keys[grepl(ENS_keys$gene_id, pattern="ENSG00000113356") ,]
# Get the interactions with POLR3G
interactions_polr3g = lapply(pcHIC, function(x) x[grepl(x$locus, pattern = "ENSG00000113356"),])

location_polr3g = lapply(interactions_polr3g, function(x) add_cCRE_location(x, cCRE_keys))
BENGI_polr3g = na.omit(select(ldply(location_polr3g), cCRE, chromosome))
print(BENGI_polr3g)

#save(BENGI_polr3g , file = "./data/BENGI_interactions_POLR3G.RData")


################################################################################
# RAW bengi polr3g interactions
load(file.path(data.dir,"/raw_BENGI_interactions_POLR3G.RData"))
raw_bengi_polr3g_df = ldply(raw_bengi_polr3g)
interactions = lapply(raw_bengi_polr3g, function(x) x$polr3g_interaction)


# intersection of GM and CD34 results 
library(UpSetR)
upset(fromList(interactions),
      sets.bar.color = "#56B4E9",
      order.by = "freq", #mb.ratio = c(0.4,0.6),
      empty.intersections = NULL, set_size.show = TRUE,
      set_size.scale_max = 30)


# Unique vs common
freq = data.frame(table(raw_bengi_polr3g_df$polr3g_interaction))
unique = as.vector(freq[freq$Freq == 1,]$Var1)
not_unique = as.vector(freq[freq$Freq >= 2,]$Var1)


# counting interactions
find_unique <- function(df, unique){
  tmp = ifelse(df$polr3g_interaction %in% unique, yes = TRUE, no = FALSE)
  return(data.frame("interactions" = c(sum(tmp==TRUE), sum(tmp==FALSE)),
                    "condition" = c("unique", "common"))) }

interactions_unique = lapply(raw_bengi_polr3g, 
                             function(x) find_unique(x, unique))

interactions_unique = ldply(interactions_unique)

p <- ggplot(data = interactions_unique, aes(x = reorder(.id,interactions), y = interactions)) +
  geom_col(aes(fill = condition), width = 0.7)+
  geom_text(data=subset(interactions_unique, interactions!=0),
            position = position_stack(vjust = 0.5),
            size = 4,
            aes(y = , label = interactions, group =condition), color = "white")+
  scale_fill_manual(values = c("#00539CFF", "orangered"))+
  theme_classic() +  
  xlab("Cell line") + ylab("Contact number") + 
  ggtitle("Number of statistically significant interactions with POLR3G (BENGI dataset)") 
p


pie_df = data.frame("value" = c(length(unique), length(not_unique)),
                    "condition" = c("unique","common"))


ggplot(pie_df, aes(x="", y=value, fill=condition)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_manual(values = c("#00539CFF", "orangered")) + 
  geom_text(
    position = position_stack(vjust = 0.5), 
    size = 5, 
    aes(y = ,label = value, group =condition), color = "white")+
  ggtitle("Number of statistically significant interactions \n with POLR3G (BENGI dataset)") 





################################################################################
# convert IDs to genomic locations
# Poubelle :
pchic_converted = lapply(pcHIC, function(x) add_ID_location(x, cCRE_keys, ENS_keys))


gm = pcHIC$GM12878
test = gm[grepl(gm$locus, pattern="ENSG00000113356"),]







# Converting data to the adequate format
pchic_converted = lapply(pchic_converted, function(x) dplyr::rename(x,
    frag1 = cCRE_location,
    frag2 = locus_location
  ))

interactions = lapply(pchic_converted, function(x) hic_format(x))


# select the possible interactions with POLR3G
polr3g_exact = GRanges(
  seqnames = Rle("chr5"),
  IRanges(89768410, 89777404))
# No fragment overlapping with this POLR3G location
interactions_polr3g = lapply(interactions, function(x) subset_by_overlap(ROI = polr3g_exact, df = x, remove_col = FALSE))
print(ldply(interactions_polr3g))





"ENSG00000113356"

gm = pcHIC$GM12878

filter(gm, locus == "ENSG00000113356" )

ENS_keys[grepl(ENS_keys$gene_id, pattern="ENSG00000113356") ,]

