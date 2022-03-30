################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : weighted_interactions                                             #                                                        
# > Fonction : Apply a weight on interactions based on their frequence         #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")


# I M P O R T
source("./scripts/fun_interactions.r")
library(UpSetR)
load(file.path(data.dir,"3div_0.05.Rdata"))

### Chromosome 5: 
#89,767,565-89,810,370
# good
#89,758,361-89,767,753
# ce doc
#89,770,680-89,777,400

# LOADING DATA 
data = File_list_filtered0.05
rm(File_list_filtered0.05)
data_05 = lapply(data, hic_format )
data_01 = lapply(data_05, function(x) filter_pvalue(df=x, seuil = 2))

# RETRIEVING polr3g location
load(file.path(data.dir,"/promoterBaitsHuman.RData"))
polr3g = baitGR[baitGR$SYMBOL == "POLR3G"]

polr3g = GRanges(
  seqnames = Rle("chr5"),
  IRanges(89768380, 89777410))

polr3g_exact = GRanges(
  seqnames = Rle("chr5"),
  IRanges(89768410, 89777404))

rs10474352 = GRanges(
  seqnames = Rle("chr5"),
  IRanges(90680578, 90680579))


# interactions lists 
interactions = lapply(data_05, function(x) subset_by_overlap(ROI = polr3g_exact, df = x))
interaction_polr3g = lapply(interactions, function(x) x$interaction_polr3g)
# as dataframe
interactions_df = ldply(interactions)
# listing of Cell_lines
cell_lines = unique(interactions_df$.id)

upset(fromList(interaction_polr3g),
      sets = cell_lines, sets.bar.color = "#56B4E9",
      order.by = "degree", mb.ratio = c(0.4,0.6),
      empty.intersections = NULL, set_size.show = TRUE,
      set_size.scale_max = 50, nintersects = NA)




# All the interactions that are shared 
freq = data.frame(table(interactions_df$interaction_polr3g))
unique = as.vector(freq[freq$Freq == 1,]$Var1)
not_unique = as.vector(freq[freq$Freq >= 2,]$Var1)
  
  
# counting unique and common interactions
find_unique <- function(df, unique){
  tmp = ifelse(df$interaction_polr3g %in% unique, yes = TRUE, no = FALSE)
  return(data.frame("interactions" = c(sum(tmp==TRUE), sum(tmp==FALSE)),
                    "condition" = c("unique", "common"))) }

interactions_unique = lapply(interactions, 
                             function(x) find_unique(x, unique))

interactions_unique = ldply(interactions_unique)


p <- ggplot(data = interactions_unique, aes(x = reorder(.id,interactions), y = interactions)) +
  geom_col(aes(fill = condition), width = 0.7)+
  geom_text(
    position = "stack", size = 4,
    hjust=1.5,
    aes(y = , label = interactions, group =condition), color = "white")+
    scale_fill_manual(values = c("#00539CFF", "orangered"))+
    theme_classic() + coord_flip() +  
    xlab("Cell line") + ylab("Contact number")  
p 



#frequence des interactions :
polr3g_3div = data.frame(xtabs(~interaction_polr3g, interactions_df))



# Nb interactions / freq moyenne
count_freq = lapply(interactions, function(x) x$freq)
count_freq = lapply(count_freq, function(x) 
  return(data.frame("mean_freq" = mean(x), "interaction_number" = length(x))))
count_freq = ldply(count_freq)
row.names(count_freq) = count_freq$.id
count_freq$.id = NULL

# Non normal
shapiro.test(count_freq$mean_freq)

# non significatif
cor.test(~mean_freq+interaction_number, count_freq, method = "spearman")


ggplot(count_freq, aes(x=mean_freq, y=interaction_number)) +
  geom_point(color = "red") + 
  ggrepel::geom_text_repel(aes(label = rownames(count_freq)),
                 size = 3.5) +
  theme_classic() +
  xlab("Mean Frequency / Capturability") + ylab("Interaction number")











