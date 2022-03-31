################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : to_bed                                                            #                                                        
# > Function : Create bed files with bengi and 3div                            #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")

# I M P O R T
source("./scripts/fun_interactions.r")
# 3DiV
load(file.path(data.dir,"interactions_polr3g.Rdata"))
# BENGI
load(file.path(data.dir,"/raw_BENGI_interactions_POLR3G.RData"))
raw_bengi_polr3g_df = ldply(raw_bengi_polr3g)
# TADs
load(file.path(data.dir,"/ENCODE_TAD.RData"))
TADs = ldply(TADs, .id=NULL)
bed_TADs = filter(data.frame(table(TADs)),Freq > 0)
bed_TADs$name = "TADs_ENCODE"
colnames(bed_TADs) = c("chr","start","end","score","name")

# cancer TADs
load(file.path(data.dir,"/cancer_TADs.RData"))
cancer_TADs = ldply(cancer_TADs)#, .id=NULL)
bed_TADs = filter(data.frame(table(TADs)),Freq > 0)
cancer_TADs$name = cancer_TADs$.id
cancer_TADs$score = 1000
cancer_TADs = cancer_TADs[,c("chr","start","end","score","name")]
cancer_TADs = full_bed_format(cancer_TADs,colorscale = "plasma")


# Split genomic location (chr5:XXXX-XXX to start XXX ; end XXX) (3DiV)
bed_3div = polr3g_3div %>% mutate(
  chr = "chr5",
  start = str_split_fixed(interaction_polr3g,'-', n=2)[,1],
  end = str_split_fixed(interaction_polr3g,'-', n=2)[,2],
  score = Freq,
  name = "3div"
) %>% select(-interaction_polr3g, -Freq)


# Compute a score value based on frequency 
bed_bengi = data.frame(table(raw_bengi_polr3g_df$polr3g_interaction ))
# Split genomic location (chr5:XXXX-XXX to start XXX ; end XXX) (BENGI)
bed_bengi = bed_bengi %>% mutate(
  chr = "chr5",
  start = str_split_fixed(Var1,'-', n=2)[,1],
  end = str_split_fixed(Var1,'-', n=2)[,2],
  score = Freq,
  name = "BENGI"
) %>% select(-Var1, -Freq)




# Translate score to colors and make it bed format (3DiV)
bed_3div_col = full_bed_format(bed_3div)

# Translate score to colors and make it bed format (bengi / Misfud)
bed_bengi_col = full_bed_format(bed_bengi)

# Translate score to colors and make it bed format (ENCODE TADs)
bed_TADs_col = full_bed_format(bed_TADs)



write.table( bed_3div_col, "./data/interaction/bed_file/3div_track_color.bed", sep="\t", quote = F, row.names = F, col.names = F )
write.table( bed_bengi_col, "./data/interaction/bed_file/bengi_track_color.bed", sep="\t", quote = F, row.names = F, col.names = F )
write.table( bed_TADs_col, "./data/interaction/bed_file/ENCODE_TADs_color.bed", sep="\t", quote = F, row.names = F, col.names = F )





write.table( cancer_TADs, "./data/interaction/bed_file/cancer_TADs_test.bed", sep="\t", quote = F, row.names = F, col.names = F )

