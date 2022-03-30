################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : BENGIvs3div                                                       #                                                        
# > Function : comparing polr3g interactions between 3div and BENGI            #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")


# I M P O R T
source("./scripts/fun_interactions.r")
load(file.path(data.dir,"interactions_polr3g.Rdata"))
load(file.path(data.dir,"/raw_BENGI_interactions_POLR3G.RData"))
raw_bengi_polr3g = ldply(raw_bengi_polr3g, .id = NULL)


## For p<0.01
load(file.path(data.dir,"interactions_polr3g_01.Rdata"))
polr3g_3div = polr3g_3div_01
####

range_3div = GRanges(
  seqnames=Rle("chr5"),
  IRanges(polr3g_3div$interaction_polr3g)
)

range_bengi = GRanges(
  seqnames = Rle("chr5"),
  IRanges(raw_bengi_polr3g$polr3g_interaction)
)

subsetByOverlaps(range_3div, range_bengi)


vir = viridis(15,direction = -1)
colors = vir[cut(polr3g_3div$Freq,15)] 


show_col(colors)

inf = seq(88000000,91000000,by = (91000000-88000000)/15)
sup = seq(88190000,91000000,by = (91000000-88000000)/15)
bornes = paste0(inf[1:15],"-",sup[1:15])

################################################################################
# G V I Z

polr3g_exact = GRanges(
  seqnames = Rle("chr5"),
  IRanges(89768410, 89777404)
)

bengi = GRanges(raw_bengi_polr3g$polr3g_interaction, seqnames = Rle("chr5"))

div3 = GRanges(
  polr3g_3div$interaction_polr3g,
  seqnames = Rle("chr5"))  

color_scale = GRanges(bornes, seqnames = Rle("chr5"))

atrack_3div <- AnnotationTrack(div3, name = "3div", fill = colors, col = NULL)
atrack_bengi = AnnotationTrack(bengi, name = "BENGI",col = NULL, fill = "black")
atrack_scale = AnnotationTrack(color_scale, col = NULL, 
                               name = "color scale", 
                               fill =viridis(15, direction = -1),
                               showFeatureId=TRUE,
                               id = c(1:15)) 
gtrack <- GenomeAxisTrack(IRanges(start = 89768410 - 2000000, end = 90000000 ), chromosome = "chr5")



itrack <- IdeogramTrack(genome="hg19", 
                        chromosome="chr5", # specify chromosome in ucsc naming
                        from =89768410 - 2000000,
                        to=89777404 + 2000000)
itrack <- IdeogramTrack(genome="hg19", 
                        chromosome="chr5", # specify chromosome in ucsc naming
                        from =89768410 - 2000000,
                        to=90000000)

ht = HighlightTrack(trackList = list(atrack_scale, itrack, gtrack,atrack_3div, atrack_bengi),
                    start = 89768410, end = 89777404, name = "POLR3G_Bait")

plotTracks(list(ht),background.title = "darkblue")
