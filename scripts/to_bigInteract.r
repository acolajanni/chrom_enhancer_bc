################################################################################                                                                            #
# > March 2022                                                                 #                                                
# > Script : to_bigInteract                                                    #                                                        
# > Function : Create a bigInteract file for TAD annotation                    #        
# @ COLAJANNI Antonin                                                          #
################################################################################


# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")

# I M P O R T
source("./scripts/fun_interactions.r")

load(file.path(data.dir,"3D_genome_TADs.RData"))

genome_TADs = ldply(genome_TADs, .id=NULL)
genome_TADs$name = "3Dgenome"
bed_TADs = filter(data.frame(table(genome_TADs)),Freq > 0)
colnames(bed_TADs) = c("chr","start","end","name","score")
bed_TADs$value = bed_TADs$score
bed_TADs$exp = "."
bed_TADs = cbind(bed_TADs, rgb_to_1col(num_to_colors(bed_TADs$score)) )

#bed_TADs$value = 1

bed_TADs$sourceChrom = bed_TADs$chr
bed_TADs$sourceStart = bed_TADs$start
bed_TADs$sourceEnd = bed_TADs$start
bed_TADs$sourceName = NA
bed_TADs$sourceStrand = "."

bed_TADs$targetChrom = bed_TADs$chr
bed_TADs$targetStart = bed_TADs$end
bed_TADs$targetEnd = bed_TADs$end
bed_TADs$targetName = NA
bed_TADs$targetStrand = "."






write.table( bed_TADs, "./data/interaction/bed_file/3Dgenome_TADs_bigInteract.bed", sep="\t", quote = F, row.names = F, col.names = F )


bed_only_TAD = select(bed_TADs, chr,start,end,name,score,sourceStrand,sourceStart,sourceEnd,color)

write.table( bed_only_TAD, "./data/interaction/bed_file/3Dgenome_TADs.bed", sep="\t", quote = F, row.names = F, col.names = F )
