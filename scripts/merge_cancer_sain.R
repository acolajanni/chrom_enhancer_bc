################################################################################                                                                            #
# > March 2022                                                                 #                                                
# > Script : merge cancer sain                                                 #                                                        
# > Function : merge TADs 3D genome + ENCODE and make 2 categories :           #
# >            cancer vs other                                                 #        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/rdata/")

# I M P O R T
source("./scripts/fun_interactions.r")

load(file.path(data.dir,"3D_genome_TADs.RData"))
load(file.path(data.dir,"ENCODE_TAD.RData"))

# Cancer cell from ENCODE : all the cancer cell line of 3Dgenome are in ENCODE TADs
cancer_cell = c("NCI-H460","SJCRH30", "RPMI7951", "T47D", "SK-MEL-5",
                "Caki2", "A549", "SK-N-MC", "Panc1", "ACHN", 
                "LNCaP", "G401", "DLD1", "K562","KBM7",
                "NCIH460","SKMEL5","SKNMC","T470","SKNDZ") #Alternative spelling

# finding cancer data (= duplicate ?)
cancer_3d = names(genome_TADs)[grepl(pattern = paste(toupper(cancer_cell),collapse="|"),toupper(names(genome_TADs))) ]
cancer_encode = names(encode_TADs)[grepl(pattern = paste(cancer_cell,collapse="|") , names(encode_TADs))]

cancer_TADs = c(genome_TADs[cancer_3d],encode_TADs[cancer_encode])


# Other cell lines (= duplicate ?)
other_3d = names(genome_TADs)[! names(genome_TADs) %in% cancer_3d]
other_encode = names(encode_TADs)[! names(encode_TADs) %in% cancer_encode]

other_TADs = c(genome_TADs[other_3d],encode_TADs[other_encode])                



#######
cancer_TADs = ldply(cancer_TADs, .id=NULL)
cancer_TADs = filter(data.frame(table(cancer_TADs)),Freq > 0)
cancer_TADs$name = "cancer_TADs"
colnames(cancer_TADs) = c("chr","start","end","score","name")

cancer_TADs = full_bed_format(cancer_TADs,colorscale = "mako")
cancer_TADs = bed_to_bigInteract(cancer_TADs)


other_TADs = ldply(other_TADs, .id=NULL)
other_TADs = filter(data.frame(table(other_TADs)),Freq > 0)
other_TADs$name = "not_cancer_TADs"
colnames(other_TADs) = c("chr","start","end","score","name")

other_TADs = full_bed_format(other_TADs,colorscale = "viridis")
other_TADs = bed_to_bigInteract(other_TADs)



save(cancer_TADs, file = "./data/interaction/rdata/cancer_TADs.RData")
save(other_TADs, file = "./data/interaction/rdata/other_TADs.RData")


write.table( cancer_TADs, "./data/interaction/bed_file/cancer_TADs_bigInteract.bed", sep="\t", quote = F, row.names = F, col.names = F )
write.table( other_TADs, "./data/interaction/bed_file/other_TADs_bigInteract.bed", sep="\t", quote = F, row.names = F, col.names = F )





