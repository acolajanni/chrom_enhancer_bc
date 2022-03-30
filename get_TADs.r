# > February 2022                                                              #                                                
# > Script : get_TADs                                                          #                                                        
# > Fonction : extract bed files for UCSC vizualisation of TADs                #                        
# @ COLAJANNI Antonin                                                          #
################################################################################


# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/TADs_ENCODE")


# I M P O R T
library(stringr)
library(dplyr)
source("./scripts/fun_interactions.r")
# Retrieve data
metaData = read.delim(file.path(data.dir,"metadata.tsv"))
metaData = filter(metaData, File.assembly == "hg19") 

# Filter nested and non nested TADs
nested_TAD = filter(metaData, grepl("nested", metaData$Output.type) )
TAD = filter(metaData, !grepl("nested", metaData$Output.type) )

# Get list of dataframe containing nested TADs
File_list_nested = list()
for (accession_num in nested_TAD$File.accession){
  
  bed_file = paste0(accession_num,".bed.gz")
  bed_location = file.path(data.dir,bed_file)
  table = read.table(gzfile(bed_location), col.names = c(
    "chr","start",'end',"name","score"))
  
  name = filter(nested_TAD, File.accession==accession_num)$Biosample.term.name
  File_list_nested[[name]] = table
}

# Get list of dataframe containing TADs
File_list_TAD = list()
for (accession_num in TAD$File.accession){
  
  bed_file = paste0(accession_num,".bed.gz")
  bed_location = file.path(data.dir,bed_file)
  table = read.table(gzfile(bed_location), col.names = c(
    "chr","start",'end',"name","score"))
  
  name = filter(TAD, File.accession==accession_num)$Biosample.term.name
  File_list_TAD[[name]] = table
}


TADs = lapply(File_list_TAD, function(x) 
                  x %>% filter(chr == "chr5"
                  ) %>% select(-name, -score) )

encode_TADs = TADs
save(encode_TADs, file = "./data/interaction/rdata/ENCODE_TAD.RData")


############################################################################

### TADs 3D genome

data.dir = file.path(main.dir,"data/interaction/3D_genome_TAD/hg19.TADs/")
files = list.files(data.dir)

TADs = list()
for (f in files){
  file = file.path(data.dir, f)
  table = read.delim(file, header = TRUE, col.names = c("chr","start","end"))
  
  # Filter on chr5 to save space
  table = filter(table, chr == "chr5")
  filename = str_replace(f, "_TADs.txt", "")
  TADs[[filename]] = table
}

genome_TADs = TADs

save(genome_TADs, file = "./data/interaction/rdata/3D_genome_TADs.RData")





############################################################################
### loops 3D genome

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/3D_genome_loop")


# I M P O R T
library(stringr)
library(dplyr)


files = list.files(data.dir)



loops = list()
for (f in files){
  file = file.path(data.dir, f)
  table = read.table(file, header = FALSE) 
  table = table[,c(1:6)]
  colnames(table) = c("sourceChr","sourceStart","sourceEnd","targetChr","targetStart","targetEnd")
  
  # Filter on chr5 to save space
  table = filter(table,  sourceChr == "chr5")#, sourceStart > 87000000, sourceStart < 92000000)
  filename = str_replace(f, ".hg19.peakachu-merged.loops", "")
  loops[[filename]] = table
}



bigInteract_loops = ldply(loops, .id=NULL)
bigInteract_loops$name = "loops"
bigInteract_loops = bigInteract_loops %>% group_by_all %>% count

names(bigInteract_loops)[names(bigInteract_loops) == 'freq'] <- 'score'


bigInteract_loops$chr = bigInteract_loops$sourceChr
bigInteract_loops$start = bigInteract_loops$sourceStart
bigInteract_loops$end = bigInteract_loops$targetEnd
bigInteract_loops$value = bigInteract_loops$score
bigInteract_loops$exp = "."
bigInteract_loops$color = rgb_to_1col(num_to_colors(bigInteract_loops$score))$color
bigInteract_loops$sourceName = "."
bigInteract_loops$targetName = "."
bigInteract_loops$sourceStrand = "."
bigInteract_loops$targetStrand = "."


bigInteract_loops = select(bigInteract_loops, 
                           chr,start,end,name,score,value,exp,color,
                           sourceChr,sourceStart,sourceEnd,sourceName,sourceStrand,
                           targetChr,targetStart,targetEnd,targetName,targetStrand)



write.table( bigInteract_loops, "./data/interaction/bed_file/loops_biginteract.bed", sep="\t", quote = F, row.names = F, col.names = F )




############################################################################
TADs = ldply(File_list_TAD) %>% 
  select(chr,start,end,#.id,score
         ) %>%
  filter(chr == "chr5")

nested_TADs = ldply(File_list_nested) %>% 
  select(chr,start,end,#.id,score
  ) %>%
  filter(chr == "chr5")

freq = data.frame(table(TADs)) %>%
  filter(Freq > 0)



nested_TADs = unique(nested_TADs)

write.table( nested_TADs, "./data/interaction/bed_file/encode_nested_TAD.bed", sep="\t", quote = F, row.names = F, col.names = F )




T47d = list("nested_TAD"=File_list_nested$T47D, "TAD" = File_list_TAD$T47D)
T47d = ldply(T47d)

test_t47d = T47d %>% 
  select(chr,start,end,.id) %>%
  filter(chr == "chr5")




write.table( test_t47d, "./data/interaction/bed_file/T47D_test.bed", sep="\t", quote = F, row.names = F, col.names = F )






    