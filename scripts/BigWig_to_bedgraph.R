################################################################################                                                                            #
# > March 2022                                                                                                                 
# > Script : BigWig_to_bedgraph                                                                                                         
# > Function : convert bigwig file into file readable by UCSC genome browser : histogram of signal                 
# @ COLAJANNI Antonin                                                          
################################################################################

## ATAC_seq
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

# I M P O R T
source("./scripts/fun_interactions.r")
data.dir = file.path(main.dir,"data/interaction/ATAC")
save_file = file.path(main.dir,"data/interaction/bed_file/bigwig" )

files = list.files(file.path(data.dir,"K562"), pattern = ".bw|bigWig", full.names = TRUE, recursive = TRUE)

## K562
K562_filename = list.files(file.path(data.dir,"K562"), pattern = ".bw|.bigWig", full.names = F, recursive = TRUE)

K562_bed_filename =ifelse(str_detect(K562_filename,"K562"),
                      yes = str_remove_all(K562_filename,".bw|.bigWig|_nan"),
                      no = paste0(str_remove(K562_filename,".bw|.bigWig"),"_K562") )

K562_bed_filename =ifelse(str_detect(K562_bed_filename,"ATAC"),
                          yes = K562_bed_filename,
                          no = paste0(K562_bed_filename,"_ATAC-seq"))



K562_bed_filenameK562_bed_filename = str_replace(string = K562_bed_filename, pattern = "/", replacement = "_")


count = 1  
ATAC_K562 = list()

paths = files[grepl(paste0(K562_filename,collapse = "|"), files)]

BigWig_paths = data.frame("path" = paths,
                          "names" = K562_bed_filename)


for(f in paths ){
  name = K562_bed_filename[count]
  path = file.path(save_file, paste0(name,".bedGraph")   )
  
  
  ATAC_K562[[ K562_bed_filename[count] ]] = BigWig_to_bedGraph(bigwig_path = f, 
                                                               file_path = path,
                                                               name = name, 
                                                               apply.log = FALSE)
  
  count = count+1
}


save(ATAC_K562, file = "./data/interaction/rdata/GROseq/ATAC_K562.RData")

## H1
files = list.files(file.path(data.dir), pattern = ".bw|bigWig", full.names = TRUE, recursive = TRUE)
H1_filename = list.files(file.path(data.dir,"H1"), pattern = ".bw|.bigWig", full.names = F, recursive = TRUE)

H1_bed_filename =str_replace(H1_filename,"_openChrmtn.bw","")
H1_bed_filename = str_replace(string = H1_bed_filename, pattern = "/", replacement = "_")

H1_bed_filename = ifelse(str_detect(H1_bed_filename, "ATAC"),
                         yes = H1_bed_filename,
                         no = paste0(H1_bed_filename, "_ATAC-seq") )


paths = files[grepl(paste0(H1_filename,collapse = "|"), files)]

BigWig_paths = rbind(BigWig_paths, data.frame("path" = paths, "names" = H1_bed_filename))

ATAC_H1 = multiple_bigWig_to_bedGraph(paths = paths,
                                      filenames = H1_bed_filename,
                                      save_file = save_file,
                                      apply.log = FALSE)


save(ATAC_H1, file = "./data/interaction/rdata/GROseq/ATAC_H1.RData")

## MB231
files = list.files(file.path(data.dir), pattern = ".bw|.bigwig", full.names = TRUE, recursive = TRUE)
MB231_filename = list.files(file.path(data.dir,"MDAMB231"), pattern = ".bw$|.bigwig$", full.names = F, recursive = TRUE)


MB231_bed_filename =str_replace(MB231_filename,".bw|.bigwig","")
MB231_bed_filename = str_replace(string = MB231_bed_filename, pattern = "/", replacement = "_")
MB231_bed_filename = ifelse(str_detect(MB231_bed_filename, "ATAC"),
                         yes = MB231_bed_filename,
                         no = paste0(MB231_bed_filename, "_ATAC-seq") )

paths = paste0(file.path(data.dir,"MDAMB231/"), MB231_filename)


BigWig_paths = rbind(BigWig_paths, data.frame("path" = paths, "names" = MB231_bed_filename))
write.table(BigWig_paths, file = "./data/interaction/ATAC/BigWig_paths.txt", quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)



ATAC_MB231 = multiple_bigWig_to_bedGraph(paths = paths,
                                      filenames = MB231_bed_filename,
                                      save_file = save_file,
                                      apply.log = FALSE)


save(ATAC_MB231, file = "./data/interaction/rdata/GROseq/ATAC_MB231.RData")











#bw = file.path("/shared/projects/chrom_enhancer_bc/data/interaction/ATAC/MDAMB231/GSE129646/GSE129646_MDAMB231_Par_atac.merge.nomt.bs5.RPKM.bw")
#gr = rtracklayer::import.bw(bw)

#track_type = "bedGraph"
#name = "test"
#descr = name
#color = "0,0,255"
#options = 'visibility=full yLineOnOff=on autoScale=on yLineMark="0.0" alwaysZero=on graphType=bar maxHeightPixels=128:75:11 windowingFunction=maximum smoothingWindow=off'

#df <- data.frame(chr=seqnames(gr),
#                 start=start(gr)-1,
#                 end=end(gr),
#                 score= gr$score)
#df = filter(df, chr == 5)



#first_line = paste0("track type=",track_type, " name=", name, " description=", descr, " color=",color, " ",options )
#headline = data.frame(chr = first_line,
#                      start = "",
#                      end = "",
#                      score = "")
  
#bedGraph = rbind(headline, df)
  


#write.table( bedGraph, "./data/interaction/bed_file/bigwig/test.bedGraph", sep="\t", quote = F, row.names = F, col.names = F )







#test3 = BigWig_to_bedGraph(bw, file_path = "./data/interaction/bed_file/bigwig/test.bedGraph",
#                   chr_filter = TRUE, chr_to_keep = "chr5", name = "test", color = "255,0,0")





############## GROseq
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

# I M P O R T
source("./scripts/fun_interactions.r")
chain_liftover = file.path("./docs/hg38ToHg19.over.chain")

data.dir = file.path(main.dir,"data/interaction/GROseq")


# Retrieving all gros seq bigwig files
files = list.files(file.path(data.dir), pattern = ".bw|bigWig", full.names = TRUE, recursive = TRUE)

# Removing files that are too far in the folder (not necessary : tar.xz for other GSE96859 files)
#files = files[unlist(lapply(files,function(x) return(length(str_extract_all(x, "/")[[1]]) < 10) ))]


# Creating subobject for different experiment
GSE96859 = files[grepl(pattern = "GSE96859", files)]
GSE158669 = files[grepl(pattern = "GSE158669", files)]
GSE64758 = files[grepl(pattern = "GSE64758", files)]


# H1 : GSE64758
GSE64758_plus = BigWig_to_bed(GSE64758[grepl("plus|Plus",GSE64758)])
GSE64758_plus = filter(GSE64758_plus, between(start,87500000,91500000))

GSE64758_minus = BigWig_to_bed(GSE64758[grepl("minus|Minus",GSE64758)])
GSE64758_minus = filter(GSE64758_minus, between(start,87500000,91500000))

H1_GSE64758 = merge_plus_minus_GROseq(plus_strand = GSE64758_plus,
                             minus_strand = GSE64758_minus, 
                             name = "H1_GSE64758_GRO-seq", color = "0,0,255")

rm(GSE64758_minus, GSE64758_plus)
save(H1_GSE64758, file = "./data/interaction/rdata/GROseq/H1_GSE64758.RData")
write.table(H1_GSE64758, file = "./data/interaction/bed_file/bigwig/H1_GSE64758_GRO-seq.bedGraph",quote = F, row.names = F, col.names = F )


# MB231 : GSE96859
GSE96859_plus = BigWig_to_bed(GSE96859[grepl("plus|Plus",GSE96859)], liftover = FALSE, lift_path = chain_liftover)
GSE96859_plus = filter(GSE96859_plus, between(start,87500000,91500000))


GSE96859_minus = BigWig_to_bed(GSE96859[grepl("minus|Minus",GSE96859)], liftover = FALSE, lift_path = chain_liftover)
GSE96859_minus = filter(GSE96859_minus, between(start,87500000,91500000))

MB231_GSE96859 = merge_plus_minus_GROseq(plus_strand = GSE96859_plus,
                             minus_strand = GSE96859_minus, 
                             name = "MB231_GSE96859_GRO-seq", color = "0,0,255")

write.table(MB231_GSE96859, file = paste0("./data/interaction/bed_file/bigwig/","MB231_GSE96859_untreated_GRO_seq",".bedGraph"),
            quote = F, row.names = F, col.names = F )


save(MB231_GSE96859, file = "./data/interaction/rdata/GROseq/MB231_GSE96859_GRO-seq.RData")
rm(GSE96859_plus, GSE96859_minus)

# MB231 : GSE158669
# 2 conditions
shNC = GSE158669[grepl("shNC",GSE158669)]
shBCDIN3D = GSE158669[grepl("shBCDIN3D",GSE158669)]

# list of paths sorted by condition and replicates
GSE158669_mb231 = list("shNC" =  list("rep1" = shNC[grepl(pattern = "Rep1", shNC)],
                                      "rep2" = shNC[grepl(pattern = "Rep2", shNC)]),
                       
                       "shBCDIN3D" = list("rep1" = shBCDIN3D[grepl(pattern = "Rep1", shBCDIN3D)],
                                          "rep2" = shBCDIN3D[grepl(pattern = "Rep2", shBCDIN3D)]))


# import and merge +/- strand GROseq data 
mb231_final = list()
for (cond in c("shNC","shBCDIN3D") ) {

  for (rep in c("rep1","rep2") ){
    message(cond,rep)
    filename = paste0("MB231_GSE158669_",cond,"_",rep,"_GRO-seq")
    
    paths = GSE158669_mb231[[cond]][[rep]]
    
    plus_strand = BigWig_to_bed(paths[grepl("plus|Plus",paths)], liftover = TRUE, lift_path = chain_liftover)
    plus_strand = filter(plus_strand, between(start,88000000,92000000))
    
    minus_strand = BigWig_to_bed(paths[grepl("minus|Minus",paths)], liftover = TRUE, lift_path = chain_liftover)
    minus_strand = filter(minus_strand, between(start,88000000,92000000))
    
    
    mb231_final[[cond]][[rep]] = merge_plus_minus_GROseq(plus_strand,minus_strand,
                            name = filename, color = "0,0,255")
    
    write.table(mb231_final[[cond]][[rep]], file = paste0("./data/interaction/bed_file/bigwig/",filename,".bedGraph"),
                quote = F, row.names = F, col.names = F )
  }
}
MB231_GSE158669 = mb231_final
rm(plus_strand,minus_strand,mb231_final)
save(MB231_GSE158669,  file = "./data/interaction/rdata/GROseq/MB231_GSE158669_GRO-seq.RData")
load("./data/interaction/rdata/GROseq/MB231_GSE158669_GRO-seq.RData")

# Now the bedgraphs file : 
bedgraph_file = dir(data.dir,full.names = TRUE,recursive = TRUE, pattern = ".bedGraph$")
HS = bedgraph_file[grepl(pattern = "HS_yhd", bedgraph_file)]
control = bedgraph_file[! bedgraph_file %in% HS]

#minus = read.table(file.path(data.dir,"K562-GROseq-C-minus.bedGraph"),skip = 1, col.names = c("chr","start","end","score"))
#plus = read.table(file.path(data.dir,"K562-GROseq-C-plus.bedGraph"),skip = 1, col.names = c("chr","start","end","score"))



K562_GSE66448 = list()
for (paths in list(HS,control) ) {
  colnames = c("chr","start","end","score")
  
  print(paths[grepl("plus|Plus",paths)])
  plus_strand = read.table(paths[grepl("plus|Plus",paths)], skip=1, col.names = colnames)
  plus_strand = filter(plus_strand, between(start,88000000,92000000))
  
  print(paths[grepl("minus|Minus",paths)])
  minus_strand = read.table(paths[grepl("minus|Minus",paths)], skip = 1, col.names = colnames)
  minus_strand = filter(minus_strand, between(start,88000000,92000000))
  
    
    
  if (str_detect(paths[1], "HS")) {cond = "HS"}
  else {cond = "control"}
  filename = paste0("K562_GSE66448_",cond,"_GRO-seq")
  
  
  K562_GSE66448[[cond]]= merge_plus_minus_GROseq(plus_strand,minus_strand,
                                        name = filename, color = "0,0,255")
  
  write.table(K562_GSE66448[[cond]], file = paste0("./data/interaction/bed_file/bigwig/",filename,".bedGraph"),
              quote = F, row.names = F, col.names = F )
  
}
save(K562_GSE66448, file = "./data/interaction/rdata/GROseq/K562_GSE66448.RData")
rm(plus_strand,minus_strand)


#total = merge_plus_minus_GROseq(plus,minus,"azertyuiop","255,0,0")
#write.table( total, "./data/interaction/bed_file/bigwig/total.bedGraph", sep="\t", quote = F, row.names = F, col.names = F )





############## GROseq_paths
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)
data.dir = file.path(main.dir,"data/interaction/GROseq")

# Retrieving all gros seq bigwig files
files = list.files(file.path(data.dir), pattern = ".bw|bigWig", full.names = TRUE, recursive = TRUE)

files_K562 = list.files(path = file.path(data.dir,"K562"), pattern = ".bedGraph$", recursive = T, full.names = T)
filenames_K562 = list.files(path = file.path(data.dir,"K562"), pattern = ".bedGraph$", recursive = T)

#files_K562
#filenames_K562


#files = c(files,files_K562)
filenames = str_remove(files, paste0(data.dir,"/") )
filenames_K562 = str_remove(files_K562, paste0(data.dir,"/K562/") )

filenames = str_remove(filenames, pattern = ".bw|.bigWig|.bedGraph")
filenames = str_replace_all(filenames, pattern= "/", replacement = "_")
paths = data.frame(files, filenames)


filenames_K562 = str_remove(filenames_K562, ".bw|.bigWig|.bedGraph")
filenames_K562 = str_replace_all(filenames_K562, pattern= "/", replacement = "_")



write.table(paths, file = "./data/interaction/GROseq/BigWig_paths.txt", quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

write.table(data.frame(files_K562,filenames_K562), file = "./data/interaction/GROseq/bedgraph_paths.txt", quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

