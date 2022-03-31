################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : get_files_list.r                                                  #                                                        
# > Fonction : Create a .Rdata file containing a list of dataframe             #
# > of hi-c file (.po.txt)  + filtering on pvalues                             #
# @ COLAJANNI Antonin                                                          #
################################################################################


# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/interaction/pcHiC/3div/")


# I M P O R T
library(stringr)


files = list.files(data.dir,full.names = TRUE)
File_list_filtered0.05 = list()
File_list_filtered0.001 = list()

for (f in files) {
  filename = str_replace(f, ".po.txt.zip", "" )
  filename = str_replace(filename, "./data/phic/", "")
  if (filename == "rdata"){
    next
  }
  
  txt_file = paste(filename,".po.txt", sep = "")
  data = read.table(unz(f,txt_file),header = TRUE)
  
  filtered_data_0.05 = data[data$X.log10.result. > -log10(0.05),]
  filtered_data_0.001 = data[data$X.log10.result. > -log10(0.001),]
  
  File_list_filtered0.05[[filename]] = filtered_data_0.05
  File_list_filtered0.001[[filename]] = filtered_data_0.001
  
}

#rm(f,filename, files, txt_file, data)
#save.image(file = "./data/phic/rdata/3div_po_txt.Rdata")
#H1 = File_list$H1
#rm(File_list)
#save.image(file = "./data/phic/rdata/H1_po_txt.Rdata")


save(File_list_filtered0.05, file = "./data/interaction/rdata/3div_0.05.Rdata")
save(File_list_filtered0.001, file = "./data/interaction/rdata/3div_0.001.Rdata")



