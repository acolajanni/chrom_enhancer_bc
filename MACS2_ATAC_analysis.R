################################################################################                                                                            #
# > March 2022                                                                                                                 
# > Script : MACS2_ATAC_analysis                                                                                                         
# > Function : figure to analyse the MACS treated bedgraphs           
# @ COLAJANNI Antonin                                                          
################################################################################

## ATAC_seq
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc/"
setwd(main.dir)

# I M P O R T
source("./scripts/fun_interactions.r")
data.dir = "./data/interaction/bedgraph_analysis"
chr5.dir = file.path (data.dir,"chr5_default/")

file_names = list.files(chr5.dir)


import_MACS2_output = function(filenames, directory){
  bed_list = list()
  columns = c("chrom","start","end","name","score","strand","signalValue","pvalue","qvalue","peak")
  for (file in filenames){
  
    file_path = paste0(directory,"/",file)
  
    tmp = readLines(file_path, n=2)
    if (length(tmp) == 1){ tmp = read.table(text = "", col.names = columns) }
    else {tmp = read.table(file_path,skip = 1, col.names = columns)}
    
      bed_list[[file]] = tmp
  }
  return(bed_list)
}


MACS_output = import_MACS2_output(filenames = file_names, directory = chr5.dir) 
MACS_output_filtered = lapply(MACS_output, function(x) filter(x, between(start, 84000000,94000000)))

counter_row = lapply(MACS_output_filtered, function(x) x = nrow(x))
counter_row = ldply(counter_row)


title ="aaaa"


ggplot(data = counter_row, aes(x = V1, y= reorder(.id,V1), fill=V1)) +
  geom_bar(stat="identity") +
  geom_text(data = subset(counter_row, V1 != 0), 
            color="black",hjust=-0.2, size=3.8, angle=0,
            aes(y = ,label = V1))+
  coord_cartesian(xlim = c(0, max(counter_row$V1)+0.1*max(counter_row$V1))) + 
  labs(title = title, x = "peak number chr5:84,000,000-chr5:94,000,000",
       y = "samples") +
  theme_minimal() 





ggplot(data=interact, aes(x=reorder(variable, -value), y=value)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=value), color="white",hjust=-0.1 , size=3.8, angle=-90) +
  scale_y_continuous(
    #breaks = c(0,100000,300000), 
    labels = comma
  ) +
  labs(title=title,
       x = "Cell line", y="Contact number")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))
