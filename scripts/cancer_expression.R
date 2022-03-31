# > February 2022                                                              #                                                
# > Script : cancer_expression                                                 #                                                        
# > Fonction : Utilisation du package depmap pour récupérer l'expression       # 
# > de gènes dans les lignées cancers                                          #                        
# @ COLAJANNI Antonin                                                          #
################################################################################

# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/expression/")


# I M P O R T
#library(stringr)
#library(dplyr)
#source("./scripts/fun_interactions.r")

#BiocManager::install("depmap")
library("depmap")
library("ExperimentHub")

eh = ExperimentHub()
query(eh, "depmap")


## CCLE
#TPM <- eh[["EH2264"]]

polr3g = filter(TPM, gene_name == "POLR3G")




cancer_cell = c("NCI-H460","SJCRH30", "RPMI7951", "T47D", "SK-MEL-5",
                "CAKI2", "A549", "SK-N-MC", "Panc1", "ACHN", 
                "LNCaP", "G401", "DLD1", "K562","KBM7",
                "NCIH460","SKMEL5","SKNMC","T470","SKNDZ")
                
cancer_cell = c("MDAMB231","LNCAPCLONEFGC","SJRH30","MDAMB157","MDAMB436",
                "MDAMB453","MDAMB468","A549","G401","ACHN",
                "K562","NCIH460","RPMI7951","SKMEL5","SKNMC",
                "T47D","SKNDZ","DLD1","CAKI2","KBM7")



other_cell = c("H1","AD2","IMR90","HUVEC","HMEC","STL011","LG1",
               "P01","NHEK","PX1","STL002","BL","SB2","GM12878",
               "PA2","PX1","STL003","RV3")


cancer_data = TPM[grepl(pattern = paste(cancer_cell,collapse="|") , TPM$cell_line),]
#cancer_data = filter(cancer_data, expression > 0)
save(cancer_data, file = "./data/expression/rdata/cancer_expression.RData")

other_data = TPM[grepl(pattern = paste(other_cell,collapse="|") , TPM$cell_line),]
  

################################################################################
# W O R K I N G  D I R E C T O R Y
main.dir = "/shared/projects/chrom_enhancer_bc"
setwd(main.dir)

data.dir = file.path(main.dir,"data/expression/rnaseq")


# I M P O R T
CCLE_expr = read.table(file.path(data.dir,"./DepMap/CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct"), header = TRUE, skip=2)
#CCLE_probes = read.table(file.path(data.dir,"./DepMap/probeMap_gencode.v19.annotation.gene.probemap"), header = TRUE )


Roadmap_expr = read.table(file.path(data.dir,"./RoadMap/57epigenomes.N.pc"), header = TRUE)
Roadmap_cell_line = read.table(file.path(data.dir,"./RoadMap/EG.name.txt"), header = FALSE, col.names = c("ID","cell_line"))
Roadmap_ensembl = read.table(file.path(data.dir,"./RoadMap/Ensembl_v65.Gencode_v10.ENSG.gene_info"), 
                             header = FALSE, col.names = c("ID","chr","start","end","transcription","type","gene_name","hgnc") )





# Select cancer cell line of interest
cancer_cell = c("MDAMB231","MDAMB436","HCC1937", "HCC1806","HCC1954","MDAMB468",# TNBC
                "T47D","MCF7","LNCAPCLONEFGC","SJRH30","A549","G401",
                "ACHN","K562","NCIH460","RPMI7951","SKMEL5","SKNMC",
                "SKNDZ","DLD1","CAKI2","KBM7", "MCF10A",
                "Name", "Description")

CCLE_expr = CCLE_expr[,grepl(pattern = paste(cancer_cell,collapse="|") , names(CCLE_expr))]


# Select non-cancer cell line of interest ==> ALL of them
#other_cell = c("H1","hESC", "Aorta", "GM12878", "HMEC",
#               "HUVEC", "liver", "lung", "NHEK", "pancreas",
#               "spleen", "thymus", "ventricule", "breast", "pancreatic",
#               "gastric","psoas", "atrium", "colon", "thymus","spleen",
#               "HUVEC")

#Roadmap_cell_line = Roadmap_cell_line[grepl(paste(toupper(other_cell), collapse="|"), toupper(Roadmap_cell_line$cell_line)),]
#Roadmap_expr = Roadmap_expr[,names(Roadmap_expr) %in% c(Roadmap_cell_line$ID,"gene_id")]

# Replace colnames of roadmap expression matrixes ( E001 ==> H1)
tmp = data.frame(t(Roadmap_cell_line))
colnames(tmp) = tmp[1,]
tmp = tmp[-1,]
tmp$gene_id = "gene_id"
row.names(tmp) = NULL

Roadmap_expr = rbind(Roadmap_expr,tmp)
colnames(Roadmap_expr) = Roadmap_expr[nrow(Roadmap_expr),]
Roadmap_expr = Roadmap_expr[-nrow(Roadmap_expr),]

# Add gene annotation in Roadmap ENS ==> gene_name
Roadmap_expr = merge(Roadmap_expr, Roadmap_ensembl[,c("ID","gene_name")] ,by.x = "gene_id", by.y = "ID", all.x = TRUE)


# Get keys to differenciate cancer cells from non cancer ones
cancer = names(CCLE_expr)[! names(CCLE_expr) %in% c("Name","Description")] 
other = names(Roadmap_expr)[! names(Roadmap_expr) %in% c("gene_id","gene_name")]

cell_keys = data.frame("cell_line" = c(cancer,other) ,
                       "state" =  c(rep("cancer",length(cancer)) , rep("other", length(other)) ) )

# Merging cancer and not cancer dataframe of expession
expr = merge(CCLE_expr, Roadmap_expr, by.x = "Description", by.y = "gene_name")
colnames(expr)[colnames(expr) == "Description"] = "gene_name"
colnames(expr)[colnames(expr) == "Name"] = "ENSEMBL_id"

# Remove unnecessary columns
cell_lines = names(expr)[! names(expr) %in% c("Name","Description", "gene_id")]
expr = select(expr, c("gene_name","ENSEMBL_id",cell_lines))

save(expr, file = "./data/expression/rdata/cancer_healthy_expression.RData")


cell_keys$cell_line[,cell_keys$state == "other"]





