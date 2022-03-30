################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : functions_vizu                                                    #                                                        
# > Fonction : Fonctions pour la visu                                          #        
# @ COLAJANNI Antonin                                                          #
################################################################################


library(tidyverse)
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)



#' Extract one feature from the df list 
#'
#' @param data list of dataframes 
#' @param feature name of one column of the dataframes
#'
#' @return list of vector for each dataframe of the original list
extract_features_from_list <- function(data,feature = "dist"){
  return(lapply(data, "[", , feature)) }


#' Creates a dataframe with 2 columns (variable and value) from a list of vector
#'
#' @param data list of vector
#'
#' @return DataFrame
list_to_2colDF <- function(data){
  df = data.frame(do.call(cbind, data))
  return (melt(df))
}

#title = "Distribution of contact distances of the cell lines from 3div database (PCHi-C filtered at p < 0.05)"
hist_contact <- function(distances, title){
  ggplot(data = distances) +
    geom_histogram(aes(x = value), bins = 30 )+
    labs(title=title) +
    facet_wrap(~ variable)
}

#title = "Number of contact of the cell lines from 3div database (PCHi-C filtered at p < 0.05)"
contact_distrib <- function(interact, title){
  ggplot(data=interact, aes(x=reorder(variable, -value), y=value)) +
    geom_bar(stat="identity")+
    geom_text(aes(label=value), color="white",hjust=-0.1 , size=3.8, angle=-90) +
    scale_y_continuous(breaks = c(0,100000,300000), labels = comma) +
    labs(title=title,
         x = "Cell line", y="Contact number")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
}


