################################################################################                                                                            #
# > February 2022                                                              #                                                
# > Script : shiny test                                                        #                                                        
# > Function : comprendre le fonctionnement de Rshiny                          #        
# @ COLAJANNI Antonin                                                          #
################################################################################

library(shiny)
library(ggplot2)

dataset <- diamonds

fluidPage(
  
  titlePanel(""),
  
  sidebarPanel(
    
    sliderInput('Window Size (kb)', 'Window Size (kb)', min=1, max=20e+3,
                value=min(0, 20e+3), step=10, round=0),
    
    sliderInput('Common interaction threhold', 'Common interaction threhold', min=1, max=40,
                value=min(0, 40), step=1, round=0),  
    
    
    selectInput('Choose Gene', 'Gene', c("POLR3G")),
    
    sliderInput('Distance to bait(kb)', 'Distance to bait (kb)', min=1, max=10000,
      value=min(0, 10000), step=10, round=0),  

    
    radioButtons("Cell lines", "Cell lines",choices = c("H1","GM12878","etc"))
    
  ),
    
  mainPanel(
    plotOutput('plot')
  )
)
