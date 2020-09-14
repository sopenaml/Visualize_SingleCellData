#### Generate Gene Plots 

library(Seurat)
library(dplyr)
library(shiny)


# read data in
# it should be copied to  www/shiny/bioinformatician/project_dir


scefilt.toplot <- readRDS( "PATH_TO_SEURAT_OBJECT.RDS")


################################################################
## this code chunk splits the page in panels 
## you can change the width to add more panels per row
## right not I have two plots of width 4 
################################################################

ui <- fluidPage(
  # This is the main title for the app
  
  titlePanel ( "PROJECT ID"),
  # this creates a box where scientist can input their gene of interest
  
  textInput(inputId = "Gene",
            label = "Gene Symbol:",
            value = "Myl3" ),
  # Top 2 panels 
  fluidRow(
    # left panel is a static view of UMAP generated for the analysis
    column(4,
           "TSNE Plot",
           plotOutput('tsne')
           ),
    ## Panel for features Plot, responds to change in input Gene box
    column(4,
           "Features Plot",
           downloadButton("downloadPlot", "Download"),
           plotOutput('plot1')
        )
    ),
  # bottom panel
  fluidRow(
    ## Panel for Violin Plot, responds to change in input Gene box
    column( 4,
            "Violin Plot",
            downloadButton("downloadVPlot", "Download"),
            plotOutput('plot2')
           )
   )
)

##########################################################
## this is the code that reads data and generates plots ##
## for static plots you just create an output
## for plots that respond to an input you need two calls
## reactive: to get the data
## renderPlot: to generate the plot including the reactive data in the code
## you'll need a reative / render plot or renderPlot (only for static plots), for each panel above 
## and it will have to match names with the names given above 
##########################################################

server <- function(input, output) {
  
 # Plot static UMAP/TSNE only uses renderPlot
  
  output$tsne <- renderPlot({ 
    DimPlot(scefilt.toplot, reduction = "umap")
    
   
  }) 
  
## Features plot responds to input Gene box

# create plot
    geneplot <- reactive ({ 
      FeaturePlot(object = scefilt.toplot, 
      features = input$Gene)

   })
    
# print plot
  output$plot1 <- renderPlot({
    print(geneplot())
  })
  
# download feature
output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$Gene, '.png', sep='') },
    content = function(file) {
      ggsave(file,geneplot(),device = "png")
    }
  )

## Violin plot responds to input Gene box
# create plot
  Vlnplot <- reactive ({
     VlnPlot(object = scefilt.toplot, 
                        features = input$Gene )
})

# print plot
  output$plot2 <- renderPlot({
      
    print(Vlnplot())
    
  })
# download feature 
 output$downloadVPlot <- downloadHandler(
    filename = function() { paste(input$Gene, '.png', sep='') },
    content = function(file) {
      ggsave(file,Vlnplot(),device = "png")
    }
  )
 
    
  
}

shinyApp(ui=ui, server = server)
