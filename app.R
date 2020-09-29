#### Generate Gene Plots 

library(Seurat)
library(dplyr)
library(shiny)


# reads data in from Seurat3 objects
# to be shared via our shiny app server
# the app.R and the seuratObject.RDS  should be copied to  www/shiny/bioinformatician/project_dir


scefilt.toplot <- readRDS( "PATH_TO_SEURAT_OBJECT.RDS")

allGenes <- rownames(scefilt.toplot)

################################################################
## this code chunk splits the page in panels 
## you can change the width to add more panels per row
## right not I have two plots of width 4 
################################################################

ui <- fluidPage(
  # This is the main title for the app
  
  titlePanel ( "TESTING"),
  # this creates a box where scientist can input their gene of interest
  # 
  # textInput(inputId = "Gene",
  #           label = "Gene Symbol:",
  #           value = "Myl3" ),
  # 
  ## this allows for entering multiple genes 
  selectizeInput(inputId = "genes",
                 label = "Enter/Select Genes to plot:",
                 choices = c(as.vector(sort(unique(allGenes)))),
                 selected = c("Myl3","Myl2", "Myl7"),
                 multiple = TRUE,
                 options = NULL ),
  # Upload file
  fileInput(
    inputId = "File1",
    label = "Upload File to use as gene signature",
    multiple = FALSE,
    accept = c("text/csv",
               "text/comma-separated-values,text/plain",
               ".csv"),
    width = NULL,
    buttonLabel = "Browse...",
    placeholder = "No file selected"
  ),
  # give list a name to appear in the Feature Plot
  textInput(inputId = "GeneListName",
            label = "Name Gene Signature",
            value = "" ),
  
  # Top 2 panels
  fluidRow(
    # left panel is a static view of UMAP generated for the analysis
    column(5,
           "TSNE Plot",
           plotOutput('tsne')
    ),
    ## Panel for features Plot, responds to change in input Gene box
    column(6,
           "Features Plot",
           downloadButton("downloadPlot", "Download"),
           plotOutput('plot1')
    )
  ),
  # middle panel
  fluidRow(
    ## Panel for Violin Plot, responds to change in input Gene box
    column( 8,
            "Violin Plot",
            downloadButton("downloadVPlot", "Download"),
            plotOutput('plot2')
    ),
    column(4,
           "Features Plot for Gene Signatures",
           downloadButton("downloadPlot3", "Download"),
           plotOutput('plot3')
    )
  ),
  # bottom panel
  fluidRow(
    ## Panel for DotPlot
    column( 8,
            "DotPlot",
            downloadButton("downloadDPlot", "Download"),
            plotOutput('dotplot'))
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
  
  geneplot <- reactive ({ 
    
    FeaturePlot(object = scefilt.toplot, 
                #features = input$Gene)
                features = input$genes)
    
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
             #features = input$Gene)
             features = input$genes,
             ncol=2)
    
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
 
 
 
 ### Features Plot for Gene signatures
 
 
 # # Read in file with list of genes
 # data1 <- reactive({
 #   req(input$File1)
 #      inFile <- input$File1
 #      read.csv(inFile$datapath, stringsAsFactors = FALSE)  
 #   
 #   
 # })
 # create plot

 genesignplot <- reactive ({ 
   req(input$File1)
   inFile <- input$File1
   data1 <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
   
   scefilt.toplot <- AddModuleScore(object= scefilt.toplot,
                                    features=  data1,
                                    name= input$GeneListName,
                                    ctrl= length(data1))
   
   FeaturePlot(object = scefilt.toplot, 
               #features = input$Gene)
               features = paste0(input$GeneListName,"1"))   
 }) 

 
 output$plot3 <- renderPlot({
   
   print(genesignplot())
   
 })
 # download feature 
 output$downloadPlot3 <- downloadHandler(
   filename = function() { paste(input$GeneListName, '.png', sep='') },
   content = function(file) {
     ggsave(file,genesignplot(),device = "png")
   }
 )
 
 ## Features plot responds to multiple gene entry
 # 
 # # create plot
 Dotplot <- reactive ({
   DotPlot(object = scefilt.toplot,
               features = input$genes)
           #cols = c("blue", "red", "green"), 
          # split.by = "genotype")  
             


 })
 # 
 # # print plot
 output$dotplot <- renderPlot({
   print(Dotplot())
 })
 # 
 # # download feature
 # output$downloadPlot <- downloadHandler(
 #   filename = function() { paste(input$Gene, '.png', sep='') },
 #   content = function(file) {
 #     ggsave(file,Dotplot(),device = "png")
 #   }
 # )
 # 
  
}

shinyApp(ui=ui, server = server)
