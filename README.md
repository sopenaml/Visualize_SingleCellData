# Visualize_SingleCellData

## Shiny app for visualization of single cell Data

This app has been writen to view single cell data analized with Seurat3. It plots an static UMAP with clusters coloured by cluster ID and two interactive plots; a FeaturesPlot and a VlnPlot for a single gene each time. Individual plots can be downloaded with the "Download" button. 

We have added 2 new features: 
  Plot multiple genes at the same time. 
  Plot gene signatures: Upload a file with a list of genes, (txt or csv format), use AddModuleScore and plot values on a FeaturePlot.
