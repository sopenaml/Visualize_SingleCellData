# Visualize_SingleCellData

## Shiny app for visualization of single cell Data

This app has been writen to view single cell data analized with Seurat3. It plots an static UMAP with clusters coloured by cluster ID and two interactive plots; a FeaturesPlot and a VlnPlot for a single gene each time. Individual plots can be downloaded with the "Download" button. 

We have added 2 new features: 
 * Plot multiple genes at the same time (We don't recommend to plot more than 4 genes, as plots will become too small).
 * Plot gene signatures: Upload a file with a list of genes, (txt or csv format), use AddModuleScore and plot values on a FeaturePlot.
 
Once the file is uploaded, the app will automatically calculate ModuleScore and generate a FeaturePlot. As default it will be called "1". A name can be added to the box "Name Gene Signature".
