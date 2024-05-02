processAndPlotGSEA <- function(fgseaResults, geneSetList, conditionDeframed, inputFilePath) {
  # Extract the base name of the input file without the path and extension
  baseName <- tools::file_path_sans_ext(basename(inputFilePath))
  
  # Collapse pathways
  collapsePath <- collapsePathways(fgseaResults[order(pval)], geneSetList, conditionDeframed)
  
  # Select main pathways
  mainPath <- fgseaResults[pathway %in% collapsePath$mainPathways][order(-NES), pathway]
  
  # Select top 20 pathways
  top20 <- fgseaResults[order(pval)]
  top20 <- top20[1:20,]
  top20_GOs <- top20$pathway
  
  # Plot GSEA tables for main pathways and save with a unique name
  plotGseaTable(geneSetList[mainPath], conditionDeframed, fgseaResults, gseaParam = 0.5)
  ggsave(paste0("figures/", baseName, "_mainPathwaysPlot.png"), width = 12, bg = "white")
  
  # Plot GSEA tables for top 20 pathways and save with a unique name
  plotGseaTable(geneSetList[top20_GOs], conditionDeframed, fgseaResults, gseaParam = 0.5)
  ggsave(paste0("figures/", baseName, "_top20PathwaysPlot.png"), width = 12, bg = "white")
  
  # Return main and top 20 pathways
  list(mainPathways = mainPath, top20Pathways = top20_GOs)
}
