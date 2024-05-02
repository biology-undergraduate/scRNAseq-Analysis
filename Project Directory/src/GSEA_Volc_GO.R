# TO RUN THIS SCRIPT, THE TERMINAL MUST BE USED. USE 'cd (directory path)" TO
# SET THE DIRECTORY AND THEN PASS THIS CODE INTO THE TERMINAL
# "Rscript src/GSEA_Volc_GO.r alpha_unstim" THIS WILL RUN THE SCRIPT USING
# THE 'alpha_unstim.csv" DATASET. WILL ONLY WORK USING THE CORRECT DIRECTORY
# FORMATTING

# Load necessary libraries
setwd("~/University/Fourth Year/RNAseq Directory")

source("functions/libraryLoader.R")
source("functions/GSEAprep.R")
source("functions/GSEAplotExport.R")
source("functions/GSEAtableExport.R")
source("functions/GOplotExport.R")
source("functions/volcanoPlotExport.R")

biocLibs <- c("fgsea","clusterProfiler","erichplot","org.Hs.eg.db")
BiocManager::install(biocLibs)

libs <- c("BiocManager","dplyr","Seurat","patchwork","tidyverse","ggplot2","tibble","msigdbr","stats","fastmatch","fgsea","clusterProfiler","erichplot","org.Hs.eg.db")

library_loader(libs)

# Accept command-line arguments for dataset name
args <- commandArgs(TRUE)
datasetName <- args[1]

# Construct file paths and export names dynamically
inputFilePath <- paste0("processedData/", datasetName, ".csv")
fgseaExportPath <- paste0("processedData/fgsea_", datasetName, ".csv")

# Preparing comparison for FGSEA
conditionDeframed <- prepare_data(inputFilePath)

# Obtaining gene lists and associated functions
humanGeneSets_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
geneSetList <- split(x = humanGeneSets_bp$gene_symbol, f = humanGeneSets_bp$gs_name)

# Ranking comparison
conditionRanked <- conditionDeframed[order(conditionDeframed, decreasing = TRUE)]
table(is.na(conditionRanked))

# Generating FGSEA
fgsea <- fgseaMultilevel(pathways = geneSetList, stats = conditionRanked, minSize = 15)
table(is.na(fgsea))

# Plotting GSEA of interest
plotDetails <- list(
  GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION = list(
    gene_set_name = "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
    title = paste(datasetName, "- Positive regulation of cytokine production")
  ),
  GOBP_REGULATION_OF_VIRAL_PROCESS = list(
    gene_set_name = "GOBP_REGULATION_OF_VIRAL_PROCESS",
    title = paste(datasetName, "- Regulation of viral process")
  )
)

generate_and_export_plots(conditionDeframed, geneSetList, plotDetails, paste0("figures/", datasetName))

# Multiplot
results <- processAndPlotGSEA(fgsea, geneSetList, conditionDeframed, inputFilePath)

mainPathways <- results$mainPathways
top20Pathways <- results$top20Pathways

# GO

GODetails <- paste(datasetName, "- GO")

generate_and_export_GO(inputFilePath,GODetails,paste0("figures/", datasetName))


# Volcano Plot

volcanoDetails <- paste(datasetName, "- Volcano")

volcanoPlotExport(inputFilePath,volcanoDetails,paste0("figures/", datasetName))

# Exporting FGSEA
write_csv(fgsea, fgseaExportPath)