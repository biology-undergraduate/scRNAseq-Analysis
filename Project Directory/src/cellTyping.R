setwd("~/University/Fourth Year/RNAseq Directory")

source("functions/libraryLoader.R")
source("functions/variableFeatureFinder.R")
source("functions/featureVisualiser.R")

libs <- c("dplyr","Seurat","patchwork","ggplot2","tibble")

library_loader(libs)

pbmc_processed <- readRDS("processedData/pbmc_umap.rds")

pbmc_metaOnly <- readRDS("processedData/pbmc_metadata.rds")

markerData <- FindAllMarkers(pbmc_processed, only.pos = TRUE)

# Plotting features to UMAP

f1 <- FeaturePlot(pbmc_processed, features = "CD4:RPA-T4-CD4-AHS0227-pAbO")+
  labs(title = "CD4")
f2 <- FeaturePlot(pbmc_processed, features = "CD8:SK1-CD8A-AHS0228-pAbO")+
  labs(title = "CD8")
f3 <- FeaturePlot(pbmc_processed, features = "CD56:NCAM16.2-NCAM1-AHS0019-pAbO")+
  labs(title = "CD56")
f4 <- FeaturePlot(pbmc_processed, features = "CustomAbSeq11-ACU7011-pAbO")+
  labs(title = "CCR6")
f5 <- FeaturePlot(pbmc_processed, features = "CD45RA:HI100-PTPRC-AHS0009-pAbO")+
  labs(title = "CD45RA")
f6 <- FeaturePlot(pbmc_processed, features = "CD45RO-PTPRC-AHS0036-pAbO")+
  labs(title = "CD45RO")
f7 <- FeaturePlot(pbmc_processed, features = "CustomAbSeq12-ACU7012-pAbO")+
  labs(title = "CD161")
f8 <- FeaturePlot(pbmc_processed, features = "TCR-Valpha7-TRAV7--AHS0282-pAbO")+
  labs(title = "TCR V-Alpha 7.2")
f9 <- FeaturePlot(pbmc_processed, features = "CD69")

CombinePlots(
  plots = list(f1,f2,f3,f4,f5,f6,f7,f8,f9),
  legend = 'right'
)
ggsave("figures/featureUMAP.png", width = 9, height = 9)

# Confirming MAIT location

FeaturePlot(pbmc_processed,features = c('KLRB1', "ZBTB16", "SLC4A10")) 
ggsave("figures/MAITmarkers.png")

# Typing MAITS

pbmc_mait <- subset(pbmc_processed, subset = KLRB1 > 0 & SLC4A10 > 0)

maitMetaData <- pbmc_mait@meta.data

maitMetaData$Cell_Type_Experimental <- "MAIT"

MAIT_barcodes = pbmc_mait@assays[["RNA"]]@counts@Dimnames[[2]]

pbmc_metaOnly[MAIT_barcodes, "Cell_Type_Experimental"] = "MAIT"
pbmc_metaOnly[MAIT_barcodes, "Cell_Type_Experimental"] 

pbmc_typed <- AddMetaData(pbmc_processed, pbmc_metaOnly)

# Visualisation

pbmc_typed_var <- varFeatFinder(pbmc_typed, 2000)
featureView(pbmc_typed_var)

all.genes_typed <- rownames(pbmc_typed_var)
pbmc_typed_scale <- ScaleData(pbmc_typed_var, features = all.genes_typed)

pbmc_typed_pca <- RunPCA(pbmc_typed_scale)

ElbowPlot(pbmc_typed_pca)

DimHeatmap(pbmc_typed_pca, dims = 1:10, balanced = TRUE)
savePlot("figures/dimplotHeatmap.png", type = "png")

pbmc_typed_neighbors <- FindNeighbors(pbmc_typed_pca, dims = 1:10)

pbmc_typed_cluster <- FindClusters(pbmc_typed_neighbors, resolution = 0.5)

head(Idents(pbmc_typed_cluster), n = 5)

pbmc_typed_umap <- RunUMAP(pbmc_typed_cluster, dims = 1:15)

DimPlot(pbmc_typed_umap, group.by = "Cell_Type_Experimental")

KLRB1_data <- FetchData(pbmc_typed_umap, vars = "KLRB1")
SLC4A10_data <- FetchData(pbmc_typed_umap, vars = "SLC4A10")
DimPlot(pbmc_typed_umap, label = TRUE,cols = "grey", cells.highlight = which(KLRB1_data > 0 & SLC4A10_data > 0), cols.highlight = "red") + 
  NoLegend()+
  labs(title = "MAIT cells")

savePlot("figures/MAITtypingUMAP.png", type = "png")

saveRDS(pbmc_typed_umap, file = "processedData/pbmc_typed_umap.rds")
