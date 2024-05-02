setwd("~/University/Fourth Year/RNAseq Directory")

source("functions/libraryLoader.R")
source("functions/PBMC_cleaner.R")
source("functions/metaGrouping.R")
source("functions/featureVisualiser.R")
source("functions/variableFeatureFinder.R")
source("functions/qcplotter.R")

libs <- c("dplyr","Seurat","patchwork")

library_loader(libs)

pbmc<- readRDS('rawData/IFN-Mar-1_Seurat.rds')

table(pbmc$Sample_Name)

## QC 

# [[]] allows for addition of data to object metadata

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

qcPlotter(pbmc)

# Try to put this into a function "dualPlotScatter"

pbmc_cleaned <- pbmc_QC(pbmc,160,12500,c("Multiplet","Undetermined"))

## Normalization

pbmc_norm <- NormalizeData(pbmc_cleaned)

# Check that standard normalisation assumptions fits with our data (Do they originally have the same number of RNA molecules)

## Feature selection

pbmc_var_feat <- varFeatFinder(pbmc_norm, 2000)

featureView(pbmc_var_feat)

## Data Scaling

all.genes <- rownames(pbmc_var_feat)
pbmc_scale <- ScaleData(pbmc_var_feat, features = all.genes)

# could be changed into a single function

## Linear dimensional reduction - PCA

pbmc_PCA <- RunPCA(pbmc_scale, features = VariableFeatures(object = pbmc_scale))

print(pbmc_PCA[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc_PCA, dims = 1:15, reduction = "pca")

DimPlot(pbmc_PCA, reduction = "pca") + NoLegend()

DimHeatmap(pbmc_PCA, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc_PCA, dims = 1:15, cells = 500, balanced = TRUE)

## Determining Dimensionality

# Discuss whether a jackstraw procedure is a better or worse determinant and see if my desktop is up to the computational load

ElbowPlot(pbmc_PCA)

## Cluster the cells

# The references may be useful to look back on for the write up later

pbmc_neighbors <- FindNeighbors(pbmc_PCA, dims = 1:10)

pbmc_cluster <- FindClusters(pbmc_neighbors, resolution = 0.5)

head(Idents(pbmc_cluster), n = 5)

## Non-linear Dimensional reduction (UMPA/tSNE)

pbmc_UMAP <- RunUMAP(pbmc_cluster, dims = 1:15)

DimPlot(pbmc_UMAP, reduction = "umap")

DimPlot(pbmc_UMAP, reduction = "umap", group.by = "Cell_Type_Experimental")

DimPlot(pbmc_UMAP, reduction = "umap", group.by = "Sample_Name")
ggsave("figures/sampleUMAP.png")

saveRDS(pbmc_UMAP, file = "processedData/pbmc_umap.rds")

metadata <- pbmc_UMAP@meta.data

metaGrouping(metadata)

saveRDS(metadata, file = "processedData/pbmc_metadata.rds")
