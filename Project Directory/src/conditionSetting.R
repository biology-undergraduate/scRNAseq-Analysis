setwd("~/University/Fourth Year/RNAseq Directory")

source("functions/libraryLoader.R")

libs <- c("dplyr","Seurat","patchwork","ggplot2","tibble")

library_loader(libs)

pbmc_cellTyped <- readRDS('processedData/pbmc_typed_umap.rds')

# Adding a conditions column

pbmc_cellTyped_meta <- pbmc_cellTyped@meta.data

pbmc_cellTyped_meta[pbmc_cellTyped_meta$Sample_Name %in%
                      c("PK32_IFNa2","PK33_IFNa2","PK34_IFNa2"), "condition"] = "IFN-a2"

pbmc_cellTyped_meta[pbmc_cellTyped_meta$Sample_Name %in%
                      c("PK32_IFNL3","PK33_IFNL3","PK34_IFNL3"), "condition"] = "IFN-l3"

pbmc_cellTyped_meta[pbmc_cellTyped_meta$Sample_Name %in% 
                      c("PK32_IFNg","PK33_IFNg","PK34_IFNg"), "condition"] = "IFN-g"

pbmc_cellTyped_meta[pbmc_cellTyped_meta$Sample_Name %in% 
                      c("PK32_US","PK33_US","PK34_US"), "condition"] = "Unstimulated"

pbmc_condition <- AddMetaData(pbmc_cellTyped,pbmc_cellTyped_meta)

# Cell type subsetting

MAIT <- pbmc_condition[, pbmc_condition@meta.data$Cell_Type_Experimental == "MAIT"]

non_MAIT <- pbmc_condition[, pbmc_condition@meta.data$Cell_Type_Experimental != "MAIT"]

CD8mem <- pbmc_condition[, pbmc_condition@meta.data$Cell_Type_Experimental == "T_CD8_memory"]

# Markers bulk

alphaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) 

gammaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) 

lambdaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
)

# Markers MAIT

MAIT_alphaMarkers <- FindMarkers(
  MAIT,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) 

MAIT_gammaMarkers <- FindMarkers(
  MAIT,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) 

MAIT_lambdaMarkers <- FindMarkers(
  MAIT,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
)

# Markers non-MAIT

nonMAIT_alphaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) 

nonMAIT_gammaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) 

nonMAIT_lambdaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
)

# Markers CD8mem

CD8mem_alphaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) 

CD8mem_gammaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) 

CD8mem_lambdaMarkers <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
)

# Identifying Upregulation bulk

alphaMarkersUP <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
  ) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

gammaMarkersUP <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

lambdaMarkersUP <- FindMarkers(
  pbmc_condition,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

# Identifying Upregulation MAIT

MAIT_alphaMarkersUP <- FindMarkers(
  MAIT,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

MAIT_gammaMarkersUP <- FindMarkers(
  MAIT,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

MAIT_lambdaMarkersUP <- FindMarkers(
  MAIT,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

# Identifying Upregulation non-MAIT

nonMAIT_alphaMarkersUP <- FindMarkers(
  non_MAIT,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

nonMAIT_gammaMarkersUP <- FindMarkers(
  non_MAIT,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

nonMAIT_lambdaMarkersUP <- FindMarkers(
  non_MAIT,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

# Identifiying Upregulation CD8-Memory

CD8mem_alphaMarkersUP <- FindMarkers(
  CD8mem,
  ident.1 = "IFN-a2", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

CD8mem_gammaMarkersUP <- FindMarkers(
  CD8mem,
  ident.1 = "IFN-g", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

CD8mem_lambdaMarkersUP <- FindMarkers(
  CD8mem,
  ident.1 = "IFN-l3", ident.2 = "Unstimulated",
  group.by = "condition"
) %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5)

# Identifying cell numbers

cellNumbers <- table(pbmc_condition$Cell_Type_Experimental)
view(cellNumbers)

cellNumbers_condition <- table(pbmc_condition$condition,pbmc_condition$Cell_Type_Experimental)
View(cellNumbers_condition)

# Visualising condition distribution on umap

DimPlot(pbmc_condition, group.by = "condition",split.by = "Cell_Type_Experimental") +
  labs(title = "Condition Distribution of Cell Types")
savePlot("figures/conditionDistributionUMAP.png", type = "png")

# Export data alpha

write.csv(alphaMarkers, "processedData/alpha_unstim.csv")
write.csv(MAIT_alphaMarkers, "processedData/MAIT_alpha_unstim.csv")
write.csv(nonMAIT_alphaMarkers, "processedData/nonMAIT_alpha_unstim.csv")
write.csv(CD8mem_alphaMarkers, "processedData/CD8mem_alpha_unstim.csv")

write.csv(alphaMarkersUP, "processedData/alpha_unstim_UP.csv")
write.csv(MAIT_alphaMarkersUP, "processedData/MAIT_alpha_unstim_UP.csv")
write.csv(nonMAIT_alphaMarkersUP, "processedData/nonMAIT_alpha_unstim_UP.csv")
write.csv(CD8mem_alphaMarkersUP, "processedData/CD8mem_alpha_unstim_UP.csv" )

# Export data gamma

write.csv(gammaMarkers, "processedData/gamma_unstim.csv")
write.csv(MAIT_gammaMarkers, "processedData/MAIT_gamma_unstim.csv")
write.csv(nonMAIT_gammaMarkers, "processedData/nonMAIT_gamma_unstim.csv")
write.csv(CD8mem_gammaMarkers, "processedData/CD8mem_gamma_unstim.csv")

write.csv(gammaMarkersUP, "processedData/gamma_unstim_UP.csv")
write.csv(MAIT_gammaMarkersUP, "processedData/MAIT_gamma_unstim_UP.csv")
write.csv(nonMAIT_gammaMarkersUP, "processedData/nonMAIT_gamma_unstim_UP.csv")
write.csv(CD8mem_gammaMarkersUP, "processedData/CD8mem_gamma_unstim_UP.csv" )

# Export data Lambda

write.csv(lambdaMarkers, "processedData/lambda_unstim.csv")
write.csv(MAIT_lambdaMarkers, "processedData/MAIT_lambda_unstim.csv")
write.csv(nonMAIT_lambdaMarkers, "processedData/nonMAIT_lambda_unstim.csv")
write.csv(CD8mem_lambdaMarkers, "processedData/CD8mem_lambda_unstim.csv")

write.csv(lambdaMarkersUP, "processedData/lambda_unstim_UP.csv")
write.csv(MAIT_lambdaMarkersUP, "processedData/MAIT_lambda_unstim_UP.csv")
write.csv(nonMAIT_lambdaMarkersUP, "processedData/nonMAIT_lambda_unstim_UP.csv")
write.csv(CD8mem_lambdaMarkersUP, "processedData/CD8mem_lambda_unstim_UP.csv" )

# Export data cell numbers

write.csv(as.data.frame(cellNumbers), "processedData/cellNumbers.csv")
write.csv(as.data.frame(cellNumbers_condition), "processedData/cellNumbersCondition.csv")

# Save RDS

saveRDS(pbmc_condition,"processedData/pbmc_condition.rds")
