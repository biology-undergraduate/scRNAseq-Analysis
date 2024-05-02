qcPlotter <- function(data) {
  
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  plot1 + plot2
  
}