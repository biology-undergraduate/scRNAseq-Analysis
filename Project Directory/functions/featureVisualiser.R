featureView <- function(data) {
  
  top10 <- head(VariableFeatures(data), 10)
  
  plot1 <- VariableFeaturePlot(data)
  
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  plot1 + plot2
  
}