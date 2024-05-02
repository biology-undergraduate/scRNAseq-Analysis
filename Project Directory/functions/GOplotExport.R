generate_and_export_GO <- function(file, details, output_dir = "plots") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Plotting GO
  plot_GO <- function(data, title) {
    GO_data <- read_csv(data)
    GO_processed <- as.character(GO_data$...1)
    
    eGO <- enrichGO(gene = GO_processed,
                    OrgDb         = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.001,
                    qvalueCutoff = 0.001,
                    readable      = FALSE)
    print(eGO)
    
    plot <- dotplot(eGO, showCategory = 20) +
      ggtitle(title) +
      theme(text = element_text(size = 9))
    
    return(plot)
  }
  
  # Generate and export plot
  output_file <- paste0(output_dir, "/", gsub(" ", "_", details), ".png")
  
  plot <- plot_GO(file, details)
  # Use ggsave to save the plot
  ggsave(filename = output_file, plot = plot, width = 8, height = 15, dpi = 300)
}

