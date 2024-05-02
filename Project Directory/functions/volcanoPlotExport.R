volcanoPlotExport <- function(file, details, output_dir = "plots") {
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Plotting Volcano
  plot_volcano <- function(data, title) {
    
    volcanoData <- read_csv(data)
    
    volcano <- ggplot(volcanoData, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = ifelse(p_val_adj < 0.05 & avg_log2FC > 0.5, "Upregulated",
                                    ifelse(p_val_adj < 0.05 & avg_log2FC < -0.5, "Downregulated", "Not Significant"))),
                 size = 1) +
      geom_text(data = subset(volcanoData, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
                aes(label = subset(volcanoData, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)$...1),
                vjust = -1, hjust = 0, size = 3) +
      labs(x = expression("Log"["2"] ~ "fold" ~ "change"),
           y = expression("-Log"["10"] ~ "adjusted" ~ "p-value"),
           title = title) +
      geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      guides(color = guide_legend(title = ""))
    
    return(volcano)
    
  }
  
  # Generate and export plot
  output_file <- paste0(output_dir, "/", gsub(" ", "_", details), ".png")
  
  plot <- plot_volcano(file, details)
  # Use ggsave to save the plot
  ggsave(filename = output_file, plot = plot, width = 10, height = 7, dpi = 300)
  
}
