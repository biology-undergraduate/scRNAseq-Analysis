generate_and_export_plots <- function(conditionDeframed, geneSetList, plotDetails, output_dir = "plots") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Plotting GSEA of interest
  plot_gsea <- function(gene_set, title, output_file) {
    plot <- plotEnrichment(gene_set, conditionDeframed) +
      labs(title = title, x = "Rank", y = "Enrichment Score") +
      geom_line(aes(x = rank, y = ES), color = "green")
    return(plot)
  }
  
  # Generate and export plots
  for (plotDetail in plotDetails) {
    gene_set <- geneSetList[[plotDetail$gene_set_name]]
    title <- plotDetail$title
    # Automate file name generation based on title
    output_file <- paste0(output_dir, "/", gsub(" ", "_", title), ".png")
    
    plot <- plot_gsea(gene_set, title, output_file)
    # Use ggsave to save the plot
    ggsave(filename = output_file, plot = plot, width = 8, height = 8, dpi = 300)
  }
}
