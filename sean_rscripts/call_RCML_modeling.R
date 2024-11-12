gene_lists <- readRDS('~/RCproject/gene_lists.rds')

print('test')

lapply(names(gene_lists), function(gene_list) {
  output_file <- paste0("~/reports/CV3.modeling_extra_final.", gene_list, ".report.html")
  
  # Print the output file name for tracking
  cat("Rendering file:", output_file, "\n")
  
  rmarkdown::render(
    "~/RCproject_code/sean/RCML_modeling.Rmd", 
    output_file = output_file, 
    params = list(gene.list = gene_list, title = paste0("public signature model CV, gene list ", gene_list))
  )
})  # End of the render reports code chunk