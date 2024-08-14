gene.list <- c("DE4", "DEall", "meta_intersect_unionDE" , "5k" , "meta") # "de_intersect"  
title <- c("4 datasets DE, scaled data", "all datasets DE, scaled data",
           "intersect meta genes with union DE 4 and all datasets, scaled data",
           "most varing 5000k genes, scaled data", "meta genes, scaled data") 
  # "intersect results for 4 datasets and all datasets DE, scaled data"

lapply(seq(1:5), function(x) { 
  rmarkdown::render("/data/pangelin/HUG/Thibaud/RC_GEOMX/analysis/scripts/modeling_extra_per_dataset.Rmd", 
                    output_file=paste0("/data/pangelin/HUG/Thibaud/RC_GEOMX/analysis/scripts/modeling_extra_per_dataset_",gene.list[x],"_scaled.html"), 
  params = list(gene.list = gene.list[x], title=title[x])) })
