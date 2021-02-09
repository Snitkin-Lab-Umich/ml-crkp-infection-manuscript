source("code/log_smk.R")
library(tidyverse)

#models <- lapply(snakemake@input[["rds"]], function(x) readRDS(x))

#methods <- sapply(snakemake@input[["rds"]], function(filename){
#       method = str_replace(filename, "^results/runs/(.*)_(.*)_model.Rds", "\\1")
#})

for(i in snakemake@input[["rds"]]){
  print(i)
  hp_perf <- readRDS(i)
  method <- unique(hp_perf$method)
  print(method)
  hp_plot_list <- lapply(hp_perf$params, function(param){
  mikropml::plot_hp_performance(hp_perf$dat, !!sym(param), !!sym(hp_perf$metric)) + theme_classic() + scale_color_brewer(palette = "Dark2")
  })
  print(hp_plot_list)
  hp_plot <- cowplot::plot_grid(plotlist = hp_plot_list)
  ggsave(snakemake@output[["plot"]][grep(paste0(method,'.png$'), snakemake@output[["plot"]])], plot = hp_plot)
}

