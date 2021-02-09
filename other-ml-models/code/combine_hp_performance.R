source("code/log_smk.R")
library(tidyverse)

models <- lapply(snakemake@input[["rds"]], function(x) readRDS(x))

methods <- sapply(snakemake@input[["rds"]], function(filename){
       method = str_replace(filename, "^results/runs/(.*)_(.*)_model.Rds", "\\1")
})



for(method in unique(methods)){
  print(method)
  hp_perf <- mikropml::combine_hp_performance(models[methods == method])
  hp_perf$method <- method
  saveRDS(hp_perf, file = snakemake@output[["rds"]][grep(paste0(method,'.Rds$'), snakemake@output[["rds"]])])
}

