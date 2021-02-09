source("code/log_smk.R")

doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

data_processed <- readRDS(snakemake@input[["rds"]])$dat_transformed
groups <- tibble::deframe(readr::read_csv(snakemake@input[["groups"]]))

print(snakemake@params[["method"]])
hps <- NULL
if(snakemake@params[["method"]] == 'rf') hps <- list(mtry = c(10, 30, 50, 70, 90, 110))
if(snakemake@params[["method"]] == 'glmnet') hps <- list(alpha = c(0, 0.001, 0.01), lambda = c(0.0001,0.001,0.01, 0.05, 0.1,0.25,0.75,1))
if(snakemake@params[["method"]] == 'svmRadial') hps <- list(sigma = c(0.00001,0.0001,0.001,0.005), C = c(0.01,0.1,1,5,10))

print(hps)

ml_results <- mikropml::run_ml(
  dataset = data_processed,
  method = snakemake@params[["method"]],
  outcome_colname = snakemake@params[['outcome_colname']],
  hyperparameters = hps,
  find_feature_importance = TRUE,
  kfold = as.numeric(snakemake@params[['kfold']]),
  groups = groups,
  seed = as.numeric(snakemake@params[["seed"]])
)

saveRDS(ml_results$trained_model, file = snakemake@output[["model"]])
readr::write_csv(ml_results$performance, snakemake@output[["perf"]])
readr::write_csv(ml_results$feature_importance, snakemake@output[["feat"]])
