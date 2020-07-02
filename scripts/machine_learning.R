# do machine learning

source('scripts/get_model_matrix.R')

# The dependinces for this script are consolidated in the first part
deps = c("dplyr", "tictoc", "caret" ,"rpart", "xgboost", "randomForest", "kernlab","LiblineaR", "pROC", "tidyverse","PRROC");
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org", dependencies=TRUE);
  }
  suppressPackageStartupMessages(library(dep, verbose=FALSE, character.only=TRUE))
}

# source functions
files = paste0('ML_pipeline_microbiome/code/R/',c('model_pipeline.R','tuning_grid.R','permutation_importance.R'))
s = sapply(files,source)
args = commandArgs(trailingOnly = T)
fname = args[[3]]
pref = gsub('.*/|.RData','',fname)
dir = gsub('_[^_]+$','',pref)
out_feat = gsub('_[^_]+$','',dir)
model = gsub('.*_','',dir)

models = c('L2_Logistic_Regression','Random_Forest','XGBoost')
names(models) = c('lr','rf','xgb')
model = models[names(models) == model]

# load data
dat = get_model_matrix(args[[1]],incl_out=T)
dat <- data.frame(dat[order(dat[,1],decreasing = T), ],stringsAsFactors=F) # make 1st outcome value 'yes'
# Remove features with near zero variance and scale remaining from 0 to 1
preProcValues <- preProcess(dat, method = c("nzv", "range"))
dataTransformed <- predict(preProcValues, dat)
dataTransformed[,1] = as.factor(ifelse(dataTransformed[,1]==1,'yes','no'))

# get groups (LTACHs)
ltachs = read.delim('data/ltachs.tsv',stringsAsFactors=F)
grp = ltachs[,1]
names(grp) = rownames(ltachs)
grp = grp[rownames(dat)]
# set seed
seed = as.numeric(gsub('.*_','',pref))
set.seed(seed)
# no permutation
perm = F
# hyperparameters
hp = args[[2]]

dir.create('data/temp/NA', showWarnings=F)
results = model_pipeline(dataTransformed, model, seed, permutation=perm, hyperparameters=hp,group=grp)

# output of pipeline is list:
# 1. cv_auc
# 2. test_auc
# 3. results_individual
# 4. feature_importance_non_cor
# 5. feature_importance_cor
# 6. trained_model

save(results,file = fname)


