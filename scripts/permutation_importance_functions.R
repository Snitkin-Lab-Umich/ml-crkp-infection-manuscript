# permutation importance functions
library(tidyverse)
library(caret) 
library(Hmisc)

source('scripts/get_model_matrix.R')

# Usage: input file - "data/input_data.csv"
#        outcome - e.g. "dx"
#        level - name of modeling experiment
#        cor_value - select correlations greater than or equal to cor_value
#        p_value - select correlation with value below p_value
compute_correlation_matrix <- function(input_file, outcome, output_file, cor_value = 1, p_value = 0.01){
  
  ############### READ IN THE INPUT DATA ###############
  #data_corr <- read_tsv(input_file)
  # load data
  dat = get_model_matrix(input_file,incl_out=T)
  dat <- data.frame(dat[order(dat[,1],decreasing = T), ],stringsAsFactors=F) # make 1st outcome value 'yes'
  # Remove features with near zero variance and scale remaining from 0 to 1
  preProcValues <- preProcess(dat, method = c("nzv", "range"))
  data_corr <- predict(preProcValues, dat)
  data_corr[,1] = as.factor(ifelse(data_corr[,1]==1,'yes','no'))
  # remove outcome, only keep the features
  data_corr <- data_corr[,!grepl(outcome, names(data_corr))]
  #######################################################
  
  
  ########### COMPUTE CORRELATION MATRIX ##################
  r <- rcorr(as.matrix(data_corr), type="spearman")
  
  adjusted <- p.adjust(r$P, method = "holm")
  r$P <- adjusted
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
  new_r <- flattenCorrMatrix(r$r, r$P) %>%
    filter(cor>=cor_value) %>%
    filter(p<p_value) %>%
    write_tsv(output_file)
  ##########################################################
}

