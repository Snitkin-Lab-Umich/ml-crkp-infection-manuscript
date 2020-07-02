# Code written by Zena Lapp
#' Calculate auroc and auprc
calc_aucs <- function(pred,outcome){
  # get the one with fewer samples to calculate auprc
  fewer_samples <- which.min(colSums(pred))
  bin_outcome <- ifelse(outcome == names(pred)[fewer_samples], 1, 0)
  auroc <- roc.curve(pred[[fewer_samples]], weights.class0 = bin_outcome)$auc
  auprc <- pr.curve(pred[[fewer_samples]], weights.class0 = bin_outcome)$auc.integral
  return(list(auroc=auroc,auprc=auprc))
}
