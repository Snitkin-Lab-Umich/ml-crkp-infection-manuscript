# looking into importance of clinical features from lr results (not permutation version of results, but needed to determine positive association with colonization vs. infection)

library(ggplot2)
library(reshape2)
library(forcats)
library(dplyr)

dirs = c('results/infection_patient-kleborate_lr','results/infection-resp_patient-kleborate_lr','results/infection-urine_patient-kleborate_lr')

dir.create('results/figures/importances', showWarnings=F)
outfiles = paste0('results/figures/importances/',gsub('results/','',dirs),'_importances.pdf')

med_imp_all = data.frame(feature=NA,median=NA,pos_frac=NA,median_raw=NA,dat=NA)
imp_all = data.frame(matrix(NA,nrow = 1, ncol = 8))
colnames(imp_all) = c('feature','imp_rank_raw','imp_rank','imp_sign','dat','seed','imp_rank_frac_raw','imp_rank_frac')

for(j in 1:length(dirs)){
  print(dirs[[j]])
  paths = list.files(dirs[[j]],full.names = T)
  names(paths) = list.files(dirs[[j]],full.names = F)
  
  load(paths[1])
  
  importances = data.frame(matrix(NA,nrow = 1, ncol = 6))
  colnames(importances) = c('feature','imp_rank_raw','imp_rank','imp_sign','dat','seed')
  
  for(i in 1:length(paths)){
    load(paths[i])
    dat = data.frame(feature=colnames(results$feature_importance_cor),
                     imp_rank_raw=rank((results$feature_importance_cor[1,])),
                     imp_rank=rank(abs(results$feature_importance_cor[1,])),
                     imp_sign=ifelse(sign(results$feature_importance_cor[1,])==1,'Positive','Negative'),
                     dat=gsub('.*/|.RData','',paths[i]) %>% sub("_[^_]+$", "", .),
                     seed=gsub('.*_|.RData','',paths[i]))
    importances = rbind(importances,dat)
  }
  
  importances = data.frame(importances[2:nrow(importances),])
  num_features = length(unique(importances$feature))
  importances = importances %>% mutate(imp_rank_frac=imp_rank/num_features)
  importances = importances %>% mutate(imp_rank_frac_raw=imp_rank_raw/num_features)
  imp_all = rbind(imp_all,importances)
  
  imp_summary = as.data.frame(importances %>% group_by(feature) %>% 
                                summarise(median = median(imp_rank_frac, na.rm = TRUE),
                                          pos_frac = mean(imp_sign=='Positive'),
                                          median_raw = median(imp_rank_frac_raw, na.rm = TRUE)))
  imp_summary = imp_summary[order(imp_summary$median),]

  imp_summary$dat = rep(gsub('results/|/.*','',paths[i]),nrow(imp_summary))
  med_imp_all = rbind(med_imp_all,imp_summary)

}

imp_all = data.frame(imp_all[2:nrow(imp_all),])
write.table(imp_all,file='results/importances.tsv',sep='\t')
med_imp_all = data.frame(med_imp_all[2:nrow(med_imp_all),])
write.table(med_imp_all,file='results/median_importances.tsv',sep='\t')






