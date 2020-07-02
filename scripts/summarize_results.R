# Summarize ml results

# output of pipeline is list:
# 1. cv_auc
# 2. test_auc
# 3. results_individual
# 4. feature_importance_non_cor
# 5. feature_importance_cor
# 6. trained_model

library(ggplot2)
library(tidyr)

dirs = list.dirs('results')
figdirs = list.dirs('results/figures')
permdirs = list.dirs('results/permutation_importance')
dirs = dirs[!dirs %in% c('results',figdirs,permdirs)]

ltachs = read.delim('data/ltachs.tsv')
names_ltachs = rownames(ltachs)
ltachs = ltachs$metadata.ltach
names(ltachs) = names_ltachs

auc_cols = c('outcome','feature','model','mean_train_auc','sd_train_auc','mean_test_auc','sd_test_auc')
aucs = data.frame(matrix(0,nrow=length(dirs),ncol=length(auc_cols)))
colnames(aucs) = auc_cols
rownames(aucs) = dirs

aucs_all = data.frame(seed=NA,dat=NA,auc=NA,cv_auc=NA,train_ltachs=NA,test_ltachs=NA,auprc=NA)

for(i in 1:length(dirs)){
  d = dirs[i]
  print(d)
  pref = gsub('results/','',d)
  fsplit = strsplit(pref,'_')[[1]]
  aucs$outcome[i] = fsplit[1]
  aucs$feature[i] = fsplit[2]
  aucs$model[i] = fsplit[3]

  files = list.files(d,full.names = T)
  if(length(files) == 0) next()
  train_aucs = rep(NA,length(files))
  test_aucs = rep(NA,length(files))
  auprcs = rep(NA,length(files))
  train_ltachs = rep(NA,length(files))
  test_ltachs = rep(NA,length(files))
  seeds = rep(NA,length(files))
  for(j in 1:length(files)){
    f = files[j]
    seed = gsub('.*_|.RData','',f)
    seeds[j] = seed
    load(f)
    train_aucs[j] = results$cv_auc
    test_aucs[j] = results$test_auc
    auprcs[j] = results$auprc
    train_ltachs[j] = paste0(unique(ltachs[rownames(results$trainData)]),collapse = '')
    test_ltachs[j] = paste0(unique(ltachs[rownames(results$testData)]),collapse = '')
    rm(results)
  }
  
  aucs$mean_train_auc[i] = mean(train_aucs)
  aucs$median_train_auc[i] = median(train_aucs)
  aucs$sd_train_auc[i] = sd(train_aucs)
  aucs$mean_test_auc[i] = mean(test_aucs)
  aucs$median_test_auc[i] = median(test_aucs)
  aucs$sd_test_auc[i] = sd(test_aucs)
  
  
  dat = data.frame(seed=seeds,dat=rep(pref,length(test_aucs)),
                   auc=test_aucs,cv_auc=train_aucs,
                   test_ltachs=test_ltachs,train_ltachs=train_ltachs,auprc=auprcs)
  aucs_all = rbind(aucs_all,dat)
  
}

rownames(aucs) = 1:nrow(aucs)
aucs_all = data.frame(aucs_all[2:nrow(aucs_all),])
write.table(aucs_all,'results/test_aucs.tsv')
aucs_all$id <- rep(1:sum(aucs_all$dat==aucs_all$dat[1]),times = length(unique(aucs_all$dat)))
head(aucs_all)
aucs_wide = spread(aucs_all[,c('dat','auc','id')],dat,auc)
head(aucs_wide)
quantiles = apply(aucs_wide[,2:ncol(aucs_wide)],2,function(x){
  c(mean=mean(x),
    median=median(x),
    quantile(x,0.25),
    quantile(x,0.75))
})
write.table(quantiles,'results/quantiles.tsv')


pdf('results/figures/test_aucs_boxplot.pdf')
ggplot(aucs_all,aes(x=reorder(dat,auc,median),y=auc)) + geom_boxplot() + theme_bw() +  theme(axis.text.x=element_text(angle=90,hjust=1)) + xlab('')
dev.off()
aucs_all_subset = aucs_all[aucs_all$dat %in% c('infection_patient_lr','infection_genomic_lr','infection_patient-genomic_lr'),]
aucs_all_subset$dat[aucs_all_subset$dat == 'infection_patient_lr'] = 'Clinical'
aucs_all_subset$dat[aucs_all_subset$dat == 'infection_genomic_lr'] = 'Genomic'
aucs_all_subset$dat[aucs_all_subset$dat == 'infection_patient-genomic_lr'] = 'Clinical & Genomic'
aucs_all_subset$dat = factor(aucs_all_subset$dat, levels = c('Clinical','Genomic','Clinical & Genomic'), ordered=T)

write.table(aucs,file = snakemake@output[[1]],quote = F,sep = '\t',row.names = T,col.names = T)
