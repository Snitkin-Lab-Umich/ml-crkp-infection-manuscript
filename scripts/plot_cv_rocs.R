# Plot cross-validation AUROCs to be able to tune the models

# load library
library(ggplot2)
library(reshape2)


# output of pipeline is list:
# 1. cv_auc
# 2. test_auc
# 3. results_individual
# 4. feature_importance_non_cor
# 5. feature_importance_cor
# 6. trained_model

dirs = list.dirs('results')
figdirs = list.dirs('results/figures')
permdirs = list.dirs('results/permutation_importance')
dirs = dirs[!dirs %in% c('results',figdirs,permdirs)]

dir.create('results/figures/cv_rocs', showWarnings=F)
outfiles = paste0('results/figures/cv_rocs/',gsub('results/','',dirs),'_cv_rocs.pdf')

roc_summary_all = data.frame(cost=NA,mean_cv_rocs=NA,lower=NA,upper=NA,dat=NA)

for(i in 1:length(dirs)){
  d = dirs[i]
  print(d)
  pref = gsub('results/','',d)
  fsplit = strsplit(pref,'_')[[1]]
  
  outcome = fsplit[1]
  feature = fsplit[2]
  model = fsplit[3]
  
  files = list.files(d,full.names = T)
  test_aucs = rep(NA,length(files))
  names(test_aucs) = 0:(length(files)-1)
  rocs = data.frame(matrix(NA,nrow=length(files),ncol=50))
  costs = c()
  for(j in 1:length(files)){
    #if(j>5) break()
    f = files[j]
    load(f)
    test_aucs[j] = results$cv_auc
    costs = results$trained_model$results$cost
    rocs = rocs[,1:length(costs)]
    rocs[j,] = results$trained_model$results$ROC
    rm(results)
  }
  
  roc_summary = data.frame(cost = costs,
                           mean_cv_rocs = apply(rocs,2,mean,na.rm=T),
                           lower = apply(rocs,2,mean,na.rm=T) - apply(rocs,2,sd,na.rm=T),
                           upper = apply(rocs,2,mean,na.rm=T) + apply(rocs,2,sd,na.rm=T))
  g = ggplot(roc_summary,aes(x=cost,y=mean_cv_rocs))+geom_point()+geom_line()+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3)+theme_bw()+ggtitle(pref)+xlab('Cost') + ylab('Mean CV ROCs')
  ggsave(outfiles[[i]],g, height = 10 , width = 10)
  roc_summary$dat = rep(gsub('results/|/.*','',d),nrow(roc_summary))
  roc_summary_all = rbind(roc_summary_all,roc_summary)
}

roc_summary_all = data.frame(roc_summary_all[2:nrow(roc_summary_all),])
write.table(roc_summary_all,file='results/cv_roc_summary.tsv',sep='\t')

