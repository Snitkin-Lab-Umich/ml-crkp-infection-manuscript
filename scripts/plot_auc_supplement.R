# auc mega plot

# plot aucs
library(tidyverse)
library(cowplot)
library(grid)
library(colorspace)

source('scripts/geom_boxplot2_base_r_look.R')
source('scripts/build_polygon.R')

#options(OutDec="??")

# get baseline for auprc
dat = read.delim('data/combined/infection_patient-kleborate.tsv')
rdat = read.delim('data/combined/infection-resp_patient-kleborate.tsv')
udat = read.delim('data/combined/infection-urine_patient-kleborate.tsv')

null_dat = mean(dat['infection',] == 'yes')
null_rdat = mean(rdat['infection',] == 'yes')
null_udat = mean(udat['infection',] == 'yes')

aucs = read.delim('results/test_aucs.tsv',sep=' ')


inf_aucs = aucs[aucs$dat %in% c('infection_genomic_lr','infection_gene_lr','infection_kleborate_lr','infection_patient_lr','infection_patient-kleborate_lr'),]
inf_aucs = inf_aucs %>% mutate(dat = factor(dat, levels=rev(c('infection_genomic_lr','infection_gene_lr','infection_kleborate_lr','infection_patient_lr','infection_patient-kleborate_lr'))))
inf_aucs = inf_aucs[order(inf_aucs$dat),]

auc_summary = aggregate(auprc~dat,inf_aucs,summary)
auc_summary
round(auc_summary$auprc[,5] - auc_summary$auprc[,2],2)

# plot things
size = 10

p1 = ggplot(data = inf_aucs, aes(x = dat, y = auprc, fill = dat)) +
  geom_boxplot2() + 
  theme(text = element_text(size=size), legend.position = 'None') + 
  xlab('Feature set') + ylab('AUPRC') + #ylim(0.4,0.8) +
  geom_hline(yintercept=null_dat, lty=3) + 
  scale_x_discrete(labels=rev(c('Uncurated\ngenomic','Uncurated\ngrouped genomic','Curated\ngenomic','Patient','Patient &\nCurated genomic'))) +
  scale_fill_manual(#breaks=rev(c('#e5f5f9','#99d8c9','#2ca25f','#fa9fb5','#756bb1')),
    values=rev(c('#e5f5f9','#99d8c9','#2ca25f','#fa9fb5','#756bb1')),#c('#e5f5f9','#99d8c9','darkseagreen4','bisque4','#756bb1')),
    labels=rev(c('Uncurated genomic','Uncurated grouped genomic','Curated genomic','Patient','Patient & Curated Genomic'))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1))

source_names = c('infection-urine_kleborate_lr',
                 'infection-urine_patient_lr',
                 'infection-urine_patient-kleborate_lr',
                 'infection-resp_kleborate_lr',
                 'infection-resp_patient_lr',
                 'infection-resp_patient-kleborate_lr')
source_aucs = aucs[aucs$dat %in% source_names,]
source_aucs = source_aucs %>% mutate(dat = factor(dat, levels=source_names)) 
source_aucs$source = gsub('infection-|_.*','',source_aucs$dat)
source_aucs$feature = gsub('infection-|_lr|resp_|urine_','',source_aucs$dat)

dummy2 <- data.frame(source = c("resp", "urine"), Z = c(null_rdat,null_udat))

p4 = ggplot(data = source_aucs, aes(x = feature, y = auprc, fill = source)) +
  geom_boxplot2() + facet_grid(cols=vars(source),
                               labeller = labeller(source=c('resp'='Respiratory','urine'='Urinary'))) +
  geom_hline(data=dummy2,aes(yintercept=Z), lty=3) + 
  theme(text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1))+#,
  #axis.text.y=element_blank(),
  #axis.ticks.y=element_blank()) + 
  xlab('Feature set') + ylab('AUPRC') + #ylim(0.35,0.85) + 
  scale_x_discrete(labels=c('Genomic','Patient','Patient & Genomic'))  +
  theme(legend.position='none') +
  scale_fill_manual(breaks=c(resp='#2b8cbe',urine='#feb24c'),
                    values=c(resp='#2b8cbe',urine='#feb24c'))


p3 = ggplot(data = source_aucs, aes(x = feature, y = auc, fill = source)) +
  geom_boxplot2() + facet_grid(cols=vars(source),
                               labeller = labeller(source=c('resp'='Respiratory','urine'='Urinary'))) +
  geom_hline(yintercept=0.5, lty=3) + 
  theme(text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1))+#,
  #axis.text.y=element_blank(),
  #axis.ticks.y=element_blank()) + 
  xlab('Feature set') + ylab('AUROC') + ylim(0.35,0.85) + 
  scale_x_discrete(labels=c('Genomic','Patient','Patient & Genomic'))  +
  theme(legend.position='none') +
  scale_fill_manual(breaks=c(resp='#2b8cbe',urine='#feb24c'),
                    values=c(resp='#2b8cbe',urine='#feb24c'))#'khaki3','cornflowerblue'),
#labels=c('resp','urine'))

l2 <- read.delim('results/test_aucs.tsv',sep = ' ', ) %>% filter(dat == 'infection_patient-kleborate_lr') %>%
  mutate(method = 'l2') #%>% 
names(l2)[names(l2) == 'auc'] <- 'AUC' 
l2 <- l2 %>% select(AUC, method, seed)
others <- read_csv('../2021-01-05_ml/results/performance_results.csv') %>% select(AUC, method, seed)
perf <- bind_rows(l2, others)

p2 <- perf %>% 
  mutate(method=replace(method, method=='glmnet', 'Elastic net'),
         method=replace(method, method=='l2', 'L2 logistic\nregression'),
         method=replace(method, method=='rf', 'Random forest'),
         method=replace(method, method=='svmRadial', 'SVM'),
         method=factor(method, levels = c('L2 logistic\nregression', 'Elastic net', 'SVM', 'Random forest'))) %>%
  ggplot(aes(x = method, y = AUC, fill = method)) + geom_boxplot2() +
  labs(x = 'Machine learning method', y = 'AUROC', fill = 'Method') + 
  scale_fill_manual(values = c('#756bb1',lighten('#756bb1',0.3),lighten('#756bb1',0.6),lighten('#756bb1',0.9))) + 
  geom_hline(yintercept=0.5, lty=3) + theme(text = element_text(size=15)) +
  theme(legend.position='none',text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1))

ggsave('results/figures/aucs/auc_supplement.tiff', plot_grid(p1, p2, p3, p4, nrow = 2, labels = 'AUTO'), width = 5, height = 6)
