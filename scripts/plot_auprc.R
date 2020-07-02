# plot aucs (Figure S3)
library(tidyverse)
library(cowplot)
library(grid)

source('scripts/geom_boxplot2_base_r_look.R')
source('scripts/build_polygon.R')

options(OutDec="·")

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
  theme(text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1)) + 
  xlab('') + ylab('AUPRC') +
  geom_hline(yintercept=null_dat, lty=3) + 
  scale_x_discrete(labels=rev(c('Uncurated\ngenomic','Uncurated\ngrouped genomic','Curated\ngenomic','Patient','Patient &\nCurated genomic'))) +
  scale_fill_manual(breaks=rev(c('#e5f5f9','#99d8c9','#2ca25f','#fa9fb5','#756bb1')),
                    values=rev(c('#e5f5f9','#99d8c9','#2ca25f','#fa9fb5','#756bb1')),
                    labels=rev(c('Uncurated genomic','Uncurated grouped genomic','Curated genomic','Patient','Patient & Curated Genomic'))) +
  coord_flip() + scale_y_continuous(expand = c(0.01, 0.01))

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

p2 = ggplot(data = source_aucs, aes(x = feature, y = auprc, fill = source)) +
  geom_boxplot2() + facet_grid(cols=vars(source),
                               labeller = labeller(source=c('resp'='Respiratory','urine'='Urinary'))) +
  geom_hline(data=dummy2,aes(yintercept=Z), lty=3) + 
  theme(text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1))+
  xlab('') + ylab('AUPRC') + 
  scale_x_discrete(labels=c('Genomic','Patient','Patient & Genomic'))  +
  theme(legend.position='none') +
  scale_fill_manual(breaks=c(resp='#2b8cbe',urine='#feb24c'),
                    values=c(resp='#2b8cbe',urine='#feb24c'))

source_auc_summary = aggregate(auprc~dat,source_aucs,summary)
source_auc_summary

pdf('results/figures/aucs/aurocs.pdf',width = 10,height = 5)
plot_grid(p1,p2,labels = 'AUTO',rel_widths = c(0.6,0.4),axis = 'b',align = 'hv')
dev.off()

