# plot aucs
library(tidyverse)
library(cowplot)
library(grid)

source('scripts/geom_boxplot2_base_r_look.R')
source('scripts/build_polygon.R')

#options(OutDec="·")

aucs = read.delim('results/test_aucs.tsv',sep=' ')

inf_aucs = aucs[aucs$dat %in% c('infection_genomic_lr','infection_gene_lr','infection_kleborate_lr','infection_patient_lr','infection_patient-kleborate_lr'),]
inf_aucs = inf_aucs %>% mutate(dat = factor(dat, levels=rev(c('infection_genomic_lr','infection_gene_lr','infection_kleborate_lr','infection_patient_lr','infection_patient-kleborate_lr'))))
inf_aucs = inf_aucs[order(inf_aucs$dat),]

auc_summary = aggregate(auc~dat,inf_aucs,summary)
auc_summary
round(auc_summary$auc[,5] - auc_summary$auc[,2],2)

# plot things
size = 10
axis_title_size = 9

p1 = ggplot(data = inf_aucs, aes(x = dat, y = auc, fill = dat)) +
  geom_boxplot2() + 
  theme(text = element_text(size=size), legend.position = 'None',
        axis.title.x = element_text(size=axis_title_size)) + #, axis.text.x=element_text(angle=45,hjust=1)) + 
  xlab('') + ylab('AUROC') + ylim(0.4,0.8) +
  geom_hline(yintercept=0.5, lty=3) + 
  scale_x_discrete(labels=rev(c('Uncurated\ngenomic','Uncurated\ngrouped genomic','Curated\ngenomic','Patient','Patient &\nCurated genomic'))) +
  scale_fill_manual(#breaks=rev(c('#e5f5f9','#99d8c9','#2ca25f','#fa9fb5','#756bb1')),
                    values=rev(c('#e5f5f9','#99d8c9','#2ca25f','#fa9fb5','#756bb1')),#c('#e5f5f9','#99d8c9','darkseagreen4','bisque4','#756bb1')),
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
p2 = ggplot(data = source_aucs, aes(x = feature, y = auc, fill = source)) +
  geom_boxplot2() + facet_grid(cols=vars(source),
                               labeller = labeller(source=c('resp'='Respiratory','urine'='Urinary'))) +
  geom_hline(yintercept=0.5, lty=3) + 
  theme(text = element_text(size=size), axis.text.x=element_text(angle=45,hjust=1))+
  xlab('') + ylab('AUROC') + ylim(0.35,0.85) + 
  scale_x_discrete(labels=c('Genomic','Patient','Patient & Genomic'))  +
  theme(legend.position='none') +
  scale_fill_manual(breaks=c(resp='#2b8cbe',urine='#feb24c'),
                   values=c(resp='#2b8cbe',urine='#feb24c'))

source_auc_summary = aggregate(auc~dat,source_aucs,summary)
source_auc_summary


compare_aucs = function(fset1,fset2){
  pk_auc = aucs[aucs$dat %in% c(fset1,fset2),]
  dif = data.frame(auc1=pk_auc$auc[pk_auc$dat == fset1],
                   auc2=pk_auc$auc[pk_auc$dat == fset2],
                   auc_diff=pk_auc$auc[pk_auc$dat == fset1] - pk_auc$auc[pk_auc$dat == fset2])
  # two-sided pval
  pval = 2*(min(mean(dif$auc_diff >= 0), mean(dif$auc_diff <= 0)))
  #print(pval)
  return(list(dif=dif,pval=pval))
}

sets = c('infection_patient-kleborate_lr','infection_patient_lr','infection_kleborate_lr','infection_genomic_lr','infection_gene_lr')

# for(i in 1:length(sets)){
#   for(j in 1:length(sets)){
#     if(i < j){
#       namei = gsub('infection_|_lr','',sets[i])
#       namej = gsub('infection_|_lr','',sets[j])
#       comp = compare_aucs(sets[i],sets[j])
#       dif = comp$dif
#       pval = comp$pval
#       p = ggplot(dif, aes(x=auc1,y=auc2)) + geom_point() + 
#         geom_abline(intercept=0) +
#         geom_hline(yintercept=0.5, lty=3) + 
#         geom_vline(xintercept=0.5, lty=3) + 
#         xlab(paste(tools::toTitleCase(namei),'AUROC')) + 
#         ylab(paste(tools::toTitleCase(namej),'AUROC')) +
#         annotation_custom(grobTree(textGrob(paste('p =',pval), x=0.65,  y=0.95, hjust=0,
#                                             gp=gpar(fontsize=12)))) +
#         theme_bw()
#       fname = paste0(namei,'_',namej,'_auc.pdf')
#       ggsave(paste0('results/figures/aucs/',fname),p,width=4,height=4)
#     }
#   }
# }


kg = compare_aucs('infection_kleborate_lr','infection_genomic_lr')
# df_poly_top <- buildPoly(range(kg$dif$auc1),range(kg$dif$auc2),
#                      slope=1,intercept=0,above=TRUE)
# df_poly_bottom <- buildPoly(range(kg$dif$auc1),range(kg$dif$auc2),
#                          slope=1,intercept=0,above=FALSE)
df_poly_top <- buildPoly(c(0.406,0.695),c(0,0.7),
                         slope=1,intercept=0,above=TRUE)
df_poly_bottom <- buildPoly(c(0.406,0.695),c(0,0.7),
                            slope=1,intercept=0,above=FALSE)
kg_auc = ggplot(kg$dif, aes(x=auc1,y=auc2)) + geom_point(alpha=0.7) + 
  geom_abline(intercept=0) +
  geom_hline(yintercept=0.5, lty=3) + 
  geom_vline(xintercept=0.5, lty=3) + 
  ylim(c(0.40,0.8)) + 
  xlim(c(0.42,0.682)) + 
  xlab('Curated genomic AUROC') + 
  ylab('Uncurated genomic\nAUROC') +
  annotate('text',x=0.45, y=0.42,label='y=x') +
  geom_polygon(aes(x=x, y=y), data=df_poly_top, fill='#e5f5f9',alpha=0.2) +
  geom_polygon(aes(x=x, y=y), data=df_poly_bottom, fill='#2ca25f',alpha=0.2) +
  theme_bw() + theme(text = element_text(size=size),
                     axis.title.x = element_text(size=axis_title_size),
                     axis.title.y = element_text(size=axis_title_size))

pk = compare_aucs('infection_kleborate_lr','infection_patient_lr')
# df_poly_top <- buildPoly(range(pk$dif$auc1),range(pk$dif$auc2),
#                          slope=1,intercept=0,above=TRUE)
# df_poly_bottom <- buildPoly(range(pk$dif$auc1),range(pk$dif$auc2),
#                             slope=1,intercept=0,above=FALSE)
pk_auc = ggplot(pk$dif, aes(x=auc1,y=auc2)) + geom_point(alpha=0.7) + 
  geom_abline(intercept=0) +
  geom_hline(yintercept=0.5, lty=3) + 
  geom_vline(xintercept=0.5, lty=3) + 
  ylim(c(0.40,0.8)) + 
  xlim(c(0.42,0.682)) +
  xlab('Curated genomic AUROC') + 
  ylab('Patient AUROC') +
  annotate('text',x=0.45, y=0.42,label='y=x') +
  geom_polygon(aes(x=x, y=y), data=df_poly_top, fill='#fa9fb5',alpha=0.2) +
  geom_polygon(aes(x=x, y=y), data=df_poly_bottom, fill='#2ca25f',alpha=0.2) +
  theme_bw() + theme(text = element_text(size=size),
                     axis.title.x = element_text(size=axis_title_size),
                     axis.title.y = element_text(size=axis_title_size))


# pdf('results/figures/aucs/aucs.tiff',width=10,height=3)
#plot_grid(p1,auc_hist, nrow=1, align = 'h', rel_widths = c(0.65,0.35))
# plot_grid(plot_grid(NULL,p1,NULL,nrow=1,rel_widths = c(0,0.92,0.08),labels=c('','A','')),
#           plot_grid(kg_auc,pk_auc,labels=c('B','C')), 
#           ncol=1)
ggsave('results/figures/aucs/aucs.pdf',plot_grid(p1,kg_auc,pk_auc,nrow=1,labels='AUTO',rel_widths = c(0.4,0.3,0.3)),width=6.87,height=3/(10/6.87), units = 'in', dpi = 300)

ggsave('results/figures/aucs/aucs-test.pdf',plot_grid(p1,kg_auc,pk_auc,nrow=1,labels='AUTO',rel_widths = c(0.4,0.3,0.3)),width=6.87,height=3/(10/6.87), units = 'in', dpi = 72)

ggsave('results/figures/aucs/aucs.tiff',plot_grid(p1,kg_auc,pk_auc,nrow=1,labels='AUTO',rel_widths = c(0.4,0.3,0.3)), device = 'tiff',width=6.87,height=3/(10/6.87), units = 'in', dpi = 300)#,width=10,height=3)
# dev.off()

# ggsave('results/figures/aucs/source_aucs.pdf',p2,width=4,height=5)


