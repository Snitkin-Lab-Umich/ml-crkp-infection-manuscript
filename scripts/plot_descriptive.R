# descriptive plots

# load packages
library(readxl)
library(tidyverse)
library(cowplot)

# outcome ----

# read in data
dat = read_excel('../../data/459 patients data for posting.xls')

ltach_aliases = LETTERS[1:length(unique(dat$LTACH))]
ltachs = unique(dat$LTACH)
dat$ltach = sapply(dat$LTACH, function(x){
  ltach_aliases[ltachs == x]
})

# subset data
dat = dat[dat$source != 'wound',]
dat = dat[order(dat$Cx_date),]
dat = dat[!duplicated(dat$Patient_ID),]
dat$infection = as.factor(dat$infection)

states = sapply(ltach_aliases, function(x) unique(dat$state[dat$ltach == x]))
                          
dat$ltach_state = paste0(dat$ltach, ' (',dat$state,')')

p1 = ggplot(dat,aes(source,fill=infection)) +
  geom_bar(color='black',position = position_dodge2(preserve='single',reverse=T,padding=0)) +
  scale_fill_manual(values=c('black','white'),
                    labels=c('Infection','Colonization'),
                    limits = c('yes','no'),
                    drop=F) + 
  labs(fill='') + 
  scale_x_discrete(labels = c('Blood','Respiratory','Urinary')) + 
  ylab('Count') + xlab('Anatomic site') +
  theme_bw() + theme(text = element_text(size=10))


#dat$ltach_state = factor(dat$ltach_state)
#levels(dat$ltach_state) = levels(dat$ltach_state)[order(gsub('^..','',levels(dat$ltach_state)))]
p2 = ggplot(dat, aes(x=infection,fill=source)) + 
  geom_bar(position = position_dodge2(preserve='single',reverse=F,padding=0)) + 
  facet_wrap(~ltach_state,nrow = 3) +
  scale_fill_manual(values=c('coral3','#2b8cbe','#feb24c'),
                    labels=c('Blood','Respiratory','Urinary')) + 
  scale_x_discrete(labels = c('Inf','Col'),
                   limits = c('yes','no')) + 
  labs(fill='Anatomic site') + xlab('') + ylab('Count') +
  theme_bw() + theme(text = element_text(size=10))

# pdf('results/figures/descriptive/outcome_barchart.pdf',width=10,height=5)
ggsave('results/figures/descriptive/outcome_barchart.tiff',p2,width=10,height=5)
#plot_grid(p1,p2,rel_widths = c(0.2),rel_heights = c(0.5))
# dev.off()

# features ----

# read in data

# files_collapsed = c("genomic/genomic_features.tsv",
#           "genomic/gene_features.tsv",
#           "genomic/kleborate_features.tsv",
#           "patient/patient_features.tsv")
# 
# # uncollapsed genomic
# snps = read.delim("data/features/genomic/snps.tsv")
# indels = read.delim("data/features/genomic/indels.tsv")
# insertions = read.delim("data/features/genomic/insertions.tsv")
# pangenome = read.delim("data/features/genomic/pangenome.tsv")
# genes = read.delim("data/features/genomic/genes.tsv")
# 
# # collapsed all
# genomic = read.delim("data/features/genomic/genomic_features.tsv")
# nrow(genomic)
# gene = read.delim("data/features/genomic/gene_features.tsv")
# nrow(gene)
# kleborate = read.delim("data/features/genomic/kleborate_features.tsv")
# 
# patient = read.delim("data/features/patient/patient_features.tsv")
# 
# to_keep = read.delim('data/to_keep.tsv',header = F)[,1]
# 
# snps = snps[!apply(snps,1,function(x) sum(x == 1) == 0 | sum(x == 1) == ncol(snps)),]
# indels = indels[!apply(indels,1,function(x) sum(x == 1) == 0 | sum(x == 1) == ncol(indels)),]
# insertions = insertions[!apply(insertions,1,function(x) sum(x == 1) == 0 | sum(x == 1) == ncol(insertions)),]
# pangenome = pangenome[!apply(pangenome,1,function(x) sum(x == 1) == 0 | sum(x == 1) == ncol(pangenome)),]
# genes = genes[!apply(genes,1,function(x) sum(x == 1) == 0 | sum(x == 1) == ncol(genes)),]
# 
# feature_cts_genomic = data.frame(Feature=
#                                    factor(c('SNVs','Small\nindels','Large\ninsertions','Pangenome','Grouped\nvariants'),
#                                           levels=c('SNVs','Small\nindels','Large\ninsertions','Pangenome','Grouped\nvariants')),
#                                  Count=c(nrow(snps),
#                                          nrow(indels),
#                                          nrow(insertions),
#                                          nrow(pangenome),
#                                          nrow(genes)))
# feature_cts_used = data.frame(Feature=factor(c('Genomic\n','Gene','Kleborate','Patient'),
#                                              levels=c('Genomic\n','Gene','Kleborate','Patient')),
#                               Count=c(sum(!duplicated(genomic)),
#                                       sum(!duplicated(gene)),
#                                       nrow(kleborate),
#                                       nrow(patient)))
# 
# print(feature_cts_used)
# 
# p3 = ggplot(feature_cts_genomic, aes(x=Feature,y=Count)) + 
#   geom_bar(stat='identity') + xlab('Individual genomic features') + 
#   theme_bw() + theme(text = element_text(size=15))
# 
# p4 = ggplot(feature_cts_used, aes(x=Feature,y=Count)) + 
#   geom_bar(stat='identity') + xlab('Feature sets') + 
#   theme_bw() + theme(text = element_text(size=15))
# 
# pdf('results/figures/descriptive/features_barchart.pdf',width=10,height=5)
# plot_grid(p3,p4,rel_widths = c(0.5))
# dev.off()
# 
# ocs = read.delim('data/outcomes/outcomes.tsv')
# sapply(ocs,table)
# table(ocs$infection,ocs$resp)
# table(ocs$infection,ocs$urine)
# table(ocs$infection,ocs$blood)
# 192+28
# 111+28
