# descriptive plots (Figure S2)

# load packages
library(readxl)
library(tidyverse)
library(cowplot)

# outcome ----

# read in data
dat = read_excel('../../Regional_KPC_transmission/data/459 patients data for posting.xls')

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

p2 = ggplot(dat, aes(x=infection,fill=source)) + 
  geom_bar(position = position_dodge2(preserve='single',reverse=F,padding=0)) + 
  facet_wrap(~ltach_state,nrow = 3) +
  scale_fill_manual(values=c('coral3','#2b8cbe','#feb24c'),
                    labels=c('Blood','Respiratory','Urinary')) + 
  scale_x_discrete(labels = c('Inf','Col'),
                   limits = c('yes','no')) + 
  labs(fill='Anatomic site') + xlab('') + ylab('Count') +
  theme_bw() + theme(text = element_text(size=10))

pdf('results/figures/descriptive/outcome_barchart.pdf',width=10,height=5)
p2
dev.off()

