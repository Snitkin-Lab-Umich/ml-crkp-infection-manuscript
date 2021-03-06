library(tidyverse)
library(cowplot)
library(grid)
library(ggplotify)
source('scripts/geom_boxplot2_base_r_look.R')
source('scripts/prettify_names.R')

# plot_diffs = function(diffs,outfile,cutoff=0.25){
#   p = diffs %>% pivot_longer(cols = everything()) %>% 
#     ggplot(aes(x=reorder(name, value, FUN=quantile,prob=0.25),y=value)) +
#     geom_boxplot2() +
#     ylab('Test AUROC - mean permuted AUROC') + xlab('') + theme_bw() +
#     coord_flip() +
#     geom_hline(yintercept = 0, col='red') +
#     geom_vline(xintercept = ncol(diffs) - n_above_zero(diffs,cutoff=cutoff) + 0.5, lty=3)
#     #theme(axis.text.x = element_text(angle = 20, hjust = 0, vjust = 0)) #hjust=1)) #,hjust = 0, vjust = 0,)) + 
#   ggsave(plot = p,filename = outfile, width = 10,height = 10)
# }

n_above_zero = function(diffs,cutoff=0.25){
  diffs %>% #pivot_longer(cols = everything()) %>% group_by(name) %>% 
    dplyr::summarize(first=quantile(value,probs = cutoff)) %>% #,median=median(value),mean=mean(value)) %>%
    arrange(desc(first)) %>% filter(first>0) %>% nrow() 
    #print(n=Inf) 
}

files = list.files('results/permutation_importance',pattern='1_1.tsv',full.names = T)
diffs_all = lapply(files, function(x){
  d = read_tsv(x)
  colnames(d) = prettify_names(colnames(d))
  return(d)
  })
names(diffs_all) = gsub('.*/|_1_1.tsv','',files)

assn = as_tibble(read.delim('results/median_importances.tsv'))
assn$feature = prettify_names(assn$feature)
assn$dat = gsub('_lr','',assn$dat)

patient = unique(c(sapply(names(diffs_all), function(x) prettify_names(colnames(get_model_matrix(paste0('data/combined/',gsub('-kleborate','',x),'.tsv')))))))

# summarize_data = function(d, cutoff=0.25){
#   diffs_all[[d]] %>%
#     pivot_longer(cols = everything()) %>% group_by(name) %>% 
#     dplyr::summarize(first=quantile(value,probs = 0, cutoff),mean=mean(value)) %>% 
#     mutate(dat=d)
# }

make_long = function(d){
  diffs_all[[d]] %>%
    pivot_longer(cols = everything()) %>% 
    group_by(name) %>% 
    mutate(dat=d)
}

alldat = bind_rows(lapply(names(diffs_all), make_long))

alldat$pos_frac = sapply(1:nrow(alldat), function(x){
  if(alldat$name[x] %in% assn$feature){
    pos_frac = assn$pos_frac[assn$feature == alldat$name[x] & assn$dat == alldat$dat[x]]
  }else if(sum(grepl('StrA',alldat$name[x]))){
    pos_frac = assn$pos_frac[assn$feature == 'Aminoglycoside res (StrA)' & assn$dat == alldat$dat[x]]
  }else if(sum(grepl('DfrA17',alldat$name[x]))){
    pos_frac = assn$pos_frac[assn$feature == 'Trimethoprim res (DfrA17)' & assn$dat == alldat$dat[x]]
  }else{
    return(alldat$name[x])
  }
  return(pos_frac)
})

alldat$assn = sapply(alldat$pos_frac, function(x){
  if(x > 0.9) y = 'Infection'
  else if(x < 0.1) y = 'Colonization'
  else y = NA
  return(y)
})

alldat$patient = ifelse(alldat$name %in% patient,'Patient','Genomic')

plot_diffs = function(alldat,dat_name, cutoff=0.25){
  outfile = paste0('results/figures/importances/',dat_name,'_perm-imp.tiff')
  dat_sub = alldat %>% filter(dat == dat_name)
  p = dat_sub %>% 
    ggplot(aes(x=reorder(name, value, FUN=quantile,prob=cutoff),y=value,fill=assn)) +
    geom_boxplot2() +
    ylab('Test AUROC - mean permuted AUROC') + xlab('') + theme_bw() +
    coord_flip(clip='off') +
    geom_hline(yintercept = 0, col='darkgrey',size=1) +
    geom_vline(xintercept = length(unique(dat_sub$name)) - n_above_zero(dat_sub,cutoff) + 0.5, lty=3) +
    scale_fill_manual(values = c(Colonization='darkslategray3',Infection='#756bb1')) +
    theme(axis.text.y = element_text(margin = margin(r=10)),legend.position = 'none') + 
    ylim(c(min(alldat$value)-0.01,max(alldat$value)+0.01))
  p_leg = get_legend(dat_sub %>% 
    ggplot(aes(x=reorder(name, value, FUN=quantile,prob=cutoff),y=value,fill=assn)) +
    geom_boxplot2() + theme_bw() + 
    scale_fill_manual(values = c(Colonization='darkslategray3',Infection='#756bb1')) +
    labs(fill = "Association"))
  p_shape = dat_sub %>%
    ggplot() + geom_point(aes(x=0,y=reorder(name, value, FUN=quantile,prob=cutoff),shape=patient)) + 
    scale_shape_manual(values = c(19,15)) + coord_cartesian(xlim = c(-0.06,0.05)) + #c(-0.0405,0.05)) +
    theme_void() + theme(legend.position = 'none') 
  p_shape_leg = get_legend(dat_sub %>%
    ggplot() + geom_point(aes(x=0,y=reorder(name, value, FUN=quantile,prob=cutoff),shape=patient),size=3) + 
    scale_shape_manual(values = c(19,15)) + theme_bw() + labs(shape = "Feature type"))
  #theme(axis.text.x = element_text(angle = 20, hjust = 0, vjust = 0)) #hjust=1)) #,hjust = 0, vjust = 0,)) 
  comb = ggdraw(p) + draw_plot(p_shape,y = 0.018, scale = 0.948)
  
  comb_and_leg = plot_grid(comb,
                           plot_grid(NULL,p_leg,p_shape_leg,NULL,ncol=1, rel_heights = c(0.38,0.12,0.12,0.38)),
                           nrow=1,rel_widths = c(0.8,0.2))
  
  ggsave(plot = comb_and_leg,filename = outfile, width = 10,height = 10)
  return(comb_and_leg)
}

# for(dat_name in names(diffs_all)){
#   print(dat_name)
#   plot_diffs(alldat,dat_name)
# }

p_all <- plot_diffs(alldat, names(diffs_all)[1]) 
p_resp <- plot_diffs(alldat, names(diffs_all)[2]) 
p_urine <- plot_diffs(alldat, names(diffs_all)[3])

ggsave('results/figures/importances/perm_imp_comb.pdf',
       plot_grid(p_all, p_resp, p_urine, labels = c('A Overall','B Resp','C Urinary'),ncol = 3), width = 25, height = 10)

ggsave('results/figures/importances/perm_imp_comb.tiff',
       plot_grid(p_all, p_resp, p_urine, labels = c('A Overall','B Resp','C Urinary'),ncol = 3), width = 25, height = 10)




alldat$name = trimws(alldat$name)

sumdat = alldat %>% group_by(name,dat) %>% 
  dplyr::summarize(
    o05=quantile(value,probs = 0.05),
    first=quantile(value,probs = 0.25),
    third=quantile(value,probs = 0.75),
    o95=quantile(value,probs = 0.95),
    mean=mean(value),
    pos_frac=unique(pos_frac),
    assn=unique(assn),
    patient=unique(patient))

feats = sumdat %>% filter(first > 0) %>% select(name) %>% unique() %>% unlist()

sumdat = sumdat %>% filter(name %in% feats) %>% filter(!is.na(assn)) 

sumdat$fill = as.factor(sapply(1:nrow(sumdat), function(n){
  infection = sumdat$pos_frac[n] > 0.5
  d1 = sumdat$dat[n] == 'infection-resp_patient-kleborate'
  d2 = sumdat$dat[n] == 'infection-urine_patient-kleborate'
  d3 = sumdat$dat[n] == 'infection_patient-kleborate'
  if(is.na(d1) | is.na(infection)) {NA}
  else if(d1 & infection) {'d1inf'}
  else if(d1 & !infection) {'d1col'}
  else if(d2 & infection) {'d2inf'}
  else if(d2 & ! infection) {'d2col'}
  else if(d3 & infection) {'d3inf'}
  else if(d3 & !infection) {'d3col'}
}))

sumdat$r_imp = (sumdat$first > 0) & sumdat$dat == 'infection-resp_patient-kleborate'
sumdat$u_imp = (sumdat$first > 0) & sumdat$dat == 'infection-urine_patient-kleborate'
sumdat$all_imp = (sumdat$first > 0) & sumdat$dat == 'infection_patient-kleborate'

sumdat$grp = sapply(sumdat$name, function(x){
  if(sum(sumdat$u_imp[sumdat$name == x]) == 1 & sum(sumdat$r_imp[sumdat$name == x]) == 0){
    return('5u')
  }else if(sum(sumdat$r_imp[sumdat$name == x]) == 1 & sum(sumdat$u_imp[sumdat$name == x]) == 0){
    return('1r')
  }else if(sum(sumdat$r_imp[sumdat$name == x]) == 1 & sum(sumdat$u_imp[sumdat$name == x]) == 1){
    if(sumdat$mean[sumdat$name == x & sumdat$r_imp] > sumdat$mean[sumdat$name == x & sumdat$u_imp]){
      return('2ru')
    }else{
      return('4ur')
    }
  }else{ #if(sum(!d$r_imp[d$feature == x]) == 3 & sum(!d$u_imp[d$feature == x]) == 3){
    return('3a')
  }
})
# %>% arrange(grp,desc(diff))

# diffs = (sumdat$mean[sumdat$dat == 'infection-resp_patient-kleborate'] - 
#            sumdat$mean[sumdat$dat == 'infection-urine_patient-kleborate'])
# sumdat$diff[sumdat$dat == 'infection-resp_patient-kleborate'] = diffs
# sumdat$diff[sumdat$dat == 'infection-urine_patient-kleborate'] = diffs

sumdat$dat = factor(sumdat$dat,levels=c('infection-resp_patient-kleborate','infection-urine_patient-kleborate','infection_patient-kleborate'))
sumdat = sumdat %>% filter(!is.na(assn)) %>% group_by(grp) %>% arrange(grp,desc(mean))#,desc(diff))
sumdat$name = factor(sumdat$name,levels=unique(sumdat$name))

sumdat = sumdat %>% mutate(imp=first > 0)


point_cols = c('#2b8cbe','#feb24c','grey')
point_labs = c('Respiratory','Urinary','All')
outcome_labels=c('Colonization','Infection')

#
p = sumdat %>% filter(first > 0) %>% ggplot(aes(mean,name,col=dat,shape=patient,fill=fill)) + #,alpha=imp)) + 
  xlab('Mean difference between test\nand permuted AUROC') + ylab('') +
  geom_vline(aes(xintercept=0)) + 
  geom_errorbarh(aes(xmin=o05, xmax=o95),height = 0,linetype = "dotted",alpha=0.5,size=1) +
  geom_errorbarh(aes(xmin=first, xmax=third),height = 0,alpha=0.5,size=1) +#,linetype = "dotted") +
  geom_point(size=3,stroke=1) +
  scale_shape_manual(values=c(21,22)) + 
  xlim(c(-0.04,0.12)) +
  theme_bw() + 
  theme(legend.position='left', axis.title.x = element_text(size=10)) +
  scale_y_discrete(position = "right") + 
  scale_color_manual(values = point_cols, labels = point_labs) +
  scale_fill_manual(values = c('white',point_cols[1],'white',point_cols[2],'white',point_cols[3]),
                    labels = outcome_labels, 
                    breaks = c('d1col','d1inf')
                    ) +
  #scale_alpha_discrete(range = c(0.5, 1)) + 
  guides(fill = guide_legend(override.aes=list(shape=21,fill=c('white','black')),order = 3, ncol = 1, title.position = 'top'),
         color = guide_legend(order=1, ncol = 1, title.position = 'top'), 
         shape = guide_legend(order=2, ncol = 1, title.position = 'top'),
         linetype = guide_legend(ncol = 1, title.position = 'top')) + 
  labs(shape='Feature type (shape)',fill='Associated with (fill)',color='Anatomic site (color)') +
  geom_hline(aes(yintercept=8.5), linetype='dashed') +
  geom_hline(aes(yintercept=18.5), linetype='dashed') + geom_line(aes(0,0,linetype=all_imp)) +
  #labs(linetype='Confidence interval') +
  scale_linetype_discrete(name = "Quantile", labels = c('25%-75%','5%-95%')) #+
  # theme(legend.position = 'top',legend.justification = "right",
  #       legend.key.size = unit(5,'mm'))
  

# Create text
ugrob <- grobTree(textGrob("Associated in\nurinary analysis", x=-0.01,  y=0.91, hjust=1,
                           gp=gpar(fontsize=11, fontface="italic")))
rgrob <- grobTree(textGrob("Associated in\nrespiratory analysis", x=-0.01,  y=0.18, hjust=1,
                           gp=gpar(fontsize=11, fontface="italic")))
agrob <- grobTree(textGrob("Associated in\noverall analysis or \nboth site-specific analyses", x=-0.01,  y=0.245, hjust=1,
                           gp=gpar(fontsize=11, fontface="italic")))
# Plot
p1 = p + annotation_custom(ugrob) + annotation_custom(rgrob) + #+ annotation_custom(agrob) 
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 100), "points"))

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(p1))
gt$layout$clip[gt$layout$name == "panel"] <- "off"

# pdf('results/figures/importances/urine-resp_perm-imps.pdf',height = 10, width = 10)
# ggsave('results/figures/importances/urine-resp_perm-imps.pdf',as.ggplot(gt),height = 5, width = 8.1)
ggsave('results/figures/importances/urine-resp_perm-imps.pdf',
       p,#as.ggplot(gt),
       height = 5, width = 6.87, units = 'in', dpi = 300)

ggsave('results/figures/importances/urine-resp_perm-imps.tiff',
       p,#as.ggplot(gt),
       height = 5, width = 6.87, units = 'in', dpi = 300)
# dev.off()

