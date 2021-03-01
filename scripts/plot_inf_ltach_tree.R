library(ape)
library(tidyverse)
library(ggtree)
library(phytools)
library(cowplot)

# read in tree
tree = midpoint.root(read.tree("../../../Sequence_data/consensus/2019_08_07_Penn_All_variant_calling/2020_02_13_16_29_09_core_results/gubbins/iqtree_results/2020_02_13_16_29_09_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile"))
ocs = read.delim('data/outcomes/outcomes.tsv')
ltachs = read.delim('data/ltachs.tsv')
to_keep = read.delim('data/to_keep.tsv',header = F)

ocs = ocs[to_keep$V1,]
ocs = ocs[order(match(rownames(ocs), to_keep$V1)),]
colnames(ocs) = c('Infection','Blood','Respiratory','Urinary')
src = sapply(1:nrow(ocs),function(x) names(ocs)[which(ocs[x,2:4] == 'yes')+1])
names(src) = rownames(ocs)
#src = src[to_keep$V1]
ltachs = ltachs[to_keep$V1,]
names(ltachs) = to_keep$V1
tree = drop.tip(tree,tree$tip.label[!tree$tip.label %in% to_keep$V1])

dat = data.frame(src=src,infection=ifelse(ocs$Infection=='yes','Infection','Colonization'))
#rownames(dat) = to_keep$V1

cols_mat = c('black','white','#2b8cbe','#feb24c','coral3')#"#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
names(cols_mat) = c('Infection','Colonization','Respiratory','Urinary','Blood')

dat2 = dat
dat2$src = NULL
dat2$infection = as.factor(dat2$infection)
th_inf = gheatmap(ggtree(tree,layout = 'fan'),dat2,width = 0.2,color = NA,colnames = F,offset = 0) + 
  scale_fill_manual(values=cols_mat) + labs(fill='Isolate type') + 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.key = element_rect(colour = 'black'))

dat3 = dat
dat3$infection = NULL
dat3$src = as.factor(dat3$src)
th_src = gheatmap(ggtree(tree,layout = 'fan'),dat3,width = 0.2,color = NA,colnames = F,offset = 0) + 
  scale_fill_manual(values=cols_mat)  + labs(fill='Anatomic site') + 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

gt = ggtree(tree,layout = 'circular') + geom_treescale(x=0.0013,offset=1.1,fontsize = 3) +
  theme_nothing()#x=-0.0003,y=-0.0003,fontsize=2,offset=1.1) #0.0013

th = gheatmap(gt,dat,width = 0.2,color = NA,colnames = F,offset = 0) + 
  scale_fill_manual(values=cols_mat) + theme(legend.position='none') + 
  theme_nothing()

leg <- plot_grid(
  NULL,
  get_legend(th_inf),
  get_legend(th_src),
  NULL,
             ncol=1, align = 'hv')

# th + annotation_custom(ggplotGrob(leg),xmin = 0.0002,xmax=0.0008,ymin=330)
             
#pdf('results/figures/heatmap/tree_inf_src.pdf',width=6,height=5)
thm = plot_grid(leg + theme(plot.margin=grid::unit(c(0,0,0,0),"mm")),
                NULL,
                th + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")),
                nrow=1,rel_widths = c(0.2,-0.15,0.8), rel_heights = c(0.1,0.9,0.9), greedy = T)

ggsave(thm,filename='results/figures/heatmap/tree_inf_src.pdf',width=6.87,height=5, units = 'in', dpi = 300)

ggsave(thm,filename='results/figures/heatmap/tree_inf_src.tiff',width=6.87,height=5, units = 'in', dpi = 300)
#dev.off()

