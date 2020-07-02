# Plot tree with genomic features (e.g. ICEKp10; Figures 4, S8, S10)

# tree + heatmap of interesting genomic features
library(phytools)
library(ggplot2)
library(ggtree)
library(cowplot)
library(grid)
library(readxl)

# read in tree
tree = midpoint.root(read.tree("../../../Sequence_data/consensus/2019_08_07_Penn_All_variant_calling/2020_02_13_16_29_09_core_results/gubbins/iqtree_results/2020_02_13_16_29_09_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile"))

# read in matrix
mat = data.frame(t(read.delim('data/combined/infection_patient-kleborate.tsv')))
mat = mat[mat$ST == 'ST258',]
tip_order = tree$tip.label[tree$tip.label %in% rownames(mat)]
tree = midpoint.root(drop.tip(tree,tree$tip.label[!tree$tip.label %in% rownames(mat)]))
mat = mat[tree$tip.label,]

# prettify certain rows for plotting
mat$Bla_Carb[mat$Bla_Carb == 'KPC-3*'] = 'KPC-3'
mat$Bla_Carb[mat$Bla_Carb != 'KPC-2' & mat$Bla_Carb != 'KPC-3'] = NA
mat$O_locus_missing_genes.O1.O2v2_07_kfoC = ifelse(mat$O_locus_missing_genes.O1.O2v2_07_kfoC == 1,'kfoC absent',NA)
mat$Yersiniabactin.ybt_0 = ifelse(mat$Yersiniabactin.ybt_0 == 1, 'ybt0 present',NA)
mat$Yersiniabactin.ybt_17 = ifelse(mat$Yersiniabactin.ybt_17 == 1, 'ybt17 present',NA)
mat$Yersiniabactin._ICEKp10 = ifelse(mat$Yersiniabactin._ICEKp10 == 1, 'ICEKp10 present',NA)
mat$Colibactin = ifelse(mat$Colibactin != '-', 'Colibactin present',NA)
mat$infection = ifelse(mat$infection == 'yes','Infection',NA)

ocs = read.delim('data/outcomes/outcomes.tsv')
src = sapply(1:nrow(ocs), function(x){
  names(ocs)[which(ocs[x,c('blood','resp','urine')] == 'yes')+1]
})
names(src) = rownames(ocs)
src = src[rownames(mat)]

ltach = read.delim('data/ltachs.tsv')
ltach = ltach[rownames(mat),]
names(ltach) = rownames(mat)

pt_info = read_excel('../../Regional_KPC_transmission/data/459_patients_data_for_posting.xls')

state = pt_info$state
names(state) = paste0('PCMP_H','',pt_info$isolate_no)
state = state[rownames(mat)]

mat_sub = data.frame(
  kpc=mat$Bla_Carb,
  src=ifelse(src=='resp','Respiratory',NA),
  inf=mat$infection,
  blank2=rep(NA,nrow(mat)),
  kfoCmiss=mat$O_locus_missing_genes.O1.O2v2_07_kfoC,
  icekp10=mat$Yersiniabactin._ICEKp10,
  ybt_0=mat$Yersiniabactin.ybt_0,
  ybt_17=mat$Yersiniabactin.ybt_17,
  colibactin=mat$Colibactin
)

cols_mat = c(blood='coral3',Respiratory='#2b8cbe',urine='#feb24c','KPC-2'='gray77','KPC-3'='gray87','kfoC absent'='violetred3','ybt17 present'='aquamarine2','ybt0 present'='lightseagreen','Infection'='black','ICEKp10 present'='lightsteelblue4','Colibactin present' = 'lightslateblue')

ltach_cols = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
names(ltach_cols) = unique(ltach)

tip_cols = ltach_cols[ltach]

tips_kpc = mat_sub$kpc
tip_cols_kpc = as.character(tips_kpc)
tip_cols_kpc[tip_cols_kpc=='KPC-2'] = 'gray87'
tip_cols_kpc[tip_cols_kpc=='KPC-3'] = 'gray57'
mat_sub$kpc = NULL

tr = ggtree(tree,layout='rectangular')  +
  geom_balance(node=335, fill='cadetblue3', color=NA, alpha=0.3, extend=0.00027) +
  geom_tippoint(color=tip_cols_kpc) 

tr_ltach = ggtree(tree,layout='rectangular')  +
  geom_balance(node=335, fill='cadetblue3', color=NA, alpha=0.3, extend=0.00001) +
  geom_tippoint(color=tip_cols)

tr_i = gheatmap(tr,data.frame(mat_sub$inf)) + scale_fill_manual(values=cols_mat)

names(mat_sub) = c('Respiratory','Infection',' ','kfoC disrupted','ICEKp10 present','ybt0 present','ybt17 present','Colibactin present')

mat_main = mat_sub[,c('Respiratory', 'Infection', 'kfoC disrupted', 'ICEKp10 present')]
mat_main$`ICEKp10 present` = NA
mat_main$`ICEKp10 present`[mat_sub$`ybt0 present` == 'ybt0 present'] = 'ybt0 present'
mat_main$`ICEKp10 present`[mat_sub$`ybt17 present` == 'ybt17 present'] = 'ybt17 present'

mat_ice = mat_sub[,c('ICEKp10 present','Colibactin present','ybt17 present','ybt0 present')]
colnames(mat_ice) = gsub(' present','',colnames(mat_ice))


leg_kpc = get_legend(ggplot(mat[!is.na(mat$Bla_Carb),],aes(x=Bla_Carb,y=as.numeric(num_resistance_classes),col=Bla_Carb)) +
  geom_point(shape=19,size=3) + theme_bw() +
  scale_color_manual(values=c('KPC-3'='gray77','KPC-2'='gray87'))+
    guides(color = guide_legend(title.position = "top",
                                title.hjust = 0.5,
                                label.position = "top")) +
    theme(legend.title = element_blank(),
          legend.position="top",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10))
  )

leg_ybt = get_legend(ggplot(mat_main[!is.na(mat_main$`ICEKp10 present`),],aes(x=`ICEKp10 present`,y=(Respiratory),col=`ICEKp10 present`)) +
                       geom_point(shape=15,size=5) + theme_bw() +
                       scale_color_manual(values=c('ybt0 present'='lightseagreen','ybt17 present'='aquamarine2'),
                                          labels=c('ybt0','ybt17'))+
                       guides(color = guide_legend(title.position = "top",
                                                   title.hjust = 0.5,
                                                   label.position = "top")) +
                       theme(legend.title = element_blank(),
                             legend.position="top",
                             legend.margin=margin(0,0,0,0),
                             legend.box.margin=margin(-10,-10,-10,-10))
)

leg_ltach = get_legend(ggplot(data.frame(LTACH=ltach),aes(x=as.numeric(LTACH),y=as.numeric(LTACH),col=LTACH)) + geom_point(size=3) +
                       theme_bw() +
                       scale_color_manual(values=tip_cols) +
                         guides(color = guide_legend(title.position = "top",
                                                     title.hjust = 0.5,
                                                     label.position = "top",
                                         ncol=3)) +
                       theme(
                             legend.position="top",
                             legend.margin=margin(0,0,0,0),
                             legend.box.margin=margin(-10,-10,-10,-10))
)

# disrupted kfoC lineage: node 335
# ybt0 lineage: node 407
tr_ann = tr +
  geom_treescale(x=0.0002,y=355,offset = 3) + #260
  geom_cladelabel(node=335, "disrupted kfoC\nlineage", offset = -0.00065, hjust = 'right',barsize = NA) #+ 

tr_ltach_ann = tr_ltach + 
  geom_treescale(x=0.0002,y=355,offset = 3) + #260
  geom_cladelabel(node=335, "disrupted kfoC\nlineage", offset = -0.00065, hjust = 'right',barsize = NA) #+ 

thm = gheatmap(tr_ann,
         mat_main,
         color = NA, offset = 0,width = 0.3, colnames = T, colnames_angle = 45,
         colnames_offset_x = 0.000001, colnames_offset_y = 10, hjust = 0,
         colnames_position = 'top') + 
  scale_fill_manual(values=cols_mat) + 
  theme(legend.position = "none") + 
  ylim(0, 390) + xlim(0, 0.0012)

names(mat_ice)[3] = 'Yersiniabactin (ybt17)'
names(mat_ice)[4] = 'Yersiniabactin (ybt0)'

thm_ice = gheatmap(tr_ann,
               mat_ice,
               color = NA, offset = 0,width = 0.3, colnames = T, colnames_angle = 45,
               colnames_offset_x = 0.000001, colnames_offset_y = 10, hjust = 0,
               colnames_position = 'top') + 
  scale_fill_manual(values=cols_mat) + 
  theme(legend.position = "none") + 
  ylim(0, 390) + xlim(0, 0.0012)

legs = plot_grid(leg_kpc,NULL,ncol=1,rel_heights = c(0.85,0.15))

pdf('results/figures/heatmap/tree_icekp10_heatmap.pdf')
thm_ice + annotation_custom(ggplotGrob(legs),xmin = 0.0002,xmax=0.0008,ymin=300)
dev.off()

legs_both = plot_grid(leg_kpc,leg_ybt,nrow=1)

pdf('results/figures/heatmap/tree_heatmap.pdf')
thm + annotation_custom(ggplotGrob(legs_both),xmin = 0.00032,xmax=0.0008,ymin=315)
dev.off()

pdf('results/figures/heatmap/tree_ltachs.pdf')
plot_grid(tr_ltach_ann, leg_ltach, rel_widths = c(0.8,0.2))
dev.off()


