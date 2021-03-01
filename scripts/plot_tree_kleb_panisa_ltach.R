library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(cowplot)

source('scripts/panISa_parser_functions.R')



tree = read.tree("../../../Sequence_data/consensus/2019_08_07_Penn_All_variant_calling/2020_02_13_16_29_09_core_results/gubbins/iqtree_results/2020_02_13_16_29_09_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile")

kleb = data.frame(t(read.delim('data/features/genomic/kleborate_features.tsv')))

tree = midpoint.root(drop.tip(tree,tree$tip.label[!tree$tip.label %in% rownames(kleb)]))

kleb = kleb[tree$tip.label,]
kfoc = kleb$O_locus_missing_genes.O1.O2v2_07_kfoC
names(kfoc) <- rownames(kleb)
oloc = kleb$O_locus

kfoc_start = 422223
kfoc_end = 423431

isfinder = read.delim('../2020-04-01_kfoC-panisa/results/ISFinder.txt', stringsAsFactors = F)
isfinder = isfinder[isfinder$Sample != 'Sample',]
isfinder = isfinder[isfinder$Chromosome == 'CP031810',]

ismat = make_ismat(isfinder)
ismat = ismat[,names(kfoc)]

ismat_annots = data.frame(t(sapply(strsplit(rownames(ismat),'_'), function(x) sort(as.numeric(x)))))
names(ismat_annots) = c('start','end')

is_kfoc = ismat[ismat_annots$start >= kfoc_start & ismat_annots$end <= kfoc_end,]

is_kfoc_one = is_kfoc[,colSums(is_kfoc) == 1]
is_kfoc_one = is_kfoc_one[rowSums(is_kfoc_one) > 0,]

# keep only ones where it's unique to a sample
is_kfoc = is_kfoc[rownames(is_kfoc_one),]

is_elem = ifelse(colSums(is_kfoc) == 1 & oloc == 'O2v2', 1, 0)
is_elem = is_elem[tree$tip.label]

# Get IS element
pos = data.frame(t(sapply(1:nrow(isfinder), function(x){
  start = min(isfinder$Start_Position[x],isfinder$Stop_Position[x])
  end =  max(isfinder$Start_Position[x],isfinder$Stop_Position[x])
  return(as.numeric(c(start,end)))
})))
colnames(pos) = c('start','end')
#isfinder_kfoc = isfinder[pos$start >= kfoc_start & pos$end <=kfoc_end,]
is_kfoc_1 = isfinder[isfinder$Start_Position == 422869 & isfinder$Stop_Position == 422873,]
rownames(is_kfoc_1) = gsub('_panISa','',is_kfoc_1$Sample)
table(is_kfoc_1$Potential_IS,is_kfoc_1$Alignment)

# is_kfoc_1_bin = sapply(names(kfoc), function(x) as.numeric(x %in% rownames(is_kfoc_1)))
potential_is = is_kfoc_1$Potential_IS
names(potential_is) = rownames(is_kfoc_1)
potential_is_all = sapply(names(kfoc), function(x){
  ifelse(x %in% names(is_elem)[is_elem==1],ifelse(x %in% names(potential_is),potential_is[x],'No identity'),NA)
} )

# blast
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
blast = read.delim('../2020-03-13_blast-kfoc/kfoc_blast_out.tsv', header = F, stringsAsFactors = F)
blast = blast[!duplicated(blast$V1),]

samps = gsub('_[0-9].*','',blast$V1)
samp_ct_subset = table(samps)
samp_ct_subset = samp_ct_subset[names(kfoc)]
names(samp_ct_subset) = names(kfoc)
nas = names(samp_ct_subset)[is.na(samp_ct_subset)]
table(kfoc,samp_ct_subset)
breaks = blast[samps %in% names(samp_ct_subset)[samp_ct_subset == 2],]
breaks_samps = unique(gsub('_[0-9].*','',breaks$V1))

blast_missing = !tree$tip.label %in% unique(gsub('_[0-9].*','',blast$V1))
blast_break = tree$tip.label %in% breaks_samps



dat = data.frame(
  cbind(as.character(kleb$O_locus),
        #ifelse(as.character(kleb$O_locus)=='O2v2','O2v2',NA),
        ifelse(kfoc==1,'kfoC \"absent\"',NA),
        #ifelse(blast_break, 'kfoC split',NA),
        ifelse(blast_missing, 'kfoC absent',ifelse(blast_break, 'kfoC split',NA)),
        potential_is_all))
#ifelse(is_elem==1,'panISa: IS element in kfoC',NA)))

cols_mat = c('kfoC \"absent\"'='violetred3',
             'kfoC absent'='darkgrey',
             'kfoC split'='lightgrey',
             'IS1S'='aquamarine3','IS1X2'='aquamarine2',
             'No identity'='aquamarine',
             'O2v2'='lightblue1','O2v1'='lightblue2','O3b'='lightblue3','OL102'='lightblue4')


# kleb$Bla_Carb[kleb$Bla_Carb == 'KPC-3*'] = 'KPC-3'
# kleb$Bla_Carb[kleb$Bla_Carb != 'KPC-2' & kleb$Bla_Carb != 'KPC-3'] = NA
# tips_kpc = kleb$Bla_Carb
# tip_cols_kpc = as.character(tips_kpc)
# tip_cols_kpc[tip_cols_kpc=='KPC-2'] = 'gray87'
# tip_cols_kpc[tip_cols_kpc=='KPC-3'] = 'gray57'

tr = ggtree(tree,layout='rectangular')  #+
#geom_balance(node=335, fill='cadetblue3', color=NA, alpha=0.3, extend=0.00035) +
#geom_tippoint(color=tip_cols_kpc) 

tr_ann = tr +
  geom_treescale(x=0.0001,y=300,offset = 3) #+ 
#geom_cladelabel(node=335, "disrupted kfoC\nlineage", offset = -0.00065, hjust = 'right',barsize = NA)

hmo_leg = get_legend(gheatmap(tr_ann,dat[,1, drop = FALSE],
                              color = NA, offset = 0,width = 0.4, colnames = F) + 
                       scale_fill_manual(values=cols_mat,breaks=names(cols_mat),name='O locus'))

hmk_leg = get_legend(gheatmap(tr_ann,dat[,2, drop = FALSE],
                              color = NA, offset = 0,width = 0.4, colnames = F) + 
                       scale_fill_manual(values=cols_mat,breaks=names(cols_mat),name='Kleborate'))

hmb_leg = get_legend(gheatmap(tr_ann,dat[,3, drop = FALSE],
                              color = NA, offset = 0,width = 0.4, colnames = F) + 
                       scale_fill_manual(values=cols_mat,breaks=names(cols_mat),name='Assembly blast'))

hmp_leg = get_legend(gheatmap(tr_ann,dat[,4, drop = FALSE],
                              color = NA, offset = 0,width = 0.4, colnames = F) + 
                       scale_fill_manual(values=cols_mat,breaks=names(cols_mat),name='IS element'))


ltach = read.delim('data/ltachs.tsv')
ltach = ltach[rownames(mat),]
names(ltach) = rownames(mat)

ltach_cols = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
names(ltach_cols) = unique(ltach)

tip_cols = ltach_cols[ltach]

tr_ltach = ggtree(tree,layout='rectangular')  +
  geom_balance(node=335, fill='cadetblue3', color=NA, alpha=0.3, extend=0.00001) +
  geom_tippoint(color=tip_cols)

leg_ltach = get_legend(ggplot(data.frame(LTACH=ltach),aes(x=as.numeric(LTACH),y=as.numeric(LTACH),col=LTACH)) + geom_point(size=3) +
                         theme_bw() +
                         scale_color_manual(values=tip_cols) +
                         guides(color = guide_legend(title.position = "top",
                                                     # hjust = 0.5 centres the title horizontally
                                                     title.hjust = 0.5,
                                                     label.position = "top",
                                                     ncol=3)) +
                         theme(#legend.title = element_blank(),
                           legend.position="top",
                           legend.margin=margin(0,0,0,0),
                           legend.box.margin=margin(-10,-10,-10,-10))
)

tr_ltach_ann = tr_ltach + 
  geom_treescale(x=0.00005,y=320,offset = 3) + #260
  geom_cladelabel(node=335, "disrupted kfoC\nlineage", offset = -0.00065, hjust = 'right',barsize = NA) #+ 
#geom_cladelabel(node=407, "ybt0\nlineage", offset = 0.00001, hjust = 'left') 

thm = gheatmap(tr_ltach_ann,dat,
               color = NA, offset = 0,width = 0.4, colnames = F) + 
  scale_fill_manual(values=cols_mat,breaks=names(cols_mat),name='') + theme(legend.position = 'none')


# pdf('results/figures/heatmap/tree_icekp10_heatmap.pdf')
thm_ice + annotation_custom(ggplotGrob(legs),xmin = 0.00075,xmax=0.00155,ymin=-100)
# dev.off()

p1 <- plot_grid(thm,plot_grid(NULL,hmo_leg,hmk_leg,hmb_leg,hmp_leg,NULL,
                              ncol=1,align = 'hv',rel_heights = c(0.2,0.15,0.15,0.15,0.15,0.2)), plot_grid(leg_ltach),
                nrow=1,rel_widths = c(0.7,0.3, 0.1))


# pdf('results/figures/heatmap/tree_kleb_panisa_ltach.pdf', height = 8, width = 8)
# ggsave('results/figures/heatmap/tree_kleb_panisa_ltach.pdf',
#   plot_grid(thm,plot_grid(NULL,hmo_leg,hmk_leg,hmb_leg,hmp_leg,NULL,
#                         ncol=1,align = 'hv',rel_heights = c(0.2,0.15,0.15,0.15,0.15,0.2)), plot_grid(leg_ltach),
#           nrow=1,rel_widths = c(0.7,0.3, 0.1)), width = 10, height = 10
# )
# dev.off()


library(ape)
library(readxl)
library(ggtree)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

source('scripts/get_isolate_locations.R')


#tree = read.nexus('../../Regional_KPC_transmission/2017-08-31_beast-plp-pure-subset-v2-starting-tree/plp-pure-subset-v2_gtr_ucln_bs_st_resampled500000.comb.mcc.tree') #read.tree('../../../Sequence_data/consensus/2018_05_29_Penn_ST258_Limbago_Patric_variant_calling/2019_04_29_09_15_09_core_results/gubbins/2018_09_11_09_18_09_KPNIH1_ref_allele_var_consensus.treefile')
tree = read.tree('../../../Sequence_data/consensus/2018_05_29_Penn_ST258_Limbago_Patric_variant_calling/2020_03_24_15_40_08_core_results/gubbins/iqtree_results/2020_03_24_15_40_08_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile')

penn_data_file = "../../data/459 patients data for posting.xls"
patric_data_file = "../../data/2016-9-1_PATRIC_kleb_genomes_meta_data.txt"
limbago_data_file = '../../data/Limbago_metadata.txt'
srr_biosample_file = '../../data/penn-limbago-patric_srr-biosample.txt'

# public data
limbago_data = read.delim(limbago_data_file)
patric_data = read.delim(patric_data_file)
#to convert srr to biosample
srr_biosample = read.delim(srr_biosample_file,header = F,sep = ' ',col.names = c('SRR','BioSamp'))
# penn data
sample_key = read_excel(penn_data_file)

sample_key$isolate_no = paste0('PCMP_H',sample_key$isolate_no)

# sample_dates = sapply(sample_key[, 'Isolate_no'], FUN = function(x){sample_key[sample_key[, 'Isolate_no'] == x , 'culture_date']});

# get locations
# collapse ltachs
sample_key$LTACH_col = sample_key$state #as.character(sample_key$LTACH)
#sample_key$LTACH_col[sample_key$state != 'CA'] = as.character(sample_key$state[sample_key$state != 'CA'])


mat = data.frame(t(read.delim('data/combined/infection_patient-kleborate.tsv')))
mat = mat[mat$ST == 'ST258',]

kleb_public = read.delim('../2020-02-23_kleborate-public/kleborate_public_genomes.txt')

ybt0_pub = ifelse(grepl('ICEKp10',kleb_public$Yersiniabactin),'red',NA)
names(ybt0_pub) = kleb_public$strain

ybt0_penn = ifelse(mat$Yersiniabactin._ICEKp10 == 1, 'red',NA)
names(ybt0_penn) = rownames(mat)

ybt0 = c(ybt0_penn,ybt0_pub)
tree = drop.tip(tree,tree$tip.label[!tree$tip.label %in% names(ybt0)])
ybt0 = ybt0[tree$tip.label]

omiss_pub = ifelse(grepl('O1/O2v2_07_kfoC',kleb_public$O_locus_missing_genes),'red',NA)
names(omiss_pub) = kleb_public$strain

omiss_penn = ifelse(mat$O_locus_missing_genes.O1.O2v2_07_kfoC == 1, 'red',NA)
names(omiss_penn) = rownames(mat)

omiss = c(omiss_penn,omiss_pub)
omiss = omiss[tree$tip.label]

locs = get_locations(tree$tip.label)
names(locs) = gsub('\\..*','',names(locs))

cols = brewer.pal(n=length(unique(locs))+2,name='Set3')
names(cols) = c(sort(unique(locs)),'kfoC disrupted','ICEKp10 present')

ybt0_mat = ifelse(ybt0=='red','ICEKp10 present',NA)
omiss_mat = ifelse(omiss=='red','kfoC disrupted',NA)

tgen = get_legend(gheatmap(ggtree(tree),data.frame(ybt0_mat,omiss_mat)) +
                    scale_fill_manual(values=cols, breaks = names(cols)) + labs(fill='Genomic feature'))

tloc = get_legend(gheatmap(ggtree(tree),data.frame(locs[tree$tip.label])) +
                    scale_fill_manual(values=cols, breaks = names(cols)) + labs(fill='Isolate location'))

gt = ggtree(tree,layout = 'circular') + geom_treescale(x=0.0015,offset=1.1,fontsize = 3) +
  theme_nothing()

trhm = gheatmap(gt,data.frame(locs[tree$tip.label],omiss_mat,ybt0_mat),color = NA,
                colnames = F,width = 0.2) +
  scale_fill_manual(values=cols, breaks = names(cols)) + theme(legend.position = "none")
#scale_fill_brewer(palette='Set3') + labs(fill='Isolate location') 

# pdf('results/figures/heatmap/ybt0_public_tree.pdf')
p2 <- plot_grid(trhm,plot_grid(NULL,tgen,tloc,NULL,ncol=1, rel_heights = c(0.2,0.25,0.25,0.3)),
          rel_widths = c(0.8,0.2), rel_heights = c(0.9,0.1))
# dev.off()

ggsave('results/figures/heatmap/tree_ltach_kfoc_public.tiff',plot_grid(p1,p2, nrow = 2, labels = 'AUTO'), width = 10, height = 12)
