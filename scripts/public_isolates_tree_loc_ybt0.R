# plot public isolates and LTACH isolates to look at icekp10 (Figure S11)

library(ape)
library(readxl)
library(ggtree)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

source('scripts/get_isolate_locations.R')


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

# get locations
# collapse ltachs
sample_key$LTACH_col = sample_key$state 

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

pdf('results/figures/heatmap/ybt0_public_tree.pdf')
plot_grid(trhm,plot_grid(NULL,tgen,tloc,NULL,ncol=1, rel_heights = c(0.2,0.25,0.25,0.3)),
          rel_widths = c(0.8,0.2), rel_heights = c(0.9,0.1))
dev.off()
