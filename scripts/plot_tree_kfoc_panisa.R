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

thm = gheatmap(tr_ann,dat,
               color = NA, offset = 0,width = 0.4, colnames = F) + 
  scale_fill_manual(values=cols_mat,breaks=names(cols_mat),name='') + theme(legend.position = 'none')

pdf('results/figures/heatmap/tree_kleb_panisa.pdf')#, height = 6, width = 6)
plot_grid(thm,plot_grid(NULL,hmo_leg,hmk_leg,hmb_leg,hmp_leg,NULL,
                        ncol=1,align = 'hv',rel_heights = c(0.2,0.15,0.15,0.15,0.15,0.2)),
          nrow=1,rel_widths = c(0.7,0.3))
dev.off()

table(dat,useNA='always')

samps = rownames(dat[is.na(dat$V2) & !is.na(dat$potential_is_all),])
kleb = read.delim('../2019-05-09_kleborate/crkp_kleborate_results.txt')
View(kleb[rownames(kleb) %in% samps,])

