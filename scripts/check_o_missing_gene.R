# statistical tests for o locus missing gene clade

library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(tidytree)
library(readxl)
library(geiger)

# read in tree
tree = midpoint.root(read.tree("../../../Sequence_data/consensus/2019_08_07_Penn_All_variant_calling/2020_02_13_16_29_09_core_results/gubbins/iqtree_results/2020_02_13_16_29_09_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile"))

mat = data.frame(t(read.delim('data/combined/infection_patient-kleborate.tsv')))
mat = mat[mat$ST == 'ST258',]
tip_order = tree$tip.label[tree$tip.label %in% rownames(mat)]
tree = midpoint.root(drop.tip(tree,tree$tip.label[!tree$tip.label %in% rownames(mat)]))
mat = mat[tip_order,]

ocs = read.delim('data/outcomes/outcomes.tsv')
src = sapply(1:nrow(ocs), function(x){
  names(ocs)[which(ocs[x,c('blood','resp','urine')] == 'yes')+1]
})
names(src) = rownames(ocs)
src = src[rownames(mat)]

# edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
# colnames(edge)=c("parent", "node", "edge_num")
# tr %<+% edge + geom_label(aes(x=branch, label=node))

# disrupted kfoC lineage: node 335
sc_kfoc = extract.clade(tree,335)

sc_kfoc_src = src[sc_kfoc$tip.label]
sc_kfoc_inf = ocs$infection[rownames(ocs) %in% sc_kfoc$tip.label]
names(sc_kfoc_inf) = rownames(ocs)[rownames(ocs) %in% sc_kfoc$tip.label]
ice_src = src[rownames(mat)[which(mat$Yersiniabactin._ICEKp10 == 1)][rownames(mat)[which(mat$Yersiniabactin._ICEKp10 == 1)] %in% sc_kfoc$tip.label]]
ice_inf = ocs$infection[rownames(ocs) %in% rownames(mat)[which(mat$Yersiniabactin._ICEKp10 == 1)][rownames(mat)[which(mat$Yersiniabactin._ICEKp10 == 1)] %in% sc_kfoc$tip.label]]

# check if disrupted kfoC lineage enriched in respiratory isolates
resp_miss = sc_kfoc_src == 'resp'
table(resp_miss)
sum(table(resp_miss))
table(resp_miss)/sum(table(resp_miss))
resp_no_miss = src[!names(src) %in% sc_kfoc$tip.label] == 'resp'
table(resp_no_miss)
sum(table(resp_no_miss))
table(resp_no_miss)/sum(table(resp_no_miss))

wilcox.test(as.numeric(resp_miss),as.numeric(resp_no_miss))


# check if icekp10 isoaltes are enriched in infection isolates (within kfoC lineage)
inf_ybt_miss = ice_inf 
table(inf_ybt_miss)
sum(table(inf_ybt_miss))
table(inf_ybt_miss)/sum(table(inf_ybt_miss))
inf_noybt_miss = sc_kfoc_inf[!names(sc_kfoc_inf) %in% names(ice_src)]
table(inf_noybt_miss)
sum(table(inf_noybt_miss))
table(inf_noybt_miss)/sum(table(inf_noybt_miss))

wilcox.test(as.numeric(inf_ybt_miss),as.numeric(inf_noybt_miss))

