# collapse gene matrices into one gene matrix

source('/nfs/esnitkin/Zena/prewas/R/collapse_snps_by_gene.R')
source('/nfs/esnitkin/Zena/prewas/R/validate.R')

# paths
dat_paths = snakemake@input
dat = lapply(dat_paths, function(x) read.delim(x,stringsAsFactors = F))
outfile = snakemake@output[[1]]

subset_mat = function(mat,to_subset){
  mat = mat[,colnames(mat) %in% to_subset]
  mat = mat[,order(colnames(mat))]
}

# combine all matrices
in_all = intersect(colnames(dat[[1]]),colnames(dat[[2]]))
in_all = intersect(in_all,colnames(dat[[3]]))
snps = subset_mat(dat[[1]],in_all)
indels = subset_mat(dat[[2]],in_all)
insertions = subset_mat(dat[[3]],in_all)

if(all(colnames(snps) == colnames(indels) & colnames(snps) == colnames(insertions))){
  genes_uncollapsed = data.frame(rbind(snps,indels,insertions))
}else{
  'Sample names don\'t match up'
}

genes_uncollapsed = genes_uncollapsed[rowSums(genes_uncollapsed) > 1 & rowSums(genes_uncollapsed) < (ncol(genes_uncollapsed) - 1),]
genes = rownames(genes_uncollapsed)
genes_collapsed = collapse_snps_into_genes(as.matrix(genes_uncollapsed),genes)
rownames(genes_collapsed) = paste0(rownames(genes_collapsed),'_',1:nrow(genes_collapsed))
# add pangenome
genes_collapsed = rbind(genes_collapsed, dat[[4]])

write.table(genes_collapsed,outfile,sep = '\t', quote=F,row.names = T)
