# merge genomic features

suppressPackageStartupMessages(library(dplyr))

print('reading in snps')
snps = read.delim(snakemake@input[[1]])
print('reading in indels')
indels = read.delim(snakemake@input[[2]])
print('reading in insertions')
insertions = read.delim(snakemake@input[[3]])
print('reading in pangenome')
pangenome = read.delim(snakemake@input[[4]])
outfile = snakemake@output[[1]]

subset_mat = function(mat,to_subset){
  mat = mat[,colnames(mat) %in% to_subset]
  mat = mat[,order(colnames(mat))]
}

in_all = intersect(colnames(snps),colnames(indels))
in_all = intersect(in_all,colnames(insertions))
in_all = intersect(in_all,colnames(pangenome))
snps = subset_mat(snps,in_all)
indels = subset_mat(indels,in_all)
insertions = subset_mat(insertions,in_all)
pangenome = subset_mat(pangenome,in_all)

if(all(colnames(snps) == colnames(indels) & colnames(snps) == colnames(insertions) & colnames(snps) == colnames(pangenome))){
  genomic = data.frame(rbind(snps,indels,insertions,pangenome))
}else{
  'Sample names don\'t match up'
}

genomic = as.matrix(genomic[rowSums(genomic) > 1 & rowSums(genomic) < (ncol(genomic) - 1),])
rownames(genomic) = paste0(rownames(genomic),'_',1:nrow(genomic))
genomic[is.na(genomic)] = 0
write.table(genomic,file=outfile,sep='\t',quote=F,row.names = T)

