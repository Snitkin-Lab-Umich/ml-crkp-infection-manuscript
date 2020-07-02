# parse snp matrix into presence/absence matrix by variant and by gene

library(ape)
library(tidyverse)

source('/nfs/esnitkin/Zena/prewas/R/collapse_snps_by_gene.R')
source('/nfs/esnitkin/Zena/prewas/R/validate.R')
source('/nfs/esnitkin/Zena/snitkitr/R/parse_snps.R')
source('/nfs/esnitkin/Zena/snitkitr/R/validate.R')
source('/nfs/esnitkin/Zena/snitkitr/R/parse.R')
source('/nfs/esnitkin/Zena/snitkitr/R/reference_alleles.R')

# paths
snpcode_path = snakemake@input[[1]] 
snpallele_path = snakemake@input[[2]]
tree_path =  snakemake@input[[3]] 
keep_path = snakemake@input[[4]] 
outfile = snakemake@output[[1]]
outfile_genes = snakemake@output[[2]]

# read in tree
tree = read.tree(tree_path)
tree = drop.tip(tree,"gi_661922017_gb_CP008827.1_")

# read in matrices
snp_mat = read.delim(snpcode_path,row.names = 1)
colnames(snp_mat) = gsub('_R1.fastq.gz','',colnames(snp_mat))
allele_mat = read.delim(snpallele_path,row.names = 1)
colnames(allele_mat) = gsub('_R1.fastq.gz','',colnames(allele_mat))

# parse snp matrix into matrix
parsed = parse_snps(snp_mat,allele_mat,tree,return_binary_matrix = T)

mat = parsed$bin$mat
mat[is.na(mat)] = 0
genes = parsed$bin$annots$locus_tag
snpeff_impact = parsed$bin$annots$snpeff_impact
rownames(mat) = paste0(parsed$bin$annots$raw_rownames,'_',1:nrow(mat))
gene_mat = collapse_snps_into_genes(as.matrix(mat[snpeff_impact != 'LOW',]),as.character(genes[snpeff_impact != 'LOW']))

# load to keep data
to_keep = read.table(keep_path,header = F)[[1]]
# keep only patients we want to include in the analysis
mat = mat[,colnames(mat) %in% to_keep]
gene_mat = gene_mat[,colnames(gene_mat) %in% to_keep]

# write variant matrix file
write.table(mat,outfile,quote = F,sep = '\t', row.names = T)

# write gene matrix file
write.table(gene_mat,outfile_genes,quote = F,sep = '\t', row.names = T)

