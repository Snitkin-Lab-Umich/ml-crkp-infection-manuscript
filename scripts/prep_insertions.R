# parse panisa output into presence/absence matrix by variant and by gene

source('/nfs/esnitkin/bin_group/pipeline/Github/scripts/panISa_parser_functions.R')
source('/nfs/esnitkin/Zena/prewas/R/collapse_snps_by_gene.R')
source('/nfs/esnitkin/Zena/prewas/R/validate.R')

# paths
isfinder_path = snakemake@input[[1]]
gb_path = snakemake@input[[2]]
keep_path = snakemake@input[[3]]
outfile = snakemake@output[[1]]
outfile_genes = snakemake@output[[2]]

# parse IS elements info into matrix
ismat = parse_is(isfinder_path,gb_path)

mat = ismat$mat
genes = sapply(strsplit(rownames(ismat$mat),'\\|'),function(x) x[2])
gene_mat = collapse_snps_into_genes(mat,genes)

# load to keep data
to_keep = read.table(keep_path,header = F)[[1]]
# keep only patients we want to include in the analysis
mat = ismat$mat[,colnames(ismat$mat) %in% to_keep]
gene_mat = gene_mat[,colnames(gene_mat) %in% to_keep]

# write variant matrix file
write.table(mat,outfile,quote = F,sep = '\t', row.names = T)

# write gene matrix file
write.table(gene_mat,outfile_genes,quote = F,sep = '\t', row.names = T)

