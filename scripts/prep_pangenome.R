# parse roary pangenome

pangenome_path = snakemake@input[[1]]
keep_path = snakemake@input[[2]]
outfile = snakemake@output[[1]]

# load pangenome
pangenome = read.delim(pangenome_path,row.names = 1)

# load to keep data
to_keep = read.table(keep_path,header = F)[[1]]
# keep only patients we want to include in the analysis
pangenome = pangenome[,colnames(pangenome) %in% to_keep]

# write to file
write.table(pangenome,outfile,quote = F,sep = '\t', row.names = T)
outfile2 = gsub('.tsv','_features.tsv',outfile)
write.table(pangenome,outfile2,quote = F,sep = '\t', row.names = T)
