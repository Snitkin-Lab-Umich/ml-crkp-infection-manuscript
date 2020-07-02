# read and prep kleborate results for ml

# command line arguments
kleborate_file = snakemake@input[[1]] 
keep_file = snakemake@input[[2]] 
output_file_vars = snakemake@output[[1]] 
output_file_genes = snakemake@output[[2]] 

kleborate = read.delim(kleborate_file, stringsAsFactors = F)

# load to keep data
to_keep = read.table(keep_file,header = F)
# keep only patients we want to include in the analysis
kleborate = kleborate[kleborate$strain %in% to_keep[[1]],]
# keep columns we care about
features = c('ST','virulence_score','resistance_score','num_resistance_classes','num_resistance_genes','Yersiniabactin','Colibactin','Aerobactin','Salmochelin','rmpA','rmpA2','wzi','K_locus','K_locus_missing_genes','O_locus','O_locus_missing_genes','AGly','Col','Fcyn','Flq','Gly','MLS','Ntmdz','Phe','Rif','Sul','Tet','Tmt','Bla','Bla_Carb','Bla_ESBL','Bla_ESBL_inhR','Bla_broad','Bla_broad_inhR')
ids = kleborate$strain
not_keeping = names(kleborate)[!names(kleborate) %in% features]

# function to subset df and make all categorical variables factors
subset_df = function(metadata,cnames){
  meta = metadata[,names(metadata) %in% cnames]
  # change to data frame
  meta = as.data.frame(meta)
  # change categorical variables to factors
  meta = data.frame(lapply(meta, function(x){
    if(is.character(x)){
      factor(x)
    }else{
      x
    }
  }))
}

kleborate = kleborate[,names(kleborate) %in% features]
rownames(kleborate) = ids
kleborate = kleborate[,!sapply(kleborate,function(x) all(x == '-'))]
# remove stars
kleborate = data.frame(sapply(kleborate,function(x) gsub('\\*|\\?','',x)))

# genes
kleb_genes = kleborate
kleb_genes[kleb_genes == '-' | kleb_genes == ''] = 0
kleb_genes[kleb_genes != 0] = 1
kleb_genes[,1:5] = kleborate[,1:5]

# variants
to_split = sapply(names(kleborate), function(x) grepl(';|,',paste(kleborate[,x],collapse='')))

kleb_list = sapply(names(kleborate), function(feat){
  if(!to_split[feat]){
    df = data.frame(kleborate[,feat])
    names(df) = feat
    return(df)
  }else{
    # if similar according to kleborate (* or ?), call it the same gene
    vars = unique(gsub('\\*|\\?','',unlist(strsplit(as.character(kleborate[,feat]),';|,'))))
    splitdf = data.frame(matrix(0,nrow=nrow(kleborate),ncol=length(vars)))
    names(splitdf) = vars
    rownames(splitdf) = rownames(kleborate)
    all_vars = sapply(1:nrow(splitdf), function(samp){
        unique(gsub('\\*|\\?','',unlist(strsplit(as.character(kleborate[samp,feat]),';|,'))))
    })
    for(i in 1:length(all_vars)){
      splitdf[i,all_vars[[i]]] = 1
    }
    splitdf = splitdf[,names(splitdf) != '-']
    names(splitdf) = gsub(' ','_',names(splitdf))#paste0(feat,'_',names(splitdf)))
    return(splitdf)
  }
})


kleb_split = do.call(cbind,kleb_list)
# any mgrb mutation
kleb_split$Col = gsub('MgrB-.*','MgrB-trunc',kleb_split$Col)

kleb = subset_df(kleb_split,names(kleb_split))
rownames(kleb) = ids

kleb = t(kleb)
kleb_genes = t(kleb_genes)

write.table(kleb,file=output_file_vars,sep='\t',quote=F,row.names=T)
write.table(kleb_genes,file=output_file_genes,sep='\t',quote=F,row.names=T)

to_remove = c('ST','virulence_score','resistance_score','num_resistance_classes','num_resistance_genes') 

kleb_subset = kleb[!rownames(kleb) %in% to_remove,]
kleb_genes_subset = kleb_genes[!rownames(kleb_genes) %in% to_remove,]

write.table(kleb_subset,file=gsub('_features.tsv','-subset_features.tsv',output_file_vars),sep='\t',quote=F,row.names=T)
write.table(kleb_genes_subset,file=gsub('_features.tsv','-subset_features.tsv',output_file_genes),sep='\t',quote=F,row.names=T)

