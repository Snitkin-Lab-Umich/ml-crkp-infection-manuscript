# merge patient and genomic features
suppressPackageStartupMessages(library(dplyr))

print('reading in patient')
patient = read.delim(snakemake@input[[1]],stringsAsFactors=F)
print('reading in genomic')
genomic = read.delim(snakemake@input[[2]],stringsAsFactors=F)
outfile = snakemake@output[[1]]

subset_mat = function(mat,to_subset){
  mat = mat[,colnames(mat) %in% to_subset]
  mat = mat[,order(colnames(mat))]
}

in_all = intersect(colnames(patient),colnames(genomic))
patient = subset_mat(patient,in_all)
genomic = subset_mat(genomic,in_all)

if(all(colnames(patient) == colnames(genomic))){
  patient_genomic = rbind(patient,genomic)
}else{
  'Sample names don\'t match up'
}

write.table(patient_genomic,outfile,quote = F,sep = '\t', row.names = T)

