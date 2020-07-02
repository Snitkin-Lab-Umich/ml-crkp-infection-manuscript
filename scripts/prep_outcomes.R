# prep outcome for ml

# load libraries
library(readxl)
library(dplyr)
library(tidyr)

# command line arguments
patient_file = snakemake@input[[1]]
outcome = scan(snakemake@input[[2]], character())
keep_file = snakemake@input[[3]]
outfile = snakemake@output[[1]]

# load crkp data
metadata = read_excel(patient_file)
rownames(metadata) = paste0('PCMP_H',metadata$isolate_no)
# load to keep data
to_keep = read.table(keep_file,header = F)
metadata = data.frame(metadata[rownames(metadata) %in% to_keep[[1]],])
# get outcome
bin_names = c()
print(dim(metadata))
mat_outcome = matrix(nrow=nrow(metadata))
for(out in outcome){
  print(out)
  meta_outcome = metadata[[out]]
  names(meta_outcome) = paste0('PCMP_H',metadata$isolate_no)
  # keep only patients we want to include in the analysis
  table_outcome = table(meta_outcome)
  if(length(table_outcome) > 2){
    meta_outcome_bin = sapply(names(table_outcome), function(x){
      ifelse(meta_outcome == x, 'yes','no')
    })
    rownames(meta_outcome_bin) = rownames(meta_outcome)
    bin_names = c(bin_names,names(table_outcome))
    mat_outcome = cbind(mat_outcome, meta_outcome_bin)
  }
  else{
    bin_names = c(bin_names,out)
    mat_outcome = cbind(mat_outcome, meta_outcome)
  }
}
print(bin_names)
mat_outcome = data.frame(mat_outcome[,2:ncol(mat_outcome)])
colnames(mat_outcome) = bin_names
write.table(mat_outcome,file=outfile,quote=F,sep='\t',row.names=T, col.names = T)

