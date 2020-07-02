# Subset isolates so only the first one from each patient is kept

# Load libraries
library(readxl)

# command line arguments
patient_file = snakemake@input[[1]]
kleb_file = snakemake@input[[2]]
output_file = snakemake@output[[1]]

# load crkp data
metadata = read_excel(patient_file)
# load kleborate data
kleb = read.delim(kleb_file)
# sort by Cx_date
metadata = metadata[order(metadata$Cx_date),]
# keep only first patient sample
metadata = metadata[!duplicated(metadata$Patient_ID),]
# remove wound isolates
metadata = metadata[!metadata$source == 'wound',]
# keep only ST258
metadata$strain = paste0('PCMP_H',metadata$isolate_no)
metadata$ST = sapply(metadata$strain, function(x){
  st = as.character(kleb$ST[kleb$strain == x])
  if(length(st) == 0){
    return('ST?')
  }else{
    return(st)
  }
})
metadata_st258 = metadata[metadata$ST == 'ST258',]

to_keep = paste0('PCMP_H',metadata$isolate_no)
to_keep_st258 = paste0('PCMP_H',metadata_st258$isolate_no)

write.table(to_keep_st258,file=output_file,quote=F,sep='\t',row.names=F, col.names = F)
