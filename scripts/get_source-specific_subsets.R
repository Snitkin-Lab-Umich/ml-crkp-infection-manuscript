# subset isolates to each source (resp and urine) to perform source-specific analysis

# Load libraries
library(readxl)

# command line arguments
patient_file = snakemake@input[[1]]
combined_files = snakemake@input[2:length(snakemake@input)]
output_files = snakemake@output

# read in patient data
patient = read_excel(patient_file)

for(i in 1:length(output_files)){
  
  output_file = output_files[[i]]
  
  # get source
  source = gsub('_.*','',output_file)
  source = gsub('data/combined/.*-','',source)
  
  combined_file = gsub(paste0('-',source),'',output_file)
  
  # read in data to be used for ml
  combined = read.delim(combined_file,stringsAsFactors=F)
  
  # subset to only include source of interest
  combined_subset = combined[,names(combined) %in% paste0('PCMP_H',patient$isolate_no)[patient$source == source]]
  combined_subset = combined_subset[apply(combined_subset,1,function(x) length(unique(x))) != 1,]
  
  print(paste('writing',output_file))
  # write subsetted data to file
  write.table(combined_subset,file=output_file,sep='\t',quote=F,row.names=T)
}


