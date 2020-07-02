# Prep patient data for ml

# Load libraries
library(readxl)

# command line arguments
patient_file = snakemake@input[[1]]
keep_file = snakemake@input[[2]]
output_file = snakemake@output[[1]]

# load crkp data
metadata = read_excel(patient_file)
# change ltach names to aliases
# NOT THE SAME AS FROM THE FIRST MANUSCRIPT 
ltach_aliases = LETTERS[1:length(unique(metadata$LTACH))]
ltachs = unique(metadata$LTACH)
metadata$ltach = sapply(metadata$LTACH, function(x){
  ltach_aliases[ltachs == x]
})
# load to keep data
to_keep = read.table(keep_file,header = F)
# keep only patients we want to include in the analysis
metadata = metadata[paste0('PCMP_H',metadata$isolate_no) %in% to_keep[[1]],]
# keep columns we care about
features = c('LOSbeforeCx','age','sex','central_line','trach','foley','CHF','resp_failure','AKI','malignancy','brain_injury','gastro','obese','underweight','txp','cirrhosis','CKD','lung_dz','VDRF','decub','amikacin30','aztreonam30', 'bactrim30', 'cefepime30', 'ceftaroline30', 'ceftazidime30', 'ctx30', 'cipro30', 'colistin30', 'dapto30', 'erta30', 'flagyl30', 'gent30', 'imipenem30', 'levo30', 'linezolid30','meropenem30', 'polymyxin30', 'tigecycline30', 'tobra30', 'zosyn30', 'vanco30', 'AG30', 'poly_colistin30', 'third_gen_ceph30', 'anti_pseudo30', 'FQ30', 'carbapenem30', 'pseudoPCN30','prior_admit')
ids = paste0('PCMP_H',metadata$isolate_no)
not_keeping = names(metadata)[!names(metadata) %in% features]

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

meta_features = subset_df(metadata,features)
rownames(meta_features) = ids

meta_features = data.frame(lapply(meta_features, trimws), stringsAsFactors = FALSE)
rownames(meta_features) = ids
write.table(t(meta_features),file=output_file,sep='\t',quote=F,row.names=T)
ltachs = data.frame(metadata$ltach)
rownames(ltachs) = ids
write.table(ltachs,file='data/ltachs.tsv',sep='\t',quote=F,row.names=T)
