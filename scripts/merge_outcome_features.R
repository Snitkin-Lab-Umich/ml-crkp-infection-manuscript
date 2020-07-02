# combine features and outcome

# load libraries
suppressPackageStartupMessages(library(dplyr))

dat_paths = snakemake@input
dat = lapply(dat_paths, function(x) read.delim(x,stringsAsFactors = F))
names(dat) = gsub('data/|outcomes/|features/|patient/|genomic/|combined/|_features|.tsv','',dat_paths)
outcomes = dat[[1]]
names_outcomes = names(outcomes)
all_features = dat[2:(length(dat)-1)]
outfiles = snakemake@output

in_all = rownames(outcomes)
for(i in all_features){
  in_all = intersect(in_all,names(i))
}

outcomes = data.frame(outcomes[in_all,])
names(outcomes) = names_outcomes
rownames(outcomes) = in_all
all_features = lapply(all_features, function(x) x[,in_all])

combine_outcome_features = function(outcome_mat,feat_mat,files){
  col_outcomes = colnames(outcome_mat)
  outcome_names = rownames(outcome_mat)[order(rownames(outcome_mat))]
  outcomes = data.frame(outcome_mat[order(rownames(outcome_mat)),])
  rownames(outcomes) = outcome_names
  colnames(outcomes) = col_outcomes
  features = feat_mat[,order(names(feat_mat))]
  
  if(!all.equal(rownames(outcomes), names(features))){
    stop('Sample names of outcome and features don\'t match!')
  }
  for(outf in files){
    outcome_name = gsub('data/combined/|_.*','',outf)
    if(!outcome_name %in% colnames(outcomes)){
      print(paste('Skipping', outcome_name, 'because not in outcomes.tsv'))
      next()
    }
    outcome_features = rbind(as.character(outcomes[,outcome_name]),features)
    rownames(outcome_features)[1] = outcome_name
    write.table(outcome_features,file=outf,quote=F,sep='\t',row.names=T)
  }
}

print(names(all_features))

for(i in 1:length(all_features)){
  files = outfiles[gsub("^[^_]*_|.tsv", "", outfiles) %in% names(all_features)[i]]
  print(files)
  print(dim(all_features[[i]]))
  # if lots of data only take unique patterns (to make able to run - will have to think about interpretability later)
  if(nrow(all_features[[i]]) > 500){
    features_unique = all_features[[i]][!duplicated(all_features[[i]],MARGIN=1),]
    print(dim(features_unique))
    combine_outcome_features(outcomes,features_unique,files)
  }else{
    combine_outcome_features(outcomes,all_features[[i]],files)
  }
}


