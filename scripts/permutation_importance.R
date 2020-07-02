# perform permutation importance

library(PRROC)
library(future.apply)

source('scripts/permutation_importance_functions.R')
source('ML_pipeline_microbiome/code/R/calc_aucs.R')

in_file = snakemake@input[[1]]
dat = snakemake@wildcards$dat 
outfile = snakemake@output[[1]]
corr_thresh = as.numeric(snakemake@wildcards$corr) 
print(corr_thresh)
p_thresh = as.numeric(snakemake@wildcards$pval) 
print(p_thresh)

  in_file = paste0('data/combined/',dat,'.tsv')
  corr_file = paste0('data/correlation_matrix/',dat,'_',corr_thresh,'_',p_thresh,'.tsv')
  # Create correlation matrix of machine learning data
  #   filters correlation >= cor_value and p values < p_value
  #   default values are cor_value = 1, and p_value = 0.1
  compute_correlation_matrix(input_file = in_file, 
                              outcome = 'infection', 
                              output_file = corr_file, cor_value = corr_thresh, p_value = p_thresh)
  
  files = list.files(paste0('results/',dat,'_lr'),pattern='.RData',recursive=T, full.names = T)
  
  load(files[1])
  pref = gsub('.*/|.RData','',files[1])
  test_data = results$testData
  
  corr = read_tsv(corr_file)
  
  all_feats = colnames(test_data)[2:ncol(test_data)]
  corr_feats = unique(c(corr$column,corr$row))
  noncorr_feats = all_feats[!all_feats %in% corr_feats]
  
  grps = as.list(noncorr_feats)
  accounted_for = c()
  c=length(grps)+1
  for(i in corr_feats){
    if(i %in% accounted_for) next
    feats = unique(c(i,corr$row[corr$column == i],corr$column[corr$row == i]))
    new_feats = T
    while(new_feats){
      len_feats = length(feats)
      for(j in feats){
        feats = unique(c(feats,j,corr$row[corr$column == j],corr$column[corr$row == j]))
      }
      new_feats = length(feats) > len_feats
    }
    grps[[c]] = feats
    accounted_for = c(accounted_for,feats)
    c = c+1
  }
  
  nums = data.frame(matrix(NA,nrow=length(files),ncol=length(grps)))
  colnames(nums) = sapply(grps,paste,collapse='|')
  rownames(nums) = gsub('.*/|.RData','',files)
  
  for(i in files){
    
    load(i)
    pref = gsub('.*/|.RData','',i)
    
    trained_model = results$trained_model
    test_data = results$testData
    outcome = 'infection'
    test_auc = results$test_auc
    
    imps <- do.call('rbind', future_lapply(colnames(nums), function(feat){
      auc_diffs <- future_sapply(0:99, function(s){
        set.seed(s)
        full_permuted <- test_data
        # PERMUTE GROUPED TOGETHER!! 
        fs = strsplit(feat, '\\|')[[1]]
        if(length(fs) == 1) full_permuted[,fs] <- sample(test_data[,fs])
        else full_permuted[,fs] <- t(sample(data.frame(t(test_data[,fs])),))
        # Predict the diagnosis outcome with the one-feature-permuted test dataset
        rpartProbs_permuted <- predict(trained_model, full_permuted, type="prob")
        # Calculate the new auc
        new_auc <- calc_aucs(rpartProbs_permuted, full_permuted[,outcome])$auroc
        # Return how does this feature being permuted effect the auc
        return(test_auc - new_auc)
      })
      auc_diff = mean(auc_diffs)
      }))
    # Save non correlated results in a dataframe.
    imps <- as.data.frame(imps) %>%
      mutate(names=factor(colnames(nums))) %>%
      rename(mean_auc_diff=V1)

    nums[pref,] = imps$mean_auc_diff
  }
  
  write_tsv(nums,outfile)
  
  nums %>% pivot_longer(cols = everything()) %>% ggplot(aes(x=reorder(name, value, mean),y=value)) +
    geom_boxplot() + geom_hline(yintercept = 0, col='red') + 
    theme(axis.text.x = element_text(angle = 20, hjust=1))+
    ylab('Mean difference between test AUC and permuted AUC') + xlab('') 
  
  ggsave(paste0('results/permutation_importance/figures/',dat,'_',corr_thresh,'_',p_thresh,'.pdf'))
  
