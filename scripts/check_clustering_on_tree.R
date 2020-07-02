# check clustering of outcome and source on tree

# load libraries ----
library(ape)
library(phytools)
library(future.apply)

# functions ----
revert_list_str <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- future_lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  future_apply(do.call(rbind, x), 2, as.list) 
}
#largest_subtree - Takes as input: 
# 1) subtrees - A list produced by subtrees, that includes all isolates of interest
# 2) isolate_labels - A named vactor of labels by which pure clusters are defined
#Returns lists containing the largest pure subtree that each isolate belongs to and the index of that subtree
largest_subtree <- function(subtrees, isolate_labels){
  largest_st_info = future_lapply(names(isolate_labels), function(i){
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. 
    #CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE EPI LABEL AND 2) HAVE BOOTSTRAP SUPPORT GREATER THAN 90 
    largest_st = max(future_sapply(subtrees, FUN = function(st){ 
      if(sum(grepl(i, st$tip.label, perl = TRUE)) > 0 && 
         length(unique(isolate_labels[intersect(st$tip.label, names(isolate_labels))])) == 1 && 
         !is.na(as.numeric(st$node.label[[1]])) &&  as.numeric(st$node.label[[1]]) > 90){
        length(intersect(names(isolate_labels), st$tip.label))
      }else{
        0
      }
    }))
    #GET THE INDEX OF THE LARGEST SUBTREE
    largest_st_i = which.max(future_sapply(subtrees, FUN = function(st){
      if(sum(grepl(i, st$tip.label, perl = TRUE)) > 0 && 
         length(unique(isolate_labels[intersect(st$tip.label, names(isolate_labels))])) == 1 &&
         !is.na(as.numeric(st$node.label[[1]])) && as.numeric(st$node.label[[1]]) > 90){
        length(intersect(names(isolate_labels), st$tip.label))
      }else{
        0
      }
    }))
    #GET EDGES BELONGING TO SUBTREES
    largest_st_edges = which.edge(subtrees[[1]], subtrees[[largest_st_i]]$tip.label)
    return(list(largest_st=largest_st,
                largest_st_i=largest_st_i,
                largest_st_edges=largest_st_edges))
  })#end for
  largest_st_info = revert_list_str(largest_st_info)
  return(largest_st_info)
}#end largest_subtree

#largest_subtree_loc - controls for location
# Takes as input: 
# 1) subtrees - A list produced by subtrees, that includes all isolates of interest
# 2) isolate_labels - A named vactor of labels by which pure clusters are defined
#Returns lists containing the largest pure subtree that each isolate belongs to and the index of that subtree
largest_subtree_loc <- function(subtrees, isolate_labels, isolate_locations){
  largest_st_info = future_lapply(names(isolate_labels), function(i){
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. 
    #CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE EPI LABEL, 2) INCLUDE MORE THAN ONE LTACH, AND 3) HAVE BOOTSTRAP SUPPORT GREATER THAN 90 
    sts = future_sapply(subtrees, FUN = function(st){ 
      if(sum(grepl(i, st$tip.label, perl = TRUE)) > 0 && # isolate is in subtree
         length(unique(isolate_labels[intersect(st$tip.label, names(isolate_labels))])) == 1 && # only one label in subtree
         length(unique(isolate_locations[intersect(st$tip.label, names(isolate_locations))])) > 1 && # more than one location in subtree
         !is.na(as.numeric(st$node.label[[1]])) &&  as.numeric(st$node.label[[1]]) > 90){ # bootstrap support > 90
        length(intersect(names(isolate_labels), st$tip.label))
      }else{
        0
      }
    })
    # GET THE LARGEST SUBTREE
    largest_st = max(sts)
    # GET THE INDEX OF THE LARGEST SUBTREE
    largest_st_i = which.max(sts)
    # GET EDGES BELONGING TO SUBTREES
    largest_st_edges = which.edge(subtrees[[1]], subtrees[[largest_st_i]]$tip.label)
    return(list(largest_st=largest_st,
                largest_st_i=largest_st_i,
                largest_st_edges=largest_st_edges))
  })#end for
  largest_st_info = revert_list_str(largest_st_info)
  return(largest_st_info)
}#end largest_subtree_loc

largest_subtr_perm = function(subtrs,label,labelname, loc=NULL){
  # find largest subtrees
  if(is.null(loc)){
    largest_subtr <- largest_subtree(subtrs, label)
  }else{
    largest_subtr <- largest_subtree_loc(subtrs, label, loc)
  }
  #GET CLUSTERS FOR RANDOMIZED LABELS
  in_cluster = numeric(1001)
  in_cluster[1] = sum(largest_subtr[[2]] > 1 & largest_subtr[[1]] > 1)
  rand_labels = label
  in_cluster_perm = future_sapply(1:1000, function(r){
    names(rand_labels) = sample(names(label))
      if(is.null(loc)){
        largest_subtr_rand_list <- largest_subtree(subtrs, rand_labels)
      }else{
        largest_subtr_rand_list <- largest_subtree_loc(subtrs, rand_labels, loc)
      }
    sum(largest_st_rand_list[[2]] > 1 & largest_st_rand_list[[1]] > 1)    
  })
  
  in_cluster[2:1001] = in_cluster_perm
  
  p = (1 + sum(in_cluster[2: length(in_cluster)] >= in_cluster[1]))/ (1 + length(in_cluster)-1);
  
  file = paste("results/figures/cluster/", labelname,'_cluster.pdf', sep = "")
  pdf(file, height = 8 , width = 8)
  h = hist(in_cluster[2:length(in_cluster)], 20, xlim = c(0, max(in_cluster + 2)), main = paste("Null distribution of number of isolates in ", 'pure ', labelname,  ' clusters ', '\n', "(p = ", round(p,3), ")", sep = ""), xlab = "Number of isolates in a pure cluster", col = "lightgray")
  par(new = TRUE)
  points(in_cluster[1],0, col = "black", bg = "red", cex = 2, axis = FALSE, pch = 23, xlim = c(0, max(in_cluster + 2)), ylim = c(0,max(h$counts) + 5))
  dev.off()
  return(in_cluster)
}

# read in & clean data ----
# read in tree
tree = midpoint.root(read.tree("../../../Sequence_data/consensus/2019_08_07_Penn_All_variant_calling/2019_09_05_16_29_09_core_results/gubbins/iqtree_results/2019_09_05_16_29_09_KPNIH1_genome_aln_w_alt_allele_unmapped.treefile"))

# read in outcomes
ocs = read.delim('data/outcomes/outcomes.tsv')
ltachs = read.delim('data/ltachs.tsv')

# subset tree to only ones used in ml
dat = data.frame(t(read.delim('data/features/genomic/kleborate_features.tsv')))
tree = drop.tip(tree,tree$tip.label[!tree$tip.label %in% rownames(dat)])

# order ocs correctly
ocs = ocs[tree$tip.label,]
ltachs = ltachs[tree$tip.label,]
names(ltachs) = tree$tip.label

# get outcomes ----
infection = ocs$infection
names(infection) = rownames(ocs)
blood = ocs$blood
names(blood) = rownames(ocs)
resp = ocs$resp
names(resp) = rownames(ocs)
urine = ocs$urine
names(urine) = rownames(ocs)

# get subtrees ----
subtrs = subtrees(tree)

# run permutation test ----
in_cluster_inf = largest_subtr_perm(subtrs,infection,'infection_st258',ltachs)
in_cluster_blood = largest_subtr_perm(subtrs,blood,'blood_st258',ltachs)
in_cluster_resp = largest_subtr_perm(subtrs,resp,'resp_st258',ltachs)
in_cluster_urine = largest_subtr_perm(subtrs,urine,'urine_st258',ltachs)
#in_cluster_ltach = largest_subtr_perm(subtrs,ltachs,'ltach_st258')

save.image('cluster_test_ltach_st258.RData')





