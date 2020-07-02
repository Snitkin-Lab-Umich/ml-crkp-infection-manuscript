# split by group (e.g. facility)
# try to get ~80-20 split

# split into train and test set while splitting by group
# group: a vector of groups whose length matches the number of rows in the overall data set.
# p: the maximum percentage of data that goes to training (maybe less depending on group sizes)
# returns row position integers corresponding to the training data.
createGroupedDataPartition = function(group, p) {
  # get indices
  indices = seq(along = group)
  # get unique groups
  grps = unlist(unique(group))
  # initialize train groups & set
  train_grps = grps
  train_set = indices
  train_set_grp = group
  # initialize fraction of data in train set
  frac_in_train = length(train_set)/length(indices)
  # keep removing data from train set until fraction in train set < 0.8
  while(frac_in_train > p){
    # randomly choose a group
    grp = sample(train_grps, size=1)
    # remove group from train groups & set
    train_grps = train_grps[train_grps != grp]
    train_set = train_set[train_set_grp != grp]
    train_set_grp = train_set_grp[train_set_grp != grp]
    # calcuate fraction of data in train set
    frac_in_train = length(train_set)/length(indices)
  }
  # get test group & set
  #test_set = indices[!(indices %in% train_set)]
  return(train_set)
} 

# like createMultiFolds but still splitting by group using groupKFold
# for splitting into folds for cross-validation
groupKMultiFolds <- function (y, k = 10, times = 5) 
{
  if (class(y)[1] == "Surv") 
    y <- y[, "time"]
  prettyNums <- paste("Rep", gsub(" ", "0", format(1:times)), 
                      sep = "")
  for (i in 1:times) {
    tmp <- groupKFold(y, k = k)
    names(tmp) <- paste("Fold", gsub(" ", "0", format(seq(along = tmp))), 
                        ".", prettyNums[i], sep = "")
    out <- if (i == 1) 
      tmp
    else c(out, tmp)
  }
  out
}
