
# Author: Begum Topcuoglu
# Date: 2019-01-14
######################################################################
# Description:
# This script trains and tests the model according to proper pipeline
######################################################################

######################################################################
# Dependencies and Outputs:
#    Model to put to function:
#       1. "L2_Logistic_Regression"
#       2. "L2_Linear_SVM"
#       3. "RBF_SVM"
#       4. "Decision_Tree"
#       5. "Random_Forest"
#       6. "XGBoost"
#    data to put to function:
#         Features: Hemoglobin levels and 16S rRNA gene sequences in the stool
#         Labels: - Colorectal lesions of 490 patients.
#                 - Defined as cancer or not.(Cancer here means: SRN)
#
# Usage:
# Call as source when using the function. The function is:
#   pipeline(data, model)

# Output:
#  A results list of:
#     1. cvAUC and testAUC for 1 data-split
#     2. cvAUC for all hyper-parameters during tuning for 1 datasplit
#     3. feature importance info on first 10 features for 1 datasplit
#     4. trained model as a caret object
######################################################################

######################################################################
#------------------------- DEFINE FUNCTION -------------------#
######################################################################
source("ML_pipeline_microbiome/code/R/tuning_grid.R")
source("ML_pipeline_microbiome/code/R/permutation_importance.R")
source("ML_pipeline_microbiome/code/R/calc_aucs.R")
source('ML_pipeline_microbiome/code/R/split_by_group.R')

model_pipeline <- function(data, model, split_number, outcome=NA, hyperparameters=NA, level=NA, permutation=TRUE,group=NULL){

    # group is a vector of group labels if the user wants to split by group when splitting the train and test sets and when performing cross-validation

  # -----------------------Get outcome variable----------------------------->
  # If no outcome specified, use first column in data
  if(is.na(outcome)){
    outcome <- colnames(data)[1]
  }else{
    # check to see if outcome is in column names of data
    if(!outcome %in% colnames(data)){
      stop(paste('Outcome',outcome,'not in column names of data.'))
    }

		# Let's make sure that the first column in the data frame is the outcome variable
		temp_data <- data.frame(outcome = data[,outcome])
		colnames(temp_data) <- outcome
		data <- cbind(temp_data, data[, !(colnames(data) %in% outcome)]) # want the outcome column to appear first
  }

  # ------------------Check data for pre-processing------------------------->
  # Data is pre-processed in code/R/setup_model_data.R
  # This removes OTUs with near zero variance and scales 0-1
  # Then generates a correlation matrix
  # Test if data has been preprocessed - range 0-1 and are not all 0s
  feature_summary <- any(c(min(data[,-1]) < 0, 
    max(data[,-1]) > 1, 
    any(apply(data[,-1], 2, sum) == 0)))
  if(feature_summary){
    stop('Data has not been preprocessed, please use "code/R/setup_model_data.R" to preprocess data')
  }

  # ------------------Randomize features----------------------------------->
  # Randomize feature order, to eliminate any position-dependent effects 
  features <- sample(colnames(data[,-1]))
  data <- select(data, one_of(outcome), one_of(features))


  # ----------------------------------------------------------------------->
  # Get outcome variables
  first_outcome = as.character(data[1,outcome])
  outcome_vals = unique(data[,outcome])
  if(length(outcome_vals) != 2) stop('A binary outcome variable is required.')
  second_outcome = as.character(outcome_vals[!outcome_vals == first_outcome])
  print(paste(c('first outcome:','second outcome:'),c(first_outcome,second_outcome)))


  # ------------------80-20 Datasplit for each seed------------------------->
  # Do the 80-20 data-split
  # Stratified data partitioning %80 training - %20 testing
  if(is.null(group)){
    inTraining <- createDataPartition(data[,outcome], p = .80, list = FALSE)
  }else{
      # ensure groups are only in train set or only in test set (while maintaining close to an 80-20 split -- maximum 80% of data in train set)
    inTraining <- createGroupedDataPartition(group, p = 0.80)
    train_group <- group[inTraining]
    test_group <- group[-inTraining]
  }
  train_data <- data[ inTraining,]
  test_data  <- data[-inTraining,]

  # ----------------------------------------------------------------------->

  # -------------Define hyper-parameter and cv settings-------------------->
  # Define hyper-parameter tuning grid and the training method
  # Uses function tuning_grid() in file ('code/learning/tuning_grid.R')
  if(is.null(group)){
    tune <- tuning_grid(train_data, model, outcome, hyperparameters)
  }else{
    tune <- tuning_grid(trainTransformed, model, outcome, hyperparameters, group=train_group)
  }
  grid <- tune[[1]]
  method <- tune[[2]]
  cv <- tune[[3]]
  # ----------------------------------------------------------------------->

  # ---------------------------Train the model ---------------------------->
  # ------------------------------- 1. -------------------------------------
  # - We train on the 80% of the full data.
  # - We use the cross-validation and hyper-parameter settings defined above to train
  # ------------------------------- 2. -------------------------------------
  # We use ROC metric for all the models
  # To do that I had to make changes to the caret package functions.
  # The files 'data/caret_models/svmLinear3.R and svmLinear5.R are my functions.
  # I added 1 line to get Decision Values for linear SVMs:
  #
  #           prob = function(modelFit, newdata, submodels = NULL){
  #             predict(modelFit, newdata, decisionValues = TRUE)$decisionValues
  #           },
  #
  # This line gives decision values instead of probabilities and computes ROC in:
  #   1. train function with the cross-validataion
  #   2. final trained model
  # using decision values and saves them in the variable "prob"
  # ------------------------------- 3. --------------------------------------
  # - If the model is logistic regression, we need to add a family=binomial parameter.
  # - If the model is random forest, we need to add a ntree=1000 parameter.
  #         We chose ntree=1000 empirically.
  # ----------------------------------------------------------------------->
  # Make formula based on outcome
  f <- as.formula(paste(outcome, '~ .'))
  print('Machine learning formula:')
  print(f)
  # Start walltime for training model
  tic("train")
  if(model=="L2_Logistic_Regression"){
  print(model)
  trained_model <-  train(f, # label
                          data=train_data, #total data
                          method = method,
                          trControl = cv,
                          metric = "ROC",
                          tuneGrid = grid,
                          family = "binomial")
  }
  else if(model=="Random_Forest"){
      print(model)
      trained_model <-  train(f,
                              data=train_data,
                              method = method,
                              trControl = cv,
                              metric = "ROC",
                              tuneGrid = grid,
                              ntree=1000) # not tuning ntree
  }
  else{
    print(model)
    trained_model <-  train(f,
                            data=train_data,
                            method = method,
                            trControl = cv,
                            metric = "ROC",
                            tuneGrid = grid)
  }
  # Stop walltime for running model
  seconds <- toc()
  # Save elapsed time
  train_time <- seconds$toc-seconds$tic
  # Save wall-time
  write.csv(train_time, file=paste0("data/temp/", level, "/traintime_", model, "_", split_number, ".csv"), row.names=F)
  # ------------- Output the cvAUC and testAUC for 1 datasplit ---------------------->
  # Mean cv AUC value over repeats of the best cost parameter during training
  cv_auc <- getTrainPerf(trained_model)$TrainROC
  # Save all results of hyper-parameters and their corresponding meanAUCs over 100 internal repeats
  results_individual <- trained_model$results
  # ---------------------------------------------------------------------------------->

  # -------------------------- Feature importances ----------------------------------->
  #   if linear: Output the weights of features of linear models
  #   else: Output the feature importances based on random permutation for non-linear models
  # Here we look at the top 20 important features
  if(permutation){
    if(model=="L1_Linear_SVM" || model=="L2_Linear_SVM" || model=="L2_Logistic_Regression"){
      # We will use the permutation_importance function here to:
      #     1. Predict held-out test-data
      #     2. Calculate ROC and AUROC values on this prediction
      #     3. Get the feature importances for correlated and uncorrelated feautures
      roc_results <- permutation_importance(trained_model, test_data, first_outcome, second_outcome, outcome, level)
      test_auc <- roc_results[[1]]  # Predict the base test importance
      feature_importance_non_cor <- roc_results[2] # save permutation results
      # Get feature weights
      feature_importance_cor <- trained_model$finalModel$W
      auprc <-roc_results[[6]]
      sensitivity <-roc_results[[4]]
      specificity <-roc_results[[5]]
    }
    else{
      # We will use the permutation_importance function here to:
      #     1. Predict held-out test-data
      #     2. Calculate ROC and AUROC values on this prediction
      #     3. Get the feature importances for correlated and uncorrelated feautures
      roc_results <- permutation_importance(trained_model, test_data, first_outcome, second_outcome, outcome, level)
      test_auc <- roc_results[[1]] # Predict the base test importance
      feature_importance_non_cor <- roc_results[2] # save permutation results of non-cor
      feature_importance_cor <- roc_results[3] # save permutation results of cor
      auprc <-roc_results[[6]]
      sensitivity <-roc_results[[4]]
      specificity <-roc_results[[5]]
    }
  }else{
    print("No permutation test being performed.")
    if(model=="L1_Linear_SVM" || model=="L2_Linear_SVM" || model=="L2_Logistic_Regression"){
      # Get feature weights
      feature_importance_non_cor <- trained_model$finalModel$W
      # Get feature weights
      feature_importance_cor <- trained_model$finalModel$W
    }else{
      # Get feature weights
      feature_importance_non_cor <- NULL
      # Get feature weights
      feature_importance_cor <- NULL
    }
    # Calculate the test-auc for the actual pre-processed held-out data
    rpartProbs <- predict(trained_model, test_data, type="prob")
    aucs <- calc_aucs(rpartProbs, test_data[,outcome])
    test_auc <- aucs$auroc
    auprc <- aucs$auprc
    bin_outcome <- ifelse(test_data[,outcome] == names(rpartProbs)[1], 1, 0)
    test_roc <- roc(bin_outcome,rpartProbs[[1]],direction='<')
    #roc(ifelse(test_data[,outcome] == first_outcome, 1, 0), rpartProbs[[1]])
    # Calculate the test auprc (area under precision-recall curve)
    #bin_outcome <- get_binary_outcome(test_data[,outcome], first_outcome)
    #auprc <- calc_auprc(rpartProbs[[1]], bin_outcome)
    # Calculate sensitivity and specificity for 0.5 decision threshold.
    p_class <- ifelse(rpartProbs[[1]] > 0.5, second_outcome, first_outcome)
    r <- confusionMatrix(as.factor(p_class), test_data[,outcome])
    sensitivity <- r$byClass[[1]]
    specificity <- r$byClass[[2]]
    # best decision threshold (youden method)
    thr <- coords(test_roc, "best", ret = "threshold")
  }

  # ---------------------------------------------------------------------------------->

  # ----------------------------Save metrics as vector ------------------------------->
  # Return all the metrics
  #results <- list(cv_auc, test_auc, results_individual, feature_importance_non_cor, feature_importance_cor, trained_model, sensitivity, specificity, auprc)
    # Return all the metrics
  results <- list(cv_auc=cv_auc, test_auc=test_auc, results_individual=results_individual, feature_importance_non_cor=feature_importance_non_cor, feature_importance_cor=feature_importance_cor, trained_model=trained_model, sensitivity=sensitivity, specificity=specificity,auprc=auprc,trainData=train_data,testData=test_data)
  return(results)
}
