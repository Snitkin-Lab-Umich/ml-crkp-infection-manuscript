ML Results
================
2021-01-06

Machine learning algorithms used include: glmnet, rf, svmRadial. Models
were trained with 100 different random partitions of the data into
training and testing sets using 5-fold cross validation.

## Model Performance

<img src="figures/performance.png" width="80%" />

## Hyperparameter Performance

<img src="figures/hp_performance_glmnet.png" width="80%" /><img src="figures/hp_performance_rf.png" width="80%" /><img src="figures/hp_performance_svmRadial.png" width="80%" />

## Memory Usage & Runtime

<img src="figures/benchmarks.png" width="80%" />

Each model training run was given 36 cores for parallelization.
