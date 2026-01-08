# LIRT-MI
Adjusting for longitudinal measurement noninvariance in longitudinal item response theory models: Comparing Bayesian regularization and alignment optimization

By accounting for longitudinal measurement noninvariance in longitudinal item response theory models, we proposed a longitudinal IRT model with time-specific item parameters and additionally incorporated the dependence among the repeated item measurements. Based on the principle of approximate invariance, Bayesian regularization and alignment optimization methods are developed for model estimation.

The Data file folder includes three files. Data-AwC.R implements our alignment optimization procedure for the real data. Data-BR.R implements our Bayesian adaptive Lasso procedure for the real data. fy11.csv is the final data file.

The NU500L1LM file folder includes four files for the condition of non-uniform DIF, I=500, Small DIF, 1/8 DIF and alpha=0.5. NU500AwC2.R implements our alignment optimization procedure for the correctly specified model under this condition. NU500AwCM2.R implements our alignment optimization procedure for the misspecified model under this condition. NU500BR.R implements our Bayesian adaptive Lasso procedure for the correctly specified model under this condition. NU500BRM.R implements our Bayesian adaptive Lasso procedure for the misspecified model under this condition.

The U500L1LM file folder includes two files for the condition of uniform DIF, I=500, Small DIF, 1/8 DIF and alpha=0.5. U500BR.R implements our Bayesian adaptive Lasso procedure for the correctly specified model under this condition. U500BRM.R implements our Bayesian adaptive Lasso procedure for the misspecified model under this condition.
