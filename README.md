# LIRT-MI
Adjusting for longitudinal measurement noninvariance in longitudinal item response theory models: Comparing Bayesian regularization and alignment optimization

By accounting for longitudinal measurement noninvariance in longitudinal item response theory models, we proposed a longitudinal IRT model with time-specific item parameters and additionally incorporated
the dependence among the repeated item measurements. Based on the principle of approximate invariance, Bayesian regularization and alignment optimization methods are developed for model estimation.

Uniform-1000SP20.R implements our Bayesian adaptive Lasso procedure for the N=1000, 20% DIF and small DIF condition in simulation study 1.

Nonuniform-1000SP20.R implements our Bayesian adaptive Lasso procedure for the N=1000, 20% DIF and small DIF condition in simulation study 2.

Uniform-example-CFA+DIF.R uses a confirmatory simple-structure MIRT model for DIF detection in the heuristic example.

Uniform-example-EML1+DIF.R first identifies the item-trait structure by the EML1 method and then uses the structure as confirmatory for DIF detection in the heuristic example.

Uniform-example-our.R uses our Bayesian adaptive Lasso procedure for simultaneously detecting item-trait relationship and DIF in the heuristic example.

Uniform-realdata.R implements real data analysis using our Bayesian adaptive Lasso procedure for the uniform DIF condition.

Parameter-MSE.pdf presents MSEs for each parameter estimate in study 2.
