# Getting Started

1. Modify line 10 in bias_cov_test.m to point to your local irt path.
2. Run bias_cov_test.m. On a 3.5Ghz desktop with 32GB RAM, this script took about 10 minutes to run. 
3. To compute PERK bias and covariance for the holdout-optimized PERK hyperparameters, comment out lines 128, 130, 132 and uncomment lines 129, 131, 133. Note that increasing the number of training samples (`rff.K`) will increase memory demands and run times. 
