# BayesHiddenPottsMixture
This repository is used for access the performance of the Bayesian hidden Potts mixture model, which was proposed in the manuscript titled "A Bayesian hidden Potts mixture model for analyzing lung cancer pathology images." Before running the code, please install R package Rcpp.

First, please run "data_loader.R" to load simulated or real data. This is a required step.

Then, run "model_fitting.R" to fit the proposed model and obtain the preliminary results, such as runtime, estimated model parameters, and some plots. 

*Note 1, we also provided 120 simulated datasets and their true model parameters on figshare. The download link is https://figshare.com/projects/Bayesian_hidden_Potts_mixture_models/29659

*Note 2, the notations in the code are followed the notations in the manuscript.

*Note 3, for data with large number of points and for model with large lattice size, it takes longer time. Please first try to run the MCMC algorithm with small lattice size (e.g. W = L = 20) or small number of iterations (e.g. iter = 1000).