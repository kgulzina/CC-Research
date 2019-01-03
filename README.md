# CC-Research

My research goal is to build an emulator (i.e., surrogate) for the Water Erosion Predictor Project (WEPP) model using Gaussian processes and to increase runtime rate required to get a spatio-temporal distribution of soil loss across one specific field. We build a surrogate on a hierarchical Gaussian processes. We chose an exponential covariance function with the highly correlated weight parameters to assess the correlation within the response given the correlation between the functional inputs. We assume that weight parameters come from autoregressive model. Additionally, the model employs an empirical Bayes approach for parameters estimation.




# Backlog
- Derive the real gradient function for maximization of loglkl_with_penalty()
- Split weight_predictor() and loglkl_with_penalty() into separate files
- Find new a distribution which takes into account facts about weights: 1. High correlation 2.Concentration around zero 3.Positiveness 
- Write a function to assess the accuracy of loglkl_with_penalty() with different pars
- Open a separate new "log" document to keep track of changes and experiments with simulated data and preliminary model assumptions



# In-progress
- Following R-style guide
- Removing constants from functions
- Writing an abstract
- Separating simulation from functions and other constants



# Done
- Increased the penalty: rho = 0.99
- Set white noise variance: sigmasq = 1
- Changed w(t) from deterministic to sampled
- Simulated w(t) from truncated distribution
