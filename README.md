# CC-Research

My research goal is to build an emulator (i.e., surrogate) for the Water Erosion Predictor Project (WEPP) model using Gaussian processes and to increase runtime rate required to get a spatio-temporal distribution of soil loss across one specific field. We build a surrogate on a hierarchical Gaussian processes. We chose an exponential covariance function with the highly correlated weight parameters to assess the correlation within the response given the correlation between the functional inputs. We assume that weight parameters come from autoregressive model. Additionally, the model employs an empirical Bayes approach for parameters estimation.




# Backlog
- Derive the real gradient function for maximization of loglkl_with_penalty()
- Split weight_predictor() and loglkl_with_penalty() into separate files
- Find new a distribution which takes into account facts about weights: 1. High correlation 2.Concentration around zero 3.Positiveness 


# In-progress
- Following R-style guide
- Removing constants from functions



# Done
- Changed w(t) from deterministic to sampled
- Simulated w(t) from truncated distribution.
