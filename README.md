# CC-Research

My research goal is to build an emulator (i.e., surrogate) for the Water Erosion Predictor Project (WEPP) model using Gaussian processes and to increase runtime rate required to get a spatio-temporal distribution of soil loss across one specific field. We build a surrogate on a hierarchical Gaussian processes. We chose an exponential covariance function with the highly correlated weight parameters to assess the correlation within the response given the correlation between the functional inputs. We assume that weight parameters come from autoregressive model. Additionally, the model employs an empirical Bayes approach for parameters estimation.




## Backlog
- Derive the real gradient function for maximization of loglkl_with_penalty()
- Split weight_predictor() and loglkl_with_penalty() into separate files
- Find a new distribution which takes into account facts about weights: 1. High correlation 2.Concentration around zero 3.Positiveness (Beta autoregressive process(?), or Gamma(?))
- Write a function to assess the accuracy of loglkl_with_penalty() with different pars
- Transform w(t) into log scale: log(w) ~ MVN centered at -1 or -2
- Matrix differentiation techniques: Book(!)/article(!)/Wolfram Alpha(?) - probably download
- Plot 3D of loglkl_with_penalty holding all but two constant





## In-progress
- Following R-style guide
- Removing constants from functions
- Writing an abstract
- Separating simulation from functions and other constants
- A separate new "log" document was created to keep track of changes and experiments with simulated data and preliminary model assumptions



## Done
- Increased the penalty: rho = 0.99
- Set white noise variance: sigmasq = 1 
- Changed w(t) from deterministic to sampled
- Simulated w(t) from truncated distribution
