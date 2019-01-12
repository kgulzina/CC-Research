# CC-Research

My research goal is to build an emulator (i.e., surrogate) for the Water Erosion Predictor Project (WEPP) model using Gaussian processes and to increase runtime rate required to get a spatio-temporal distribution of soil loss across one specific field. We build a surrogate on a hierarchical Gaussian processes. We chose an exponential covariance function with the highly correlated weight parameters to assess the correlation within the response given the correlation between the functional inputs. We assume that weight parameters come from autoregressive model. Additionally, the model employs an empirical Bayes approach for parameters estimation.




## Backlog
- Derive the real gradient function for maximization of loglkl_with_penalty()
- Matrix differentiation techniques: Book(!)/article(!)/Wolfram Alpha(?) - probably download
- Write a function to assess the runtime of optim()




## In-progress
- Following R-style guide
- Writing an abstract
- Writing log process of optimization

- Find appropriate proposal distribution for Metropolis-Hastings algorithm/Block sampling method
- Transform w into log(w), sample from induced distribution using Stan
- Found a new distribution which takes into account facts about weights: 1. High correlation 2.Concentration around zero 3.Positiveness (Beta autoregressive process(yes), or Gamma(?)) >>>> But how to integrate it ??? <<<<




## Done
- Increased the penalty: rho = 0.99
- Set white noise variance: sigmasq = 1 
- Changed w(t) from deterministic to sampled
- Simulated w(t) from truncated distribution
- Splitted weight_predictor() and loglkl_with_penalty() into separate files
- Separated simulation from functions and other constants
- A separate new "log" document was created to keep track of changes and experiments with simulated data and preliminary model assumptions
- Removed constants from functions
- Plotted 3D of loglkl_with_penalty holding all but two constant (too lengthy)
- Wrote a function to assess the accuracy of loglkl_with_penalty() with different pars, gradients
- Changed w(t) into log scale: log(w) ~ MVN centered at -1 or -2
