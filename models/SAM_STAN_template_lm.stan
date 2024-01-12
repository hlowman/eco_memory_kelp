

// The input data are vectors 'CN' and 'Temp' of length 'nrows'.
data {
  
  int<lower=0> nrows; // number of months of data
  vector[nrows] CN; // observations - nutritional content
  vector[nrows] Temp; // predictor - raw monthly seawater bottom temperature
  
}

// The parameters accepted by the model. Our model
// accepts three parameters 'b', 'b0', and 'tau'.
parameters {
  
  real b0; // intercept
  real b; // slope
  real<lower=0> tau; // observation error must be positive
  
}

// nothing in the transformed parameters block

transformed parameters{
  
}

// the linear model goes in the model block

model {
  
  // Likelihood
  CN ~ normal(b0 + b*Temp, tau);
  
  // model priors - keeping fairly uninformative for now
  b0 ~ normal(0, 1E-2); // intercept
  b ~ normal(0, 1E-2); // covariate parameter prior
  tau ~ gamma(0.001, 0.001); // error
  // remember, script MUST end in a blank line
  
}

