

// The input data are vector 'CN' of length 'nrows', matrix 'Temp' of dimensions
// 'nrows' by 'nlag', and vector 'alpha' of length 'nlag'.
data {
  
  int<lower=0> nrows; // number of observations
  int<lower=0> nlag; // number of months into the past
  vector [nrows] CN; // observations - nutritional content
  matrix [nrows, nlag] Temp; // predictor - raw monthly seawater bottom temperature
  vector [nlag] alpha; // needed for dirichilet prior below - necessary here for syntax
  // based on the number of months alpha should be c(rep(1/months, months))
  
}

// This model estimates four parameters: 'b', 'b0', 'w', and 'sigma'.
parameters {
  
  real b0; // intercept
  real<lower=0> b; // combined weight of prior months
  simplex [nlag] wA; // separate weights of prior months
  real<lower=0> sigma; // observation error must be positive
  
}

// Antecedent summing occurs in the transformed parameters block.
// This model also estimates two transformed parameters: 'TempAnt' and 'TempTemp'.
transformed parameters{
  
  vector  [nrows] TempAnt; // estimate of monthly weighted sum of antcedent temperatures
  
for (i in 1:nrows){ // for every observation (i)

  matrix  [nrows,nlag] TempTemp; // estimate of transformed temperature (based on weight)

  for (t in 1:nlag) { // for every month of lag into the past (t),
   // estimate monthly weights and transformed temperatures
      TempTemp[i,t] = wA[t]*Temp[i,t];
    }
  
  TempAnt[i] = sum(TempTemp[i,]); // sum all antecedent values for observation i
  
  }
  
}

// Linear model goes in the model block
model {
  
  for (i in 1:nrows){
    // nutritional content of giant kelp as a function of antecedent seawater temperature
    CN[i] ~ normal(b0 + b*TempAnt[i], sigma); // likelihood
  }
  
  // model priors - keeping fairly uninformative for now
  
  // note in JAGS, normal dist. is (mu, tau) where mu is mean and tau is precision
  // tau is equal to 1/sigma^2 or 1/variance 
  // so take care to re-calculate priors when converting between JAGS and STAN
  
  b0 ~ normal(0, 10); // intercept - in JAGS, listed as (0, 1E-2)
  b ~ normal(0, 10); // covariate parameter prior - in JAGS, listed as (0, 1E-2)
  wA ~ dirichlet(alpha); // give vector of weights (antecedent climate) dirichlet priors
  // see data - will be fed in with data list
  // sigma ~ gamma(0.001, 0.001); // error - keep same scale/value as JAGS script (alpha, beta)
  
}

