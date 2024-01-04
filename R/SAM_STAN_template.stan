

// The input data are vectors 'CN' and 'Temp' of length 'N'.
data {
  
  int<lower=0> nrows; // number of months of data
  int<lower=0> n.lag; // number of months into the past
  vector[N] CN; // nutritional content
  vector[N] Temp; // raw monthly seawater bottom temperature
  vector [n.lag] alpha;
  
}

// The parameters accepted by the model. Our model
// accepts four parameters 'b', 'b0', 'w', and 'sigma'.
parameters {
  
  real<lower=0> b; // combined weight of prior months
  real b0; // intercept
  simplex [n.lag] wA; // separate weights of prior months
  real<lower=0> tau; // observation error
  
}

// the antecedent summing goes in the transformed parameters block

transformed parameters{
  
  for (t in 1:n.lag) { // lag identity
      TempTemp[i,t] = wA[t]*Temp[i,t]; // Temp is input data as a matrix of measurement by month
    }
  
  TempAnt[i] = sum(TempTemp[i,]); // sum all the antecedent monthly values for measurement i
  
  }
  
}

// the linear model goes in the model block

model {
  
  for (i in 13:nrows){
    CN[i] ~ normal(b0 + b*TempAnt[i], tau); // likelihood
  }
  
  // model priors - keeping fairly uninformative for now
  b0 ~ normal(0, 1E-2); // intercept
  b ~ normal(0, 1E-2); // covariate parameter prior
  wA ~ dirichlet(alpha); // give vector of weights (antecedent climate) dirichlet priors
  tau ~ gamma(0.001, 0.001); // error
  
}

