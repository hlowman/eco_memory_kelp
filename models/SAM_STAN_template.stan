

// The input data are vectors 'CN' and 'Temp' of length 'N'.
data {
  
  int<lower=0> nrows; // number of months of data
  int<lower=0> nlag; // number of months into the past (nweight)
  vector[nrows] CN; // observations - nutritional content
  vector[nrows] Temp; // predictor - raw monthly seawater bottom temperature
  vector [nlag] alpha; // needed for dirichilet prior below
  
}

// The parameters accepted by the model. Our model
// accepts four parameters 'b', 'b0', 'w', and 'sigma'.
parameters {
  
  real b0; // intercept
  real<lower=0> b; // combined weight of prior months
  simplex [nlag] wA; // separate weights of prior months
  real<lower=0> tau; // observation error must be positive
  
}

// Antecedent summing goes in the transformed parameters block
transformed parameters{
  
  vector  [nrows] TempAnt;
  
for (i in 1:nrows){ // need to skip over missing data at some point...

  matrix  [nrows,nlag] TempTemp;

  for (t in 1:nlag) { // for every month of lag into the past, 1 being present month
  // calculate the weighted Temperature (TempTemp)
      TempTemp[i,t] = wA[t]*Temp[i,t]; // Temp is input data as a matrix of measurement by month
    }
  
  TempAnt[i] = sum(TempTemp[i,]); // sum all the antecedent monthly values for observation i
  
  }
  
}

// Linear model goes in the model block
model {
  
  for (i in 1:nrows){
    CN[i] ~ normal(b0 + b*TempAnt[i], tau); // likelihood
  }
  
  // model priors - keeping fairly uninformative for now
  b0 ~ normal(0, 1E-2); // intercept
  b ~ normal(0, 1E-2); // covariate parameter prior
  wA ~ dirichlet(alpha); // give vector of weights (antecedent climate) dirichlet priors
  tau ~ gamma(0.001, 0.001); // error
  
}

