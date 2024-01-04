model{
  
# Likelihood --------------------------------------------------------------

  #our data are normal (or lognormal) C:N ratios
  #at a monthly scale
  for(i in 1:nrows){ # for the number of rows in the dataset of C:n ratios

    #C:N is normally distributed with mean mu[i], and variance tau
    y[i] ~ dnorm(mu[i], tau)
    
    #regression: what might impact mu
    mu[i] <- b0 + #slope
      b[1]*TempAnt[i] #antecedent temperature term
      #could consider adding to this regression:
    #random effects of site
    #other covariates (e.g., NPP)
# Antecedent summing ------------------------------------------------------

    #sum all the antecedent monthly values for measurement i
    TempAnt[i] <- sum(TempTemp[i,])
    
    #getting each month's weight to sum above
    for(t in 1:n.lag){ #number of months into past
      #Temp is input data as a matrix of measurement by month
      TempTemp[i,t] <- Temp[i,t]*wA[t]
    }
    


# Missing temp data imputing (for now) ------------------------------------

    for(t in 1:n.lag){
      Temp[i,t] ~ dnorm(mu.tmp[t], tau.tmp[t])
    }
    
# Goodness of fit parameters ----------------------------------------------

    #replicated data
    yrep[i] ~ dnorm(mu[i], tau)
    
    #residuals
    resid[i] <- y[i] - mu[i]
    
  }

# Priors ------------------------------------------------------------------

  #All priors are relatively uninformative, vague priors
  #intercept and overall variance
  b0 ~ dnorm(0, 1E-2)
  tau ~ dgamma(0.001, 0.001)
  
  #priors for covariates
  for(i in 1:n.covs){
    b[i] ~ dnorm(0, 1E-2)
  }
  
  #Antecedent climate priors
  #sum of weights for climate lag to divide so that they equal 1
  sumA <- sum(deltaA[])
  
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.lag){
    wA[t] <- deltaA[t]/sumA
    deltaA[t] ~ dgamma(1,1)
  }
  
  #MIssing data priors
  
  for(t in 1:n.lag){
    mu.tmp[t] ~ dnorm(0, 1E-2)
    sig.tmp[t] ~ dunif(0,500)
    tau.tmp[t] <- pow(sig.tmp[t],-2)
  }
  

  
}