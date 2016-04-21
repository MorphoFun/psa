 data {
   // Define variables in data
   // Number of observations (an integer)
   int<lower=0> nObs;
   // Number of parameters
   int<lower=0> p;
   // Variables
   real w[nObs];
   real TotalL[nObs];
   real AlarExtent[nObs];
   real Weight[nObs];
   real SkullL[nObs];
   real HumerusL[nObs];
   real FemurL[nObs];
   real TibioTarL[nObs];
   real SkullW[nObs];
   real SternumL[nObs];
 }
 
 parameters {
   // Define parameters to estimate
   real beta[p];
   
   // standard deviation (a positive real number)
   real<lower=0> sigma;
 }
 
 transformed parameters  {
   // Mean
   real mu[nObs];
   for (i in 1:nObs) {
     mu[i] <- beta[1] + beta[2]*TotalL[i] + beta[3]*AlarExtent[i] + beta[4]*Weight[i] + beta[5]*SkullL[i] + beta[6]*HumerusL[i] + beta[7]*FemurL[i] + beta[8]*TibioTarL[i] + beta[9]*SkullW[i] + beta[10]*SternumL[i]; 
   }

 }
 
 model {
   // Prior part of Bayesian inference (flat if unspecified)
   // Phenotypic data should be scaled to a mean of zero and unit variance
   // Taking lognormal of phenotypic data to linearize the data that should be normally distributed
     sigma ~ cauchy(0,1);
     beta ~ normal(0,1);
 
   // Likelihood part of Bayesian inference
     w ~ normal(mu, sigma);  

 }
