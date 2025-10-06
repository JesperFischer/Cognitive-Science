

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] rating1;
  vector[N] rating2;
  vector[N] group;
  
}



transformed data{
  vector<lower=0, upper = 1>[N] rating11;
  vector<lower=0, upper = 1>[N] rating22;
  vector<lower=0, upper = 1>[N] groupp;
  
  rating11 = rating1/9;
  rating22 = rating2/9;
  groupp = group/9;
  
  
  
  
  
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0, upper = 1> bias;
  real<lower=0> kappa;
  real<lower=0> kappa2;

}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += beta_proportion_lpdf(bias |0.5,2);
  target += lognormal_lpdf(kappa |3,0.5);
  target += lognormal_lpdf(kappa2 |3,0.5);
  
  
  for(i in 1:N){
  target += beta_proportion_lpdf(rating11[i] | bias, kappa);
  target += beta_proportion_lpdf(rating22[i] | inv_logit(0.5*logit(rating11[i])+0.5*logit(groupp[i])), kappa2);
}
}


generated quantities{
  real prior_bias;
  real prior_kappa;
  real prior_kappa2;
  

  vector[N] sim_rating1;
  vector[N] sim_rating2;



  prior_bias = beta_proportion_rng(0.5,2);
  prior_kappa = lognormal_rng(3,0.5);
  prior_kappa2 = lognormal_rng(3,0.5);
    
  
  for(i in 1:N){
    sim_rating1[i] = beta_proportion_rng(bias, kappa);
    sim_rating2[i] = beta_proportion_rng(inv_logit(0.5*logit(rating11[i])+0.5*logit(groupp[i])), kappa2);
  }
}

