
data {
  int<lower=0> year1;
  int<lower=0> year2;
  int<lower=0> year3;
  int<lower=0> year3_1;
  int<lower=0> N;
  int<lower=0> N1;
  real<lower=0> m0;
  real Decimal_age[N];
  real Decimal_age_under1[N1];  
  int n_2006[N];
  int z_2006[N];
  int n_2011[N];
  int z_2011[N];
  int n_2017[N];
  int z_2017[N];
  int n1_2006[N1];
  int z1_2006[N1];
  int n1_2011[N1];
  int z1_2011[N1];
  int n1_2017[N1];
  int z1_2017[N1];
}

// The parameters accepted by the model. 
parameters {
  real<lower=0> lambda; 
  real<lower=0> omega;
}

transformed parameters{
  vector[N] pr_pos_2006; 
  vector[N] pr_pos_2011; 
  vector[N] pr_pos_2017; 
  vector[N1] pr_pos_1_2006; 
  vector[N1] pr_pos_1_2011; 
  vector[N1] pr_pos_1_2017;
  
  vector[N] pr_neg_2006; 
  vector[N] pr_neg_2011; 
  vector[N] pr_neg_2017;
  vector[N1] pr_neg_1_2006; 
  vector[N1] pr_neg_1_2011; 
  vector[N1] pr_neg_1_2017;
  
  vector[N] ma_2006;
  vector[N] ma_2011;
  vector[N] ma_2017;
  vector[N1] ma_1_2006;
  vector[N1] ma_1_2011;
  vector[N1] ma_1_2017;

  // Age classes >1  
  for(i in 1:N){
    ma_2006[i] = m0*exp(-omega*Decimal_age[i]);
    pr_neg_2006[i] = m0*(omega*(exp(-omega*Decimal_age[i])-exp(-lambda*Decimal_age[i])))/(lambda-omega) + exp(-Decimal_age[i]*lambda)*(1-m0);
    pr_pos_2006[i] = 1 - pr_neg_2006[i];
    
    ma_2011[i] = m0*exp(-omega*Decimal_age[i]);
    pr_neg_2011[i] = m0*(omega*(exp(-omega*Decimal_age[i])-exp(-lambda*Decimal_age[i])))/(lambda-omega) + exp(-Decimal_age[i]*lambda)*(1-m0);
    pr_pos_2011[i] = 1 - pr_neg_2011[i];
    
    ma_2017[i] = m0*exp(-omega*Decimal_age[i]);
    pr_neg_2017[i] = m0*(omega*(exp(-omega*Decimal_age[i])-exp(-lambda*Decimal_age[i])))/(lambda-omega) + exp(-Decimal_age[i]*lambda)*(1-m0);
    pr_pos_2017[i] = 1 - pr_neg_2017[i];
  }
  
  // Age classes <1
  for(i in 1:N1){
    ma_1_2006[i] = m0*exp(-omega*Decimal_age_under1[i]);
    pr_neg_1_2006[i] = m0*(omega*(exp(-omega*Decimal_age_under1[i])-exp(-lambda*Decimal_age_under1[i])))/(lambda-omega) + exp(-Decimal_age_under1[i]*lambda)*(1-m0);
    pr_pos_1_2006[i] = 1 - pr_neg_1_2006[i];
    
    ma_1_2011[i] = m0*exp(-omega*Decimal_age_under1[i]);
    pr_neg_1_2011[i] = m0*(omega*(exp(-omega*Decimal_age_under1[i])-exp(-lambda*Decimal_age_under1[i])))/(lambda-omega) + exp(-Decimal_age_under1[i]*lambda)*(1-m0);
    pr_pos_1_2011[i] = 1 - pr_neg_1_2011[i];
    
    ma_1_2017[i] = m0*exp(-omega*Decimal_age_under1[i]);
    pr_neg_1_2017[i] = m0*(omega*(exp(-omega*Decimal_age_under1[i])-exp(-lambda*Decimal_age_under1[i])))/(lambda-omega) + exp(-Decimal_age_under1[i]*lambda)*(1-m0);
    pr_pos_1_2017[i] = 1 - pr_neg_1_2017[i];
  }
  
}

// The model to be estimated. 
model {
  lambda ~ exponential(1);
  omega ~ exponential(1); 
  
  
  for (i in 1:N)
    z_2006[i] ~ binomial(n_2006[i], pr_pos_2006[i]);

  for (i in 1:N)
    z_2011[i] ~ binomial(n_2011[i], pr_pos_2011[i]);

  for (i in 1:N)
    z_2017[i] ~ binomial(n_2017[i], pr_pos_2017[i]);
  
  for (i in 1:N1)
    z1_2006[i] ~ binomial(n1_2006[i], pr_pos_1_2006[i]);

  for (i in 1:N1)
    z1_2011[i] ~ binomial(n1_2011[i], pr_pos_1_2011[i]);

  for (i in 1:N1)
    z1_2017[i] ~ binomial(n1_2017[i], pr_pos_1_2017[i]);
  
}

generated quantities{
  int z_sim_2006[N];
  int z_sim_2011[N];
  int z_sim_2017[N];
  int z1_sim_2006[N1];
  int z1_sim_2011[N1];
  int z1_sim_2017[N1];
  vector[3*(N+N1)] log_likelihood;
  
  // Generate posterior predictive checks
  for(i in 1:N)
      z_sim_2006[i] = binomial_rng(n_2006[i], pr_pos_2006[i]);
  
  for(i in 1:N)
      z_sim_2011[i] = binomial_rng(n_2011[i], pr_pos_2011[i]);
      
  for(i in 1:N)
      z_sim_2017[i] = binomial_rng(n_2017[i], pr_pos_2017[i]);
      
  for(i in 1:N1)
      z1_sim_2006[i] = binomial_rng(n1_2006[i], pr_pos_1_2006[i]);

  for(i in 1:N1)
      z1_sim_2011[i] = binomial_rng(n1_2011[i], pr_pos_1_2011[i]);

  for(i in 1:N1)
      z1_sim_2017[i] = binomial_rng(n1_2017[i], pr_pos_1_2017[i]);
      
  // Compute log-likelihood
  for(i in 1:N)
    log_likelihood[i] = binomial_lpmf(z_2006[i] | n_2006[i], pr_pos_2006[i]);

  for(i in 1:N)
    log_likelihood[N+i] = binomial_lpmf(z_2011[i] | n_2011[i], pr_pos_2011[i]);

  for(i in 1:N)
    log_likelihood[2*N+i] = binomial_lpmf(z_2017[i] | n_2017[i], pr_pos_2017[i]);
    
  for(i in 1:N1)
    log_likelihood[3*N+i] = binomial_lpmf(z1_2006[i] | n1_2006[i], pr_pos_1_2006[i]);
    
  for(i in 1:N1)
    log_likelihood[3*N+N1+i] = binomial_lpmf(z1_2011[i] | n1_2011[i], pr_pos_1_2011[i]);
        
  for(i in 1:N1)
    log_likelihood[3*N+2*N1+i] = binomial_lpmf(z1_2017[i] | n1_2017[i], pr_pos_1_2017[i]);

}
