
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
  real<lower=0> sigma;
  real<lower=0> lambda[year3-year1+N+1]; 
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
  
  vector[N] lambda_sum_2006;
  vector[N] lambda_sum_2011;
  vector[N] lambda_sum_2017;  
  vector[N1] lambda_sum_1_2006;
  vector[N1] lambda_sum_1_2011;
  vector[N1] lambda_sum_1_2017;

  // Age classes >1  
  for(i in 1:N){
    // Deriving expressions for 2017 - Correcting contribution of FOI in 2017
    ma_2017[i] = m0*exp(-omega*Decimal_age[i]);
    pr_neg_2017[i] = (1-m0)*exp(-sum(lambda[2:(i+1)])-0.5*lambda[1]); 
    // Terms j<i
    for(j in 0:(i-1)){
      lambda_sum_2017[j+1] = 0;
      if(j<(i-1)){
        for(k in (j+1):(i-1)){
          lambda_sum_2017[j+1] += lambda[i+1-k]; 
        }
      }
      lambda_sum_2017[j+1] += 0.5*lambda[i+1-i];
      pr_neg_2017[i] += m0*exp(-j*omega)*(omega*(exp(-omega)-exp(-lambda[i+1-j]))/(lambda[i+1-j]-omega))*exp(-lambda_sum_2017[j+1]);
    }
    // Term j=i
    pr_neg_2017[i] += m0*exp(-i*omega)*(omega*(exp(-0.5*omega)-exp(-0.5*lambda[i+1-i]))/(lambda[i+1-i]-omega));
    pr_pos_2017[i] = 1 - pr_neg_2017[i];

    // Deriving expressions for 2011
    ma_2011[i] = m0*exp(-omega*Decimal_age[i]);
    pr_neg_2011[i] = (1-m0)*exp(-sum(lambda[(year3-year2+1):(year3-year2+1+i)])-0.5*lambda[year3-year2+1]);
    // Terms j<i
    for(j in 0:(i-1)){
      lambda_sum_2011[j+1] = 0;
      if(j<(i-1)){
        for(k in (j+1):(i-1)){
          lambda_sum_2011[j+1] += lambda[year3-year2+i+1-k]; 
        }
      }
      lambda_sum_2011[j+1] += 0.5*lambda[year3-year2+i+1-i]; 
      pr_neg_2011[i] += m0*exp(-j*omega)*(omega*(exp(-omega)-exp(-lambda[year3-year2+i+1-j]))/(lambda[year3-year2+i+1-j]-omega))*exp(-lambda_sum_2011[j+1]);
    }
    // Term j=i
    pr_neg_2011[i] += m0*exp(-i*omega)*(omega*(exp(-0.5*omega)-exp(-0.5*lambda[year3-year2+i+1-i]))/(lambda[year3-year2+i+1-i]-omega));
    pr_pos_2011[i] = 1 - pr_neg_2011[i];
    
    // Deriving expressions for 2006
    ma_2006[i] = m0*exp(-omega*Decimal_age[i]);
    pr_neg_2006[i] = (1-m0)*exp(-sum(lambda[(year3-year1+1):(year3-year1+1+i)])-0.5*lambda[year3-year1+1]);
    // Terms j<i
    for(j in 0:(i-1)){
      lambda_sum_2006[j+1] = 0;
      if(j<(i-1)){
        for(k in (j+1):(i-1)){
          lambda_sum_2006[j+1] += lambda[year3-year1+i+1-k]; 
        }
      }
      lambda_sum_2006[j+1] += 0.5*lambda[year3-year1+i+1-i]; 
      pr_neg_2006[i] += m0*exp(-j*omega)*(omega*(exp(-omega)-exp(-lambda[year3-year1+i+1-j]))/(lambda[year3-year1+i+1-j]-omega))*exp(-lambda_sum_2006[j+1]);
    }
    // Term j=i
    pr_neg_2006[i] += m0*exp(-i*omega)*(omega*(exp(-0.5*omega)-exp(-0.5*lambda[year3-year1+i+1-i]))/(lambda[year3-year1+i+1-i]-omega));
    pr_pos_2006[i] = 1 - pr_neg_2006[i];
  }
  
  // Age classes <1
  for(i in 1:N1){
    // Deriving expressions for 2016 - I need to control when to use lambda[2016] or lambda[2015]
    ma_1_2017[i] = m0*exp(-omega*Decimal_age_under1[i]);
    pr_neg_1_2017[i] = (1-m0)*exp(-(i-0.5)/12.0*lambda[year3-year3_1+1]);
    // Terms j<i
    if(i>1){
      for(j in 1:(i-1)){
        lambda_sum_1_2017[j] = 0;
        if(j<(i-1)){
          for(k in (j+1):(i-1)){
            lambda_sum_1_2017[j] += 1.0/12*lambda[year3-year3_1+i+1-i];
          }
        }
        lambda_sum_1_2017[j] += 0.5/12*lambda[year3-year3_1+i+1-i];
        pr_neg_1_2017[i] += m0*exp(-(j-1)/12.0*omega)*(omega*(exp(-1.0/12*omega)-exp(-1.0/12*lambda[year3-year3_1+i+1-i]))/(lambda[year3-year3_1+i+1-i]-omega))*exp(-lambda_sum_1_2017[j]);
      }
    }
    // Term j=i
    pr_neg_1_2017[i] += m0*exp(-(i-1)/12.0*omega)*(omega*(exp(-0.5/12*omega)-exp(-0.5/12*lambda[year3-year3_1+i+1-i]))/(lambda[year3-year3_1+i+1-i]-omega));
    pr_pos_1_2017[i] = 1 - pr_neg_1_2017[i];

    // Deriving expressions for 2011 - I need to control when to use lambda[2011] or lambda[2010]
    ma_1_2011[i] = m0*exp(-omega*Decimal_age_under1[i]);
    pr_neg_1_2011[i] = (1-m0)*exp(-(i-0.5)/12.0*lambda[year3-year2+1]);
    // Terms j<i
    if(i>1){
      for(j in 1:(i-1)){
        lambda_sum_1_2011[j] = 0;
        if(j<(i-1)){
          for(k in (j+1):(i-1)){
            lambda_sum_1_2011[j] += 1.0/12*lambda[year3-year2+i+1-i];
          }
        }
        lambda_sum_1_2011[j] += 0.5/12*lambda[year3-year2+i+1-i];
        pr_neg_1_2011[i] += m0*exp(-(j-1)/12.0*omega)*(omega*(exp(-1.0/12*omega)-exp(-1.0/12*lambda[year3-year2+i+1-i]))/(lambda[year3-year2+i+1-i]-omega))*exp(-lambda_sum_1_2011[j]);
      }
    }
    // Term j=i
    pr_neg_1_2011[i] += m0*exp(-(i-1)/12.0*omega)*(omega*(exp(-0.5/12*omega)-exp(-0.5/12*lambda[year3-year2+i+1-i]))/(lambda[year3-year2+i+1-i]-omega));
    pr_pos_1_2011[i] = 1 - pr_neg_1_2011[i];
  
    // Deriving expressions for 2006 - I need to control when to use lambda[2006] or lambda[2005]
    ma_1_2006[i] = m0*exp(-omega*Decimal_age_under1[i]);
    pr_neg_1_2006[i] = (1-m0)*exp(-(i-0.5)/12.0*lambda[year3-year1+1]);
    // Terms j<i
    if(i>1){
      for(j in 1:(i-1)){
        lambda_sum_1_2006[j] = 0;
        if(j<(i-1)){
          for(k in (j+1):(i-1)){
            lambda_sum_1_2006[j] += 1.0/12*lambda[year3-year1+i+1-i];
          }
        }
        lambda_sum_1_2006[j] += 0.5/12*lambda[year3-year1+i+1-i];
        pr_neg_1_2006[i] += m0*exp(-(j-1)/12.0*omega)*(omega*(exp(-1.0/12*omega)-exp(-1.0/12*lambda[year3-year1+i+1-i]))/(lambda[year3-year1+i+1-i]-omega))*exp(-lambda_sum_1_2006[j]);
      }
    }
    // Term j=i
    pr_neg_1_2006[i] += m0*exp(-(i-1)/12.0*omega)*(omega*(exp(-0.5/12*omega)-exp(-0.5/12*lambda[year3-year1+i+1-i]))/(lambda[year3-year1+i+1-i]-omega));
    pr_pos_1_2006[i] = 1 - pr_neg_1_2006[i];
  }
  
}

// The model to be estimated. 
model {
  sigma ~ exponential(1);
  lambda[1] ~ normal(0, 0.5); 
  for(i in 2:(year3-year1+N+1))
    lambda[i] ~ normal(lambda[i-1], sigma);
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
