data {
  int<lower=1> N;              
  int<lower=1> N_g;            
  int<lower=1> N_exp;           
  int<lower=1> K;             
  int<lower=1> S;              
  
  array[N] int<lower=1, upper=K> y;         
  array[N] int<lower=1, upper=S> subject;   
  array[N] int<lower=1, upper=N_g> group;    
  array[N] int<lower=0, upper=1> is_old;   
  array[S] int<lower=1, upper=N_exp> subj_exp; 
}

parameters {
  array[N_exp, N_g - 1] ordered[2] mu_ordered;
  matrix[N_g-1, N_exp] mu_theta;             
  
  vector[N_exp] beta_thresh1;               
  matrix[K-2, N_exp] beta_log_gaps;        
  
  matrix[S, N_g-1] u_dprime_raw;        
  matrix[S, N_g-1] u_dprime_2_raw;         
  matrix[S, N_g-1] u_theta_raw;      
  matrix[S, K-1] u_thresh_raw;     
  
  vector<lower=0>[N_g-1] sigma_dprime;
  vector<lower=0>[N_g-1] sigma_dprime_2;
  vector<lower=0>[N_g-1] sigma_theta;
  vector<lower=0>[K-1] sigma_thresh;
  
  cholesky_factor_corr[N_g-1] L_corr_dprime;
  cholesky_factor_corr[N_g-1] L_corr_dprime_2;
  cholesky_factor_corr[N_g-1] L_corr_theta;
  cholesky_factor_corr[K-1] L_corr_thresh;
}

transformed parameters {
  matrix[N_g-1, N_exp] mu;
  matrix[N_g-1, N_exp] mu_2;
  for (exp_idx in 1:N_exp) {
    for (g in 1:(N_g-1)) {
      mu[g, exp_idx]   = mu_ordered[exp_idx, g][2];
      mu_2[g, exp_idx] = mu_ordered[exp_idx, g][1];
    }
  }
  
  matrix[S, N_g-1] u_dprime = u_dprime_raw * diag_pre_multiply(sigma_dprime, L_corr_dprime)';
  matrix[S, N_g-1] u_dprime_2 = u_dprime_2_raw * diag_pre_multiply(sigma_dprime_2, L_corr_dprime_2)';
  matrix[S, N_g-1] u_theta = u_theta_raw * diag_pre_multiply(sigma_theta, L_corr_theta)';
  matrix[S, K-1] u_thresh = u_thresh_raw * diag_pre_multiply(sigma_thresh, L_corr_thresh)';
  
  matrix[S, N_g-1] dprime;
  matrix[S, N_g-1] dprime_2;
  matrix[S, N_g-1] theta; 
  array[S] vector[K - 1] thresh;
  
  for (s in 1:S) {
    int exp_idx = subj_exp[s];            
    
    dprime[s] = mu[,exp_idx]' + u_dprime[s];
    dprime_2[s] = mu_2[,exp_idx]' + u_dprime_2[s];
    theta[s] = inv_logit(mu_theta[,exp_idx]' + u_theta[s]);
    
    thresh[s][1] = beta_thresh1[exp_idx] + u_thresh[s,1];
    for (k in 2:(K-1)) {
      thresh[s][k] = thresh[s][k-1] + exp(beta_log_gaps[k-1, exp_idx] + u_thresh[s,k]);
    }
  }
}

model {
  for (exp_idx in 1:N_exp) {
    for (g in 1:(N_g-1)) {
      mu_ordered[exp_idx, g][2] ~ normal(2, 1);
      mu_ordered[exp_idx, g][1] ~ normal(1, 1);
    }
  }
  
  to_vector(mu_theta) ~ normal(0, 0.5);
  beta_thresh1 ~ normal(0, 1);
  to_vector(beta_log_gaps) ~ normal(-1, 0.7);
  
  sigma_dprime ~ normal(0, 0.5);
  sigma_dprime_2 ~ normal(0, 0.5);
  sigma_theta ~ normal(0, 0.5);
  sigma_thresh ~ normal(0, 0.5);
  
  to_vector(u_dprime_raw) ~ std_normal();
  to_vector(u_dprime_2_raw) ~ std_normal();
  to_vector(u_theta_raw) ~ std_normal();
  to_vector(u_thresh_raw) ~ std_normal();
  
  L_corr_dprime ~ lkj_corr_cholesky(4);
  L_corr_dprime_2 ~ lkj_corr_cholesky(4);
  L_corr_theta ~ lkj_corr_cholesky(4);
  L_corr_thresh ~ lkj_corr_cholesky(4);
  
  for (n in 1:N) {
    int s = subject[n];
    
    if (is_old[n] && group[n] < N_g) {
      target += log_mix(
        theta[s, group[n]], 
        ordered_probit_lpmf(y[n] | dprime[s, group[n]], thresh[s]),
        ordered_probit_lpmf(y[n] | dprime_2[s, group[n]], thresh[s])
      );
    } else {
      target += ordered_probit_lpmf(y[n] | 0, thresh[s]);
    }
  }
}

generated quantities {
  matrix[N_g-1, N_g-1] corr_dprime = multiply_lower_tri_self_transpose(L_corr_dprime);
  matrix[N_g-1, N_g-1] corr_dprime_2 = multiply_lower_tri_self_transpose(L_corr_dprime_2);
  matrix[N_g-1, N_g-1] corr_theta = multiply_lower_tri_self_transpose(L_corr_theta);
  matrix[K-1, K-1] corr_thresh = multiply_lower_tri_self_transpose(L_corr_thresh);
  
  array[N_exp] vector[K-1] pop_thresh;
  for (exp_idx in 1:N_exp) {
    pop_thresh[exp_idx][1] = beta_thresh1[exp_idx];
    for (k in 2:(K-1)) {
      pop_thresh[exp_idx][k] = pop_thresh[exp_idx][k-1] + exp(beta_log_gaps[k-1, exp_idx]);
    }
  }
  
  array[N] real log_lik;
  for (n in 1:N) {
    int s = subject[n];

    if (is_old[n] && group[n] < N_g) {
      log_lik[n] = log_mix(
        theta[s, group[n]],
        ordered_probit_lpmf(y[n] | dprime[s, group[n]], thresh[s]),
        ordered_probit_lpmf(y[n] | dprime_2[s, group[n]], thresh[s])
      );
    } else {
      log_lik[n] = ordered_probit_lpmf(y[n] | 0, thresh[s]);
    }
  }
}
