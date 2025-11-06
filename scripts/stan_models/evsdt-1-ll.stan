data {
  int<lower=1> N;
  int<lower=1> N_g;
  array[N] int<lower=1, upper=N_g> group;
  int<lower=1> K;
  array[N] int<lower=1, upper=K> y;
  int<lower=1> S;
  array[N] int<lower=1, upper=S> subject;
  int<lower=1> ind;
  array[N] int<lower=0, upper=1> is_old;
}

parameters {
  vector[N_g-1] mu;
  real beta_thresh1;
  vector[K-2] beta_log_gaps;
  
  matrix[S, N_g-1] u_dprime_raw;
  matrix[S, K-1] u_thresh_raw;
  
  vector<lower=0>[N_g-1] sigma_dprime;
  vector<lower=0>[K-1] sigma_thresh;
  
  cholesky_factor_corr[N_g-1] L_corr_dprime;
  cholesky_factor_corr[K-1] L_corr_thresh;
}

transformed parameters {
  matrix[S, N_g-1] u_dprime = u_dprime_raw * diag_pre_multiply(sigma_dprime, L_corr_dprime)';
  matrix[S, K-1] u_thresh = u_thresh_raw * diag_pre_multiply(sigma_thresh, L_corr_thresh)';
  
  matrix[S, N_g-1] dprime;
  array[S] vector[K-1] thresh;
  
  for (s in 1:S) {
    dprime[s] = mu' + u_dprime[s];
    
    thresh[s][1] = beta_thresh1 + u_thresh[s,1];
    for (k in 2:(K-1)) {
      thresh[s][k] = thresh[s][k-1] + exp(beta_log_gaps[k-1] + u_thresh[s,k]);
    }
  }
}

model {
  mu ~ normal(1, 1);
  beta_thresh1 ~ normal(0, 1);
  beta_log_gaps ~ normal(-1, 0.7);
  sigma_dprime ~ normal(0, 0.5);
  sigma_thresh ~ normal(0, 0.5);
  
  to_vector(u_dprime_raw) ~ std_normal();
  to_vector(u_thresh_raw) ~ std_normal();
  
  L_corr_dprime ~ lkj_corr_cholesky(4);
  L_corr_thresh ~ lkj_corr_cholesky(4);
  
  for (n in 1:N) {
    int s = subject[n];
    
    if (is_old[n]) {
      y[n] ~ ordered_probit(dprime[s, group[n]], thresh[s]);
    } else {
      y[n] ~ ordered_probit(0, thresh[s]);
    }
  }
}

generated quantities {
  matrix[N_g-1, N_g-1] corr_dprime = multiply_lower_tri_self_transpose(L_corr_dprime);
  matrix[K-1, K-1] corr_thresh = multiply_lower_tri_self_transpose(L_corr_thresh);
  
  vector[K-1] pop_thresh;
  pop_thresh[1] = beta_thresh1;
  for (k in 2:(K-1)) {
    pop_thresh[k] = pop_thresh[k-1] + exp(beta_log_gaps[k-1]);
  }
  
  // array[N] real log_lik;
  // for (n in 1:N) {
  //   int s = subject[n];
  //   
  //   if (is_old[n]) {
  //     log_lik[n] = ordered_probit_lpmf(y[n] | dprime[s, group[n]], thresh[s]);
  //   } else {
  //     log_lik[n] = ordered_probit_lpmf(y[n] | 0, thresh[s]);
  //   }
  // }
}
