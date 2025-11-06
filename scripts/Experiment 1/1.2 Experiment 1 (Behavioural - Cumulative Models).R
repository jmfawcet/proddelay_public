
# Load Experiment 1 Data
source('scripts/newscripts/Experiment 1/1.0 Experiment 1 (Read in Data).R')

#
# Compile Stan Models
#

# Equal Variance
m.ev.1 <- cmdstan_model('scripts/newscripts/stan_models/evsdt-1-ll.stan')

# Unequal Variance
m.uv.1 <- cmdstan_model('scripts/newscripts/stan_models/uvsdt-1-ll.stan')

# Dual Process
m.dp.1 <- cmdstan_model('scripts/newscripts/stan_models/dpsdt-1-ll.stan')

# Mixture Signal Detection Model
m.mix.1 <- cmdstan_model('scripts/newscripts/stan_models/mixsdt-1-ll.stan')

# Variable Recollection Dual Process Model
m.fullmix.1 <- cmdstan_model('scripts/newscripts/stan_models/fullmixsdt-1-ll.stan')

# Prepare Data for Stan
e1_test_dat = e1_test_dat %>% 
  mutate(comb_condition = interaction(condition, soa)) %>%
  mutate(instruction = case_when(
    comb_condition == 'aloud.catch' ~ 1,
    comb_condition == 'aloud.early' ~ 2,
    comb_condition == 'aloud.late' ~ 3,
    comb_condition == 'silent.catch' ~ 4,
    comb_condition == 'silent.early' ~ 5,
    comb_condition == 'silent.late' ~ 6,
    comb_condition == 'foil.foil' ~ 7
  ),
  subject = as.numeric(factor(sid)),
  is_old = ifelse(condition == 'foil', 0, 1),
  response = key
  )

# Create Stan Data
sdata_e1 <- list(
  N = nrow(e1_test_dat),
  group = e1_test_dat$instruction,
  K = length(unique(e1_test_dat$response)),
  y = e1_test_dat$response,
  subject = e1_test_dat$subject,
  S = length(unique(e1_test_dat$subject)),
  ind = 7,
  N_g = length(unique(e1_test_dat$instruction)),
  is_old = e1_test_dat$is_old
)

# Fit Equal Variance Model
fit_ev_e1 <- m.ev.1$sample(
  data = sdata_e1,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_ev_e1$summary(variables = c('mu', 'pop_thresh',
                                'sigma_dprime', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# Fit Unequal Variance Model
fit_uv_e1 <- m.uv.1$sample(
  data = sdata_e1,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_uv_e1$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# Fit Dual Process SDT Model
fit_dp_e1 <- m.dp.1$sample(
  data = sdata_e1,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_dp_e1$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

print(fit_dp_e1, pars = c('mu', 'mu_theta', 'pop_thresh',
                          'sigma_dprime', 'sigma_theta', 'sigma_thresh'))

# Fit Guessing/Contamination SDT Mixture Model
fit_mix_e1 <- m.mix.1$sample(
  data = sdata_e1,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_mix_e1$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# Fit Full SDT Mixture Model/Continuous DP
fit_fullmix_e1 <- m.fullmix.1$sample(
  data = sdata_e1,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_fullmix_e1$summary(variables = c('mu', 'mu_2', 'mu_theta', 'pop_thresh',
                               'sigma_dprime', 'sigma_dprime_2', 'sigma_theta', 
                               'sigma_thresh')) %>% 
  print(n=100)
