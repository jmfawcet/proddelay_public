
# Load Experiment 2 Data
source('scripts/newscripts/Experiment 2/2.0 Experiment 2 (Read in Data).R')

#
# Compile Stan Models (same as Experiment 1)
#

m.ev.1 <- cmdstan_model('scripts/newscripts/stan_models/evsdt-1-ll.stan')

m.uv.1 <- cmdstan_model('scripts/newscripts/stan_models/uvsdt-1-ll.stan')

m.dp.1 <- cmdstan_model('scripts/newscripts/stan_models/dpsdt-1-ll.stan')

m.mix.1 <- cmdstan_model('scripts/newscripts/stan_models/mixsdt-1-ll.stan')

m.fullmix.1 <- cmdstan_model('scripts/newscripts/stan_models/fullmixsdt-1-ll.stan')

# set up data
e2_test_dat = e2_test_dat %>% 
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

# create stan data for e2
sdata_e2 <- list(
  N = nrow(e2_test_dat),
  group = e2_test_dat$instruction,
  K = length(unique(e2_test_dat$response)),
  y = e2_test_dat$response,
  subject = e2_test_dat$subject,
  S = length(unique(e2_test_dat$subject)),
  ind = 7,
  N_g = length(unique(e2_test_dat$instruction)),
  is_old = e2_test_dat$is_old
)

# fit equal variance SDT model
fit_ev_e2 <- m.ev.1$sample(
  data = sdata_e2,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_ev_e2$summary(variables = c('mu', 'pop_thresh',
                                'sigma_dprime', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit unequal variance SDT model
fit_uv_e2 <- m.uv.1$sample(
  data = sdata_e2,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_uv_e2$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit dual process SDT mixture model
fit_dp_e2 <- m.dp.1$sample(
  data = sdata_e2,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_dp_e2$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit guessing/contamination SDT mixture model
fit_mix_e2 <- m.mix.1$sample(
  data = sdata_e2,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_mix_e2$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                 'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                 'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit full SDT mixture model/continuous DP
fit_fullmix_e2 <- m.fullmix.1$sample(
  data = sdata_e2,
  chains = 14,
  parallel_chains = 14,
  iter_warmup = 8000,
  iter_sampling = 10000,
  init = 0,
  refresh = 100
)

fit_fullmix_e2$summary(variables = c('mu', 'mu_2', 'mu_theta', 'pop_thresh',
                                     'sigma_dprime', 'sigma_dprime_2', 'sigma_theta', 
                                     'sigma_thresh')) %>% 
  print(n=100)
