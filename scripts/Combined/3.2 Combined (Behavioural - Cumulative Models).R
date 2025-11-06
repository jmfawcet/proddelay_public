
# Load Experiment 1 and 2 Data
source('scripts/Experiment 1/1.0 Experiment 1 (Read in Data).R')
source('scripts/Experiment 2/2.0 Experiment 2 (Read in Data).R')

#
# Compile Stan Models (modified for multiple experiments)
#

m.ev.2 <- cmdstan_model('scripts/stan_models/evsdt-2-ll.stan')

m.uv.2 <- cmdstan_model('scripts/stan_models/uvsdt-2-ll.stan')

m.dp.2 <- cmdstan_model('scripts/stan_models/dpsdt-2-ll.stan')

m.mix.2 <- cmdstan_model('scripts/stan_models/mixsdt-2-ll.stan')

m.fullmix.2 <- cmdstan_model('scripts/stan_models/fullmixsdt-2-ll.stan')

# set up data

combdat <- bind_rows(
  e1_test_dat, e2_test_dat) %>%
  mutate(sid=paste(exp, sid)) %>%
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

# create stan data for combined data

subj_exp <- combdat %>% 
  group_by(subject) %>% 
  summarize(exp = first(exp)) %>% 
  pull(exp)

sdata_comb <- list(
  N = nrow(combdat),
  group = combdat$instruction,
  exp_condition = as.numeric(factor(combdat$exp)),
  N_exp = length(unique(combdat$exp)),
  K = length(unique(combdat$response)),
  y = combdat$response,
  subject = combdat$subject,
  S = length(unique(combdat$subject)),
  ind = 7,
  N_g = length(unique(combdat$instruction)),
  is_old = combdat$is_old,
  subj_exp = as.numeric(factor(subj_exp))
)

# fit equal variance SDT model

fit_ev_comb <- m.ev.2$sample(
  data = sdata_comb,
  chains =      5,
  parallel_chains =      5,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_ev_comb$summary(variables = c('mu', 'pop_thresh',
                                'sigma_dprime', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit unequal variance SDT model

fit_uv_comb <- m.uv.2$sample(
  data = sdata_comb,
  chains =      5,
  parallel_chains =      5,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_uv_comb$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# Example of how to extract data
mu1_df <- fit_uv_comb$draws("mu", format = "df") %>%
  select(!starts_with(".")) %>%
  set_names(outer(c('aloud.catch', 'aloud.early', 'aloud.late',
                    'silent.catch', 'silent.early', 'silent.late'), c('e1', 'e2'), paste, sep = "_") %>% as.vector())
summary(mu1_df)


# fit dual process SDT mixture model

fit_dp_comb <- m.dp.2$sample(
  data = sdata_comb,
  chains =      5,
  parallel_chains =      5,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_dp_comb$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit guessing/contamination SDT mixture model

fit_mix_comb <- m.mix.2$sample(
  data = sdata_comb,
  chains =      5,
  parallel_chains =      5,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_mix_comb$summary(variables = c('mu', 'mu_theta', 'pop_thresh',
                                 'sigma_dprime', 'sigma_theta', 'sigma_thresh',
                                 'beta_thresh1', 'beta_log_gaps')) %>%
  print(n=100)

# fit full SDT mixture model/continuous DP

fit_fullmix_comb <- m.fullmix.2$sample(
  data = sdata_comb,
  chains =      5,
  parallel_chains =      5,
  iter_warmup = 4000,
  iter_sampling = 5000,
  init = 0,
  refresh = 100
)

fit_fullmix_comb$summary(variables = c('mu', 'mu_2', 'mu_theta', 'pop_thresh',
                                     'sigma_dprime', 'sigma_dprime_2', 'sigma_theta', 
                                     'sigma_thresh')) %>% 
  print(n=100)

