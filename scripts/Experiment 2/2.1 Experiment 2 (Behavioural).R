
#
#
# Frequentist Models ------------------------------------------------------
#
#

# Accuracy Summary Statistics ------------------------------------------------------

# Create a 'flattened' IV for use in the Bayesian d' model
e2_test_dat = e2_test_dat %>% mutate(comb_condition = interaction(condition, soa))

# Calculate summaries
e2_test_dat %>%
  group_by(sid, condition, soa) %>%
  summarize(yes=mean(said_yes), n=n()) -> e2_test_sum

# Create a wide data frame with d'
e2_test_sum_wide = e2_test_sum %>%
  pivot_wider(id_cols=sid, names_from=c(condition, soa), values_from=yes) %>%
  rename('new'='foil_foil') %>%
  mutate(across(
    matches("^(aloud|silent)"),  # Apply to all hit rate columns
    ~ qnorm(d_prime_correct(., 25)) - qnorm(d_prime_correct(new, 75)),  # Compute d'
    .names = "d_{.col}"  # Name new columns with "d_"
  ))

# Create a long d' data frame
e2_test_d_sum = e2_test_sum_wide %>%
  select(sid, starts_with('d_')) %>%
  pivot_longer(cols=c(d_silent_early:d_aloud_catch), names_to='condition', values_to='d') %>%
  mutate(condition = str_remove(condition, "^d_")) %>%
  separate(condition, into = c("condition", "soa"), sep = "_")

# Preprocess data for DPSDT
e2_dpsdt_dat = preprocess_dpsdt_delayed(e2_test_dat) %>%
  separate(comb_condition, into = c("condition", "soa"))

# Combining into a single test phrase
e2_test_sum %>%
  left_join(e2_test_d_sum) %>%
  left_join(e2_dpsdt_dat) %>%
  mutate(soa=factor(soa, levels=c('foil', 'early', 'late', 'catch')), condition=factor(condition, levels=c('foil', 'silent', 'aloud'))) -> e2_test_sum_total

# ANOVAs --------------------------------------------------------

# Analysis of the "hits"
e2_anova_acc = ezANOVA(e2_test_sum_total %>% filter(condition!='foil'), .(yes), .(sid), .(condition, soa))
e2_anova_acc$Descriptives = ezStats(e2_test_sum_total %>% filter(condition!='foil'), .(yes), .(sid), .(condition, soa)) %>%
  mutate(SE = SD / sqrt(N))

# Analysis of the d'
e2_anova_d = ezANOVA(e2_test_sum_total %>% filter(condition != 'foil'), .(d), .(sid), .(condition, soa))
e2_anova_d$Descriptives = ezStats(e2_test_sum_total %>% filter(condition != 'foil'), .(d), .(sid), .(condition, soa))

# Analysis of the Familiarity
e2_anova_fam = ezANOVA(e2_test_sum_total %>% filter(condition != 'foil'), .(df), .(sid), .(condition, soa))
e2_anova_fam$Descriptives = ezStats(e2_test_sum_total %>% filter(condition != 'foil'), .(df), .(sid), .(condition, soa))

# Analysis of the Rho
e2_anova_rho = ezANOVA(e2_test_sum_total %>% filter(condition != 'foil'), .(rho), .(sid), .(condition, soa))
e2_anova_rho$Descriptives = ezStats(e2_test_sum_total %>% filter(condition != 'foil'), .(rho), .(sid), .(condition, soa))

#
#
# Bayesian Models ------------------------------------------------------
#
#

# NOTE: These models are included in case anyone wishes to compare the
# current results against Fawcett et al. (2025). However, we now use a
# more advanced ROC curve based statistical approach detailed in the 
# following script.

# Model d'
e2_test_dat = e2_test_dat %>% mutate(comb_condition = interaction(condition, soa))
m1_e2_d = brm(said_yes ~ comb_condition-1 + (comb_condition-1 | sid) + (comb_condition-1 | word),
                   family=brms::bernoulli(link='probit'),
                   backend = 'cmdstan',
                   cores = ncores_brm,
                   chains = nchains_brm,
                   prior = c(e1_acc_priors[c(-1:-2),], prior(normal(.5,.5), class='b'), prior(normal(-.5,.5), class='b', coef='comb_conditionfoil.foil')),
                   iter = niter_brm,
                   sample_prior = 'yes',
                   threads = threading(4),
                   control = list(adapt_delta = .95),
                   file = 'models/e2_bayes_d',
                   file_refit = 'on_change',
                   data=e2_test_dat
)

# d' is the difference between foil.foil and each other condition; comparisons
# between conditions are differences in d'
pairs(emmeans(m1_e2_d, ~comb_condition))

# Analysis of the Rho
temp = e2_test_sum_total %>% filter(condition != 'foil')
e2_bayes_rho1 = brm(rho ~ condition:soa-1 + (1 | sid),
                    backend = 'cmdstan',
                    cores = ncores_brm,
                    chains = nchains_brm,
                    prior = e1_rho_priors,
                    iter = niter_brm,
                    sample_prior = 'yes',
                    threads = threading(4),
                    control = list(adapt_delta = .999, max_treedepth=20),
                    file = 'models/e2_bayes_rho1',
                    file_refit = 'on_change',
                    data=temp)

e2_bayes_rho1_means = emmeans(e2_bayes_rho1, ~condition|soa)
e2_bayes_rho1_comps = pairs(e2_bayes_rho1_means)

# Analysis of the Familiarity
temp = e2_test_sum_total %>% filter(condition != 'foil')
e2_bayes_f1 = brm(df ~ condition:soa-1 + (1 | sid),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e1_f_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .999, max_treedepth=20),
                  file = 'models/e2_bayes_f1',
                  file_refit = 'on_change',
                  data=temp)

e2_bayes_f1_means = emmeans(e2_bayes_f1, ~condition|soa)
e2_bayes_f1_comps = pairs(e2_bayes_f1_means)
