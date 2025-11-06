
#
# 
# Combine Data ------------------------------------------------------------
#
#

comb_test_dat = bind_rows(
    e1_test_dat, e2_test_dat
  ) %>%
  mutate(sid=paste(exp, sid))

e1_test_sum_total %>%
  mutate(exp='e1', sid=paste(exp, sid)) %>%
  bind_rows(e2_test_sum_total %>% mutate(exp='e2', sid=paste(exp, sid))) -> comb_test_sum_total

#
#
# Frequentist Models ------------------------------------------------------
#
#

# Accuracy Summary Statistics ------------------------------------------------------

# Calculate summaries
comb_test_dat %>%
  group_by(sid, condition, soa, exp) %>%
  summarize(yes=mean(said_yes), n=n()) -> comb_test_sum

# Create a wide data frame with d'
comb_test_sum_wide = comb_test_sum %>%
  pivot_wider(id_cols=sid, names_from=c(condition, soa), values_from=yes) %>%
  rename('new'='foil_foil') %>%
  mutate(across(
    matches("^(aloud|silent)"),  # Apply to all hit rate columns
    ~ qnorm(d_prime_correct(., 25)) - qnorm(d_prime_correct(new, 75)),  # Compute d'
    .names = "d_{.col}"  # Name new columns with "d_"
  ))

# Create a long d' data frame
comb_test_d_sum = comb_test_sum_wide %>%
  select(sid, starts_with('d_')) %>%
  pivot_longer(cols=c(d_silent_early:d_aloud_catch), names_to='condition', values_to='d') %>%
  mutate(condition = str_remove(condition, "^d_")) %>%
  separate(condition, into = c("condition", "soa"), sep = "_")

# ANOVAs --------------------------------------------------------

# Analysis of the "hits"
comb_anova_acc = ezANOVA(comb_test_sum_total %>% filter(condition!='foil', exp=='e2'), .(yes), .(sid), .(condition, soa))
comb_anova_acc$Descriptives = ezStats(comb_test_sum_total %>% filter(condition!='foil', exp=='e2'), .(yes), .(sid), .(condition, soa)) %>%
  mutate(SE = SD / sqrt(N))

# Analysis of the d'
comb_anova_d = ezANOVA(comb_test_sum_total %>% filter(condition != 'foil'), .(d), .(sid), .(condition, soa), between=exp)
comb_anova_d$Descriptives = ezStats(comb_test_sum_total %>% filter(condition != 'foil'), .(d), .(sid), .(condition, soa), between=exp)

# Analysis of the Familiarity
comb_anova_fam = ezANOVA(comb_test_sum_total %>% filter(condition != 'foil'), .(df), .(sid), .(condition, soa), between=exp)
comb_anova_fam$Descriptives = ezStats(comb_test_sum_total %>% filter(condition != 'foil'), .(df), .(sid), .(condition, soa), between=exp)

# Analysis of the Rho
comb_anova_rho = ezANOVA(comb_test_sum_total %>% filter(condition != 'foil'), .(rho), .(sid), .(condition, soa), between=exp)
comb_anova_rho$Descriptives = ezStats(comb_test_sum_total %>% filter(condition != 'foil'), .(rho), .(sid), .(condition, soa), between=exp)

#
#
# Bayesian Models ------------------------------------------------------
#
#

# NOTE: These models are included in case anyone wishes to compare the
# current results against Fawcett et al. (2025). However, we now use a
# more advanced ROC curve based statistical approach detailed in the 
# following script.

# Model d' by Experiment
comb_test_dat = comb_test_dat %>% mutate(comb_condition = interaction(condition, soa))
m2_comb_d = brm(said_yes ~ comb_condition:exp-1 + (comb_condition-1 | sid) + (comb_condition-1 | word),
                family=brms::bernoulli(link='probit'),
                backend = 'cmdstan',
                cores = ncores_brm,
                chains = nchains_brm,
                prior = c(e1_acc_priors[c(-1:-2),], prior(normal(.5,.5), class='b'), prior(normal(-.5,.5), class='b', coef='comb_conditionfoil.foil:expe1'), prior(normal(-.5,.5), class='b', coef='comb_conditionfoil.foil:expe2')),
                iter = niter_brm,
                sample_prior = 'yes',
                init = 0,
                threads = threading(4),
                #control = list(adapt_delta = .95),
                file = 'models/comb_exp_bayes_d',
                file_refit = 'on_change',
                data=comb_test_dat
)

# Estimates and Contrasts combining Experiments
emmeans(m2_comb_d, ~comb_condition)

# Analysis of Rho
temp = comb_test_sum_total %>% 
  filter(condition != 'foil') %>%
  mutate(comb_condition = interaction(condition, soa))

comb_bayes_rho1 = brm(rho ~ comb_condition:exp-1 + (1 | sid),
                    backend = 'cmdstan',
                    cores = ncores_brm,
                    chains = nchains_brm,
                    prior = e1_rho_priors,
                    iter = niter_brm,
                    sample_prior = 'yes',
                    threads = threading(4),
                    #control = list(adapt_delta = .999, max_treedepth=20),
                    file = 'models/comb_bayes_rho1',
                    file_refit = 'on_change',
                    data=temp)

# Estimates and Contrasts combining Experiments
emmeans(comb_bayes_rho1, ~comb_condition|exp)

# Analysis of the Familiarity
temp = comb_test_sum_total %>% 
  filter(condition != 'foil') %>%
  mutate(comb_condition = interaction(condition, soa))

comb_bayes_f1 = brm(df ~ comb_condition:exp-1 + (1 | sid),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e1_f_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  #control = list(adapt_delta = .999, max_treedepth=20),
                  file = 'models/comb_bayes_f1',
                  file_refit = 'on_change',
                  data=temp)

# Estimates and Contrasts combining Experiments
emmeans(comb_bayes_f1, ~comb_condition)
