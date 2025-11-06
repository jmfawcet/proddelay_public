
#
# Windowed analysis ----------------------------------------------------------------
# 

# Preprocessing
# Early = 500 - 3500, Late = 3500 - 7000 Windows; this roughly corresponds to
# word/instruction onset and response onset
e2_study_pd %>%
  filter(bin >= 500 & bin <= 7000) %>%
  mutate(window = ifelse(bin <= 3500, 'early', 'late')) %>%
  group_by(sid, condition, soa, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e2_study_pd %>%
      filter(bin >= 500 & bin <= 7000) %>%
      group_by(sid, condition, soa) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% arrange(sid, condition, window)-> e2_study_pd_win

e2_study_pd %>%
  filter(bin >= 500 & bin <= 7000) %>%
  mutate(window = ifelse(bin <= 3500, 'early', 'late')) %>%
  group_by(sid, word, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e2_study_pd %>%
      filter(bin >= 500 & bin <= 7000) %>%
      group_by(sid, word) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% pivot_wider(names_from='window', values_from='pupil_z') -> e2_study_pd_per_trial

# Analysis of the windowed pupil
e2_study_pd_win_both = e2_study_pd_win %>% filter(window != 'complete')
e2_anova_win_both = ezANOVA(e2_study_pd_win_both , .(pupil_z), .(sid), .(condition, soa, window))
e2_anova_win_both$Descriptives = ezStats(e2_study_pd_win_both, .(pupil_z), .(sid), .(condition, soa, window))

e2_study_pd_win_comp = e2_study_pd_win %>% filter(window == 'complete')
e2_anova_win = ezANOVA(e2_study_pd_win_comp , .(pupil_z), .(sid), .(condition, soa))
e2_anova_win$Descriptives = ezStats(e2_study_pd_win_comp, .(pupil_z), .(sid), .(condition, soa))

#
#
# Functional Data Analysis ------------------------------------------------
#
#

# Fit the mass univariate analysis to each delay condition
e2_a_m_s_e = fit_massu_t(e2_study_pd_sum %>% filter(soa=='early'), 'Early', c('aloud', 'silent'))
e2_a_m_s_l = fit_massu_t(e2_study_pd_sum %>% filter(soa=='late'), 'Late', c('aloud', 'silent'))
e2_a_m_s_c = fit_massu_t(e2_study_pd_sum %>% filter(soa=='catch'), 'Catch', c('aloud', 'silent'))

# For each delay, calculate cluster level p-values; this is done by generating
# time series with similar autocorrelation reflecting the t-statistics
# assuming no real effects and then calculating the summed t values across 
# clusters that emerge. This simulated distributation is used to calculate p.
e2_massu_a_m_s_e = ar1_cluster_correction_t(e2_a_m_s_e,
                                        alpha          = 0.05,
                                        n_sim          = 6000,
                                        two_sided      = TRUE,
                                        df = length(unique(e2_study_pd_sum$sid))-1
)
e2_massu_a_m_s_e$observed_clusters = e2_massu_a_m_s_e$observed_clusters %>% filter(significant)

e2_a_m_s_e = e2_a_m_s_e %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e2_massu_a_m_s_e$observed_clusters$start_tp <= .x & .x <= e2_massu_a_m_s_e$observed_clusters$end_tp)))


e2_massu_a_m_s_l = ar1_cluster_correction_t(e2_a_m_s_l,
                                        alpha          = 0.05,
                                        n_sim          = 6000,     
                                        two_sided      = TRUE,
                                        df = length(unique(e2_study_pd_sum$sid))-1
)
e2_massu_a_m_s_l$observed_clusters = e2_massu_a_m_s_l$observed_clusters %>% filter(significant)

e2_a_m_s_l = e2_a_m_s_l %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e2_massu_a_m_s_l$observed_clusters$start_tp <= .x & .x <= e2_massu_a_m_s_l$observed_clusters$end_tp)))


e2_massu_a_m_s_c = ar1_cluster_correction_t(e2_a_m_s_c,
                                        alpha          = 0.05,
                                        n_sim          = 6000,    
                                        two_sided      = TRUE,
                                        df = length(unique(e2_study_pd_sum$sid))-1
)
e2_massu_a_m_s_c$observed_clusters = e2_massu_a_m_s_c$observed_clusters %>% filter(significant)

e2_a_m_s_c = e2_a_m_s_c %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e2_massu_a_m_s_c$observed_clusters$start_tp <= .x & .x <= e2_massu_a_m_s_c$observed_clusters$end_tp)))

# Combine into a single data frame
e2_massu = e2_a_m_s_e %>%
  bind_rows(e2_a_m_s_l) %>%
  bind_rows(e2_a_m_s_c)

#
#
# GAMM --------------------------------------------------------------------
#
#

# Analysis ----------------------------------------------------------------

# Fit the GAMM model
e2_study_pd = e2_study_pd %>% mutate(comb_condition=interaction(condition, soa))
e2_study_pgam_ww = fit_gam_study_delayed(e2_study_pd, 'models/e2_study_pgam.rds')

#
#
# Predicting Study Phase Accuracy -----------------------------------------
#
#

# Basic Accuracy Model
temp = e2_test_dat %>%
  filter(condition != 'foil') %>%
  left_join(e2_study_pd_per_trial) %>%
  filter(!is.na(early) & !is.na(late)) %>%
  mutate(early = scale(early)[,1], late = scale(late)[,1])

e2_bayes_m1_pupilpred = brm(said_yes ~ condition:soa + condition:soa:early + condition:soa:late-1 + (condition:soa-1 | sid) + (condition:soa-1 | word),
                            family=brms::bernoulli(link='probit'),
                            backend = 'cmdstan',
                            cores = ncores_brm,
                            chains = nchains_brm,
                            prior = c(e1_acc_priors[c(-1:-2),], prior(normal(.5,.5), class='b')),
                            iter = niter_brm,
                            sample_prior = 'yes',
                            init=0,
                            threads = threading(4),
                            #control = list(adapt_delta = .95),
                            file = 'models/e2_bayes_m1_pupilpred',
                            file_refit = 'on_change',
                            data=temp)

# Basic Accuracy Model using overall pupil size, not windowed
temp = e2_test_dat %>%
  filter(condition != 'foil') %>%
  left_join(e2_study_pd_per_trial) %>%
  filter(!is.na(early) & !is.na(late)) %>%
  mutate(early = scale(early)[,1], late = scale(late)[,1], complete = scale(complete)[,1])

e2_bayes_m1_pupilpred_comp = brm(said_yes ~ condition:soa + condition:soa:complete-1 + (condition:soa-1 | sid) + (condition:soa-1 | word),
                                 family=brms::bernoulli(link='probit'),
                                 backend = 'cmdstan',
                                 cores = ncores_brm,
                                 chains = nchains_brm,
                                 prior = c(e1_acc_priors[c(-1:-2),], prior(normal(.5,.5), class='b')),
                                 iter = niter_brm,
                                 init=0,
                                 sample_prior = 'yes',
                                 threads = threading(4),
                                 control = list(adapt_delta = .95),
                                 file = 'models/e2_bayes_m1_pupilpred_comp',
                                 file_refit = 'on_change',
                                 data=temp)

#
#
# Rolling Window ----------------------------------------------------------
#
#

# Fit a rolling window model where 100 ms windows are averaged and used to predict
# later memory performance via MLM.
e2_rolling_window = fit_model_across_timepoints_cluster_delayed(e2_study_pd, seq(50, 6950,10), 'models/e2_rolling.rds')

# Apply a similar clustering approach to that used above, only assuming z rather than t
e2_rolling_window_cluster_cse <- ar1_cluster_correction(
  results_df     = e2_rolling_window,
  effect_z_col   = "silent_early_z",
  effect_p_col   = "silent_early_p",
  alpha          = 0.05,
  n_sim          = 6000,
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster_csl <- ar1_cluster_correction(
  results_df     = e2_rolling_window,
  effect_z_col   = "silent_late_z",
  effect_p_col   = "silent_late_p",
  alpha          = 0.05,
  n_sim          = 6000,
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster_csc <- ar1_cluster_correction(
  results_df     = e2_rolling_window,
  effect_z_col   = "silent_catch_z",
  effect_p_col   = "silent_catch_p",
  alpha          = 0.05,
  n_sim          = 6000,
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster_cae <- ar1_cluster_correction(
  results_df     = e2_rolling_window,
  effect_z_col   = "aloud_early_z",
  effect_p_col   = "aloud_early_p",
  alpha          = 0.05,
  n_sim          = 6000,
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster_cal <- ar1_cluster_correction(
  results_df     = e2_rolling_window,
  effect_z_col   = "aloud_late_z",
  effect_p_col   = "aloud_late_p",
  alpha          = 0.05,
  n_sim          = 6000,
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster_cac <- ar1_cluster_correction(
  results_df     = e2_rolling_window,
  effect_z_col   = "aloud_catch_z",
  effect_p_col   = "aloud_catch_p",
  alpha          = 0.05,
  n_sim          = 6000,
  two_sided      = TRUE,
  z_thresh       = 1.96
)

#
#
# Correlating Behavioural and Pupillometric PE ----------------------------
#
#

# Calculate the production effect for early, late and catch
e2_study_pd_win %>%
  pivot_wider(id_cols=sid, names_from=c(condition, soa, window), values_from=pupil_z) %>%
  mutate(ppe_e = aloud_early_late-silent_early_late, ppe_l = aloud_late_late-silent_late_late, ppe_c = aloud_catch_late-silent_catch_late) %>%
  select(sid, ppe_e, ppe_l, ppe_c) -> e2_pupillary_pe_sum

# Create a data frame containing behavioural and pupillary production effects
e2_test_sum_total %>%
  filter(condition != 'foil') %>%
  select(-n, -yes) %>%
  pivot_wider(id_cols=sid, names_from=c(condition, soa), values_from=c(rho, df, d)) %>%
  mutate(dpe_e = d_aloud_early-d_silent_early, dpe_l = d_aloud_late-d_silent_late, dpe_c = d_aloud_catch-d_silent_catch,
         rpe_e = rho_aloud_early-rho_silent_early, rpe_l = rho_aloud_late-rho_silent_late, rpe_c = rho_aloud_catch-rho_silent_catch,
         fpe_e = df_aloud_early - df_silent_early, fpe_l = df_aloud_late - df_silent_late, fpe_c = df_aloud_catch - df_silent_catch) %>%
  select_at(vars(contains('pe'))) %>%
  right_join(e2_pupillary_pe_sum) -> e2_behav_pe_sum

# Create a data frame suitable for correlations
temp = e2_behav_pe_sum %>%
  ungroup() %>% 
  select(-sid)

# Comparisons we care about
comparisons <- list(c("dpe_e", "ppe_e"), c("dpe_l", "ppe_l"), c("dpe_c", "ppe_c"))

# For each correlation model, check for outliers using both approaches and
# fit correlations using brms to remaining data; due to our focus
# on the response period, correlations are fit using the PPE in the 
# late window
cor_results1 <- compute_correlations_with_outlier_removal(temp, comparisons)
cor_results2 <- compute_correlations_with_outlier_removal_w_pastpriors(temp, comparisons)
