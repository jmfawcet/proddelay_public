
# Frequentist Figures -----------------------------------------------------

comb_anova_d$Descriptives %>%
  mutate(dv = "d'") %>%
  bind_rows(
    comb_anova_fam$Descriptives %>% mutate(dv = "Fam."),
    comb_anova_rho$Descriptives %>% mutate(dv = "Rec.")
  ) -> dat

plot_dat <- dat %>%
  mutate(
    exp = factor(exp,
                 levels = c("e1", "e2"),
                 labels = c("Experiment 1", "Experiment 2")),
    dv  = factor(dv,  levels = c("d'", "Fam.", "Rec.")),
    condition = factor(condition, 
                       levels = c("silent", "aloud"),
                       labels = c('Silent', 'Aloud')),
    soa = factor(soa,
                 levels = c("early", "late", "catch"),
                 labels = c("Early", "Late", "Catch")),
    err = 0.5 * FLSD
  )

p = ggplot(
  plot_dat,
  aes(x = soa, y = Mean, colour = condition, group = condition, shape = condition)
) +
  geom_errorbar(
    aes(ymin = Mean - err, ymax = Mean + err),
    position = position_dodge(width = 0.45),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.6
  ) +
  facet_grid(rows = vars(dv), cols = vars(exp), scales = "free_y") +
  labs(x = "Delay", y = "Estimate", colour = "", shape = "") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "black", colour = "black"),
    strip.text = element_text(colour = "white", face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  ) + 
  scale_color_manual(values = condition_colors) +   # Apply custom colors
  scale_shape_manual(values = condition_shapes) 

ggsave('figures/Figure S1.pdf', p, width=6, height=6, bg='transparent')
ggsave('figures/Figure S1.tiff', p, width=6, height=6, bg='transparent')


# Cumulative Models - Figure ----------------------------------------------

condition_labels = c('aloud.catch', 'aloud.early', 'aloud.late',
                     'silent.catch', 'silent.early', 'silent.late')

# Weights for combining the experiments
w1 = 32/(32+39)
w2 = 39/(32+39)

uv_sum <- fit_uv_comb$draws(variables = c('mu', 'mu_theta')) %>%
  as_draws_df() %>%
  select(.draw, `mu[1,1]`:`mu_theta[6,2]`) %>%
  mutate(across(`mu_theta[1,1]`:`mu_theta[6,2]`, exp)) %>%
  pivot_longer(
    cols = - .draw,
    names_to = c("param_base", "index", "exp"),
    names_pattern = "(mu(?:_theta)?)(?:\\[)(\\d+),(\\d+)\\]"
  ) %>%
  mutate(
    index = as.integer(index),
    exp = as.integer(exp),
    param_type = case_when(
      param_base == "mu" ~ "d'",
      param_base == "mu_theta" ~ "sigma"
    ),
    condition = condition_labels[index]
  ) %>%
  separate_wider_delim(condition, delim = ".", names = c("Condition", "SOA")) %>%
  mutate(
    Experiment = case_when(
      exp == 1 ~ "Experiment 1",
      exp == 2 ~ "Experiment 2"
    )
  ) %>%
  group_by(.draw, param_base, index, param_type, Condition, SOA) %>%
  mutate(
    value_combined = w1 * value[exp == 1] + w2 * value[exp == 2]
  ) %>%
  ungroup() %>%
  bind_rows(
    filter(., exp == 1) %>%
      mutate(Experiment = "Experiment 1"),
    
    filter(., exp == 2) %>%
      mutate(Experiment = "Experiment 2"),
    
    filter(., exp == 1) %>%
      mutate(
        value = value_combined,
        Experiment = "Combined"
      )
  ) %>%
  group_by(Condition, SOA, param_type, Experiment) %>%
  mean_qi(value, .width = 0.95) %>%
  select(-c(.width:.interval)) %>%
  mutate(
    model = "Unequal Variance",
    condition = Condition,
    estimate = value,
    lb = .lower,
    ub = .upper,
    exp = Experiment,
    type = param_type,
    soa = SOA
  ) %>%
  select(condition:soa)

uv_con_sum <- fit_uv_comb$draws(variables = c('mu', 'mu_theta')) %>%
  as_draws_df() %>%
  select(.draw, `mu[1,1]`:`mu_theta[6,2]`) %>%
  mutate(across(`mu_theta[1,1]`:`mu_theta[6,2]`, exp)) %>%
  mutate(
    mu_contrast_1_e1 = `mu[1,1]` - `mu[4,1]`,
    mu_contrast_2_e1 = `mu[2,1]` - `mu[5,1]`,
    mu_contrast_3_e1 = `mu[3,1]` - `mu[6,1]`,
    mu_theta_contrast_1_e1 = `mu_theta[1,1]` - `mu_theta[4,1]`,
    mu_theta_contrast_2_e1 = `mu_theta[2,1]` - `mu_theta[5,1]`,
    mu_theta_contrast_3_e1 = `mu_theta[3,1]` - `mu_theta[6,1]`,
    mu_contrast_1_e2 = `mu[1,2]` - `mu[4,2]`,
    mu_contrast_2_e2 = `mu[2,2]` - `mu[5,2]`,
    mu_contrast_3_e2 = `mu[3,2]` - `mu[6,2]`,
    mu_theta_contrast_1_e2 = `mu_theta[1,2]` - `mu_theta[4,2]`,
    mu_theta_contrast_2_e2 = `mu_theta[2,2]` - `mu_theta[5,2]`,
    mu_theta_contrast_3_e2 = `mu_theta[3,2]` - `mu_theta[6,2]`
  ) %>%
  mutate(
    mu_contrast_1_comb = w1 * mu_contrast_1_e1 + w2 * mu_contrast_1_e2,
    mu_contrast_2_comb = w1 * mu_contrast_2_e1 + w2 * mu_contrast_2_e2,
    mu_contrast_3_comb = w1 * mu_contrast_3_e1 + w2 * mu_contrast_3_e2,
    mu_theta_contrast_1_comb = w1 * mu_theta_contrast_1_e1 + w2 * mu_theta_contrast_1_e2,
    mu_theta_contrast_2_comb = w1 * mu_theta_contrast_2_e1 + w2 * mu_theta_contrast_2_e2,
    mu_theta_contrast_3_comb = w1 * mu_theta_contrast_3_e1 + w2 * mu_theta_contrast_3_e2
  ) %>%
  select(.draw, ends_with("_e1"), ends_with("_e2"), ends_with("_comb")) %>%
  pivot_longer(
    cols = - .draw,
    names_to = c("param_base", "index", "Experiment"),
    names_pattern = "(mu(?:_theta)?_contrast)_(\\d)_(.*)"
  ) %>%
  mutate(
    index = as.integer(index),
    type = case_when(
      param_base == "mu_contrast" ~ "d' Contrasts",
      param_base == "mu_theta_contrast" ~ "Sigma Contrasts"
    ),
    soa = case_when(
      index == 1 ~ "Catch",
      index == 2 ~ "Early",
      index == 3 ~ "Late"
    ),
    exp = case_when(
      Experiment == "e1" ~ "Experiment 1",
      Experiment == "e2" ~ "Experiment 2",
      Experiment == "comb" ~ "Combined"
    )
  ) %>%
  group_by(soa, type, exp) %>%
  mean_qi(value, .width = 0.95) %>%
  select(-c(.width:.interval)) %>%
  mutate(
    model = "Unequal Variance",
    condition = 'PE',
    estimate = value,
    lb = .lower,
    ub = .upper
  ) %>%
  select(soa:exp, condition:ub)

dp_sum <- fit_dp_comb$draws(variables = c('mu', 'mu_theta')) %>%
  as_draws_df() %>%
  select(.draw, `mu[1,1]`:`mu_theta[6,2]`) %>%  
  mutate(across(`mu_theta[1,1]`:`mu_theta[6,2]`, plogis)) %>%
  pivot_longer(
    cols = - .draw,
    names_to = c("param_base", "index", "exp"),
    names_pattern = "(mu(?:_theta)?)(?:\\[)(\\d+),(\\d+)\\]"
  ) %>%
  mutate(
    index = as.integer(index),
    exp = as.integer(exp),
    param_type = case_when(
      param_base == "mu" ~ "Fam.",
      param_base == "mu_theta" ~ "Rec."
    ),
    condition = condition_labels[index]
  ) %>%
  separate_wider_delim(condition, delim = ".", names = c("Condition", "SOA")) %>%
  mutate(
    Experiment = case_when(
      exp == 1 ~ "Experiment 1",
      exp == 2 ~ "Experiment 2"
    )
  ) %>%
  group_by(.draw, param_base, index, param_type, Condition, SOA) %>%
  mutate(
    value_combined = w1 * value[exp == 1] + w2 * value[exp == 2]
  ) %>%
  ungroup() %>%
  bind_rows(
    filter(., exp == 1) %>%
      mutate(Experiment = "Experiment 1"),
    
    filter(., exp == 2) %>%
      mutate(Experiment = "Experiment 2"),
    
    filter(., exp == 1) %>%
      mutate(
        value = value_combined,
        Experiment = "Combined"
      )
  ) %>%
  group_by(Condition, SOA, param_type, Experiment) %>%
  mean_qi(value, .width = 0.95) %>%
  select(-c(.width:.interval)) %>%
  mutate(
    model = "Dual Process",
    condition = Condition,
    estimate = value,
    lb = .lower,
    ub = .upper,
    exp = Experiment,
    type = param_type,
    soa = SOA
  ) %>%
  select(condition:soa)

dp_con_sum <- fit_dp_comb$draws(variables = c('mu', 'mu_theta')) %>%
  as_draws_df() %>%
  select(.draw, `mu[1,1]`:`mu_theta[6,2]`) %>%
  mutate(across(`mu_theta[1,1]`:`mu_theta[6,2]`, plogis)) %>%
  mutate(
    mu_contrast_1_e1 = `mu[1,1]` - `mu[4,1]`,
    mu_contrast_2_e1 = `mu[2,1]` - `mu[5,1]`,
    mu_contrast_3_e1 = `mu[3,1]` - `mu[6,1]`,
    mu_theta_contrast_1_e1 = `mu_theta[1,1]` - `mu_theta[4,1]`,
    mu_theta_contrast_2_e1 = `mu_theta[2,1]` - `mu_theta[5,1]`,
    mu_theta_contrast_3_e1 = `mu_theta[3,1]` - `mu_theta[6,1]`,
    mu_contrast_1_e2 = `mu[1,2]` - `mu[4,2]`,
    mu_contrast_2_e2 = `mu[2,2]` - `mu[5,2]`,
    mu_contrast_3_e2 = `mu[3,2]` - `mu[6,2]`,
    mu_theta_contrast_1_e2 = `mu_theta[1,2]` - `mu_theta[4,2]`,
    mu_theta_contrast_2_e2 = `mu_theta[2,2]` - `mu_theta[5,2]`,
    mu_theta_contrast_3_e2 = `mu_theta[3,2]` - `mu_theta[6,2]`
  ) %>%
  mutate(
    mu_contrast_1_comb = w1 * mu_contrast_1_e1 + w2 * mu_contrast_1_e2,
    mu_contrast_2_comb = w1 * mu_contrast_2_e1 + w2 * mu_contrast_2_e2,
    mu_contrast_3_comb = w1 * mu_contrast_3_e1 + w2 * mu_contrast_3_e2,
    mu_theta_contrast_1_comb = w1 * mu_theta_contrast_1_e1 + w2 * mu_theta_contrast_1_e2,
    mu_theta_contrast_2_comb = w1 * mu_theta_contrast_2_e1 + w2 * mu_theta_contrast_2_e2,
    mu_theta_contrast_3_comb = w1 * mu_theta_contrast_3_e1 + w2 * mu_theta_contrast_3_e2
  ) %>%
  select(.draw, ends_with("_e1"), ends_with("_e2"), ends_with("_comb")) %>%
  pivot_longer(
    cols = - .draw,
    names_to = c("param_base", "index", "Experiment"),
    names_pattern = "(mu(?:_theta)?_contrast)_(\\d)_(.*)"
  ) %>%
  mutate(
    index = as.integer(index),
    type = case_when(
      param_base == "mu_contrast" ~ "Fam. Contrasts",
      param_base == "mu_theta_contrast" ~ "Rec. Contrasts"
    ),
    soa = case_when(
      index == 1 ~ "Catch",
      index == 2 ~ "Early",
      index == 3 ~ "Late"
    ),
    exp = case_when(
      Experiment == "e1" ~ "Experiment 1",
      Experiment == "e2" ~ "Experiment 2",
      Experiment == "comb" ~ "Combined"
    )
  ) %>%
  group_by(soa, type, exp) %>%
  mean_qi(value, .width = 0.95) %>%
  select(-c(.width:.interval)) %>%
  mutate(
    model = "Dual Process",
    condition = 'PE',
    estimate = value,
    lb = .lower,
    ub = .upper
  ) %>%
  select(soa:exp, condition:ub)

combined_df <- bind_rows(uv_sum %>% filter(type!= 'sigma'), 
                         dp_sum, 
                         uv_con_sum %>% filter(type!= 'Sigma Contrasts'), 
                         dp_con_sum) %>%
  mutate(
    type = factor(type, levels = c("d'", "d' Contrasts", "Rec.", "Rec. Contrasts", "Fam.", "Fam. Contrasts")),
    condition = factor(condition, levels = c("Silent", "Aloud", "PE")),
    soa = factor(soa, levels=c('Early', 'Late', 'Catch')),
    experiment = factor(exp, levels = c('Experiment 1', 'Experiment 2', 'Combined'))
  )

p <- combined_df %>% ggplot(aes(x = soa, y = estimate, shape = condition, color = condition)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = 'grey', alpha = 0.2) +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1) +
  facet_grid(experiment~type, scales='free_x') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Estimate", x = "") +
  theme(
    strip.background = element_rect(fill = "black", color = NA),
    strip.text.x     = element_text(color = "white"),
    strip.text.y     = element_text(color = "white"),
    
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    
    legend.position = "top", 
    legend.title = element_blank(),
    
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.background = element_blank()
  ) + scale_color_manual(values = condition_colors) + 
  scale_shape_manual(values = condition_shapes) 

ggsave('figures/Figure 2 Behavioural Data.pdf', p, width=12, height=7, bg='transparent')
ggsave('figures/Figure 2 Behavioural Data.tiff', p, width=12, height=7, bg='transparent')
