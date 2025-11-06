
# Pupil Plots -------------------------------------------------------------


# From the GAMM model, derive smoothed predictions
study_pgam_plot <- plot_smooth(
  e2_study_pgam_ww, 
  view      = "bin",
  main      = "Estimated Pupil Response Over Time",
  rug       = TRUE,
  n.grid    = 1000,
  plot_all  = "comb_condition"
)$fv %>%
  dplyr::select(comb_condition, bin, fit, CI, ll, ul) %>%
  dplyr::mutate(exp = 'E2', model = 'GAMM')

# Summarize the raw pupil data
#    (Mean, standard error, confidence interval, etc.)
massu_summary <- e2_study_pd %>%
  dplyr::group_by(bin, comb_condition) %>%
  dplyr::summarize(
    m  = mean(pupil_z),
    se = sd(pupil_z) / sqrt(dplyr::n()),
    ci = se * 1.96,
    .groups = 'drop'
  ) %>%
  dplyr::mutate(model = 'Mass Univariate')

# Combine the raw data summary with the GAMM predictions
study_pd_sum_figure <- massu_summary %>%
  dplyr::bind_rows(
    study_pgam_plot %>%
      dplyr::mutate(se = CI / 1.96) %>%
      dplyr::select(
        bin,
        comb_condition,
        m   = fit,
        se,
        ci  = CI,
        model
      )
  ) %>%
  dplyr::mutate(
    model = factor(model, levels = c('Mass Univariate', 'GAMM')),
    exp   = 'E2'
  )

# Get significant bins from GAMM comparisons
gamm_plot_sig <- plot_diff(
  e2_study_pgam_ww,
  view  = 'bin',
  comp  = list(comb_condition = c('aloud.early', 'silent.early')),
  n.grid = 1000,
  plot   = FALSE
) %>%
  dplyr::select(est, bin, CI) %>%
  dplyr::mutate(
    comparison = 'Early',
    sig        = est + CI < 0 | est - CI > 0
  ) %>%
  dplyr::bind_rows(
    plot_diff(
      e2_study_pgam_ww,
      view  = 'bin',
      comp  = list(comb_condition = c('aloud.late', 'silent.late')),
      n.grid = 1000,
      plot   = FALSE
    ) %>%
      dplyr::select(est, bin, CI) %>%
      dplyr::mutate(
        comparison = 'Late',
        sig        = est + CI < 0 | est - CI > 0
      )
  ) %>%
  dplyr::bind_rows(
    plot_diff(
      e2_study_pgam_ww,
      view  = 'bin',
      comp  = list(comb_condition = c('aloud.catch', 'silent.catch')),
      n.grid = 1000,
      plot   = FALSE
    ) %>%
      dplyr::select(est, bin, CI) %>%
      dplyr::mutate(
        comparison = 'Catch',
        sig        = est + CI < 0 | est - CI > 0
      )
  ) %>%
  dplyr::filter(sig) %>%
  dplyr::mutate(
    y_offset = dplyr::case_when(
      comparison == "Early"   ~ -0.6,
      comparison == "Late"  ~ -0.6,
      comparison == "Catch" ~ -0.6
    ),
    exp   = 'E2',
    model = 'GAMM'
  ) %>%
  dplyr::select(bin, comparison, y_offset, exp, model)

# Filter the massu data for significant bins and tag them as 'Mass Univariate'
massu_sig <- e2_massu %>%
  dplyr::filter(sig_clust) %>%
  dplyr::mutate(
    y_offset = dplyr::case_when(
      comparison == "Early"   ~ -0.6,
      comparison == "Late"  ~ -0.6,
      comparison == "Catch" ~ -0.6
    ),
    exp   = 'E2',
    model = 'Mass Univariate'
  ) %>%
  tibble::rownames_to_column("original_names") %>%
  dplyr::mutate(original_names = dplyr::row_number()) %>%
  tibble::column_to_rownames("original_names") %>%
  dplyr::select(bin, comparison, y_offset, exp, model)

# Combine GAMM significance with Mass Univariate significance
plot_sig <- massu_sig %>%
  dplyr::bind_rows(gamm_plot_sig) %>%
  dplyr::mutate(
    model = factor(model, levels = c('Mass Univariate', 'GAMM'))
  )

ecomb_study_pd_sum_figure = study_pd_sum_figure %>% 
  separate(comb_condition, into=c('condition', 'soa')) %>%
  mutate(soa = factor(tolower(soa), levels=c('early', 'late', 'catch')))

ecomb_plot_sig = plot_sig %>%
  mutate(y_offset = case_when(
    y_offset == -.6 ~ -.63,
    y_offset == -.65 ~ -.7,
    y_offset == -.7 ~ -.77
  )) %>%
  rename(soa=comparison) %>%
  mutate(soa = factor(tolower(soa), levels=c('early', 'late', 'catch')))


arrow_data <- data.frame(
  condition = 'aloud',
  soa      = c("early", "late", 'catch', 'early', 'late'),  # matches facet levels
  x        = c(0, 0, 0, 3500, 4500),        # where arrow starts horizontally
  xend     = c(0, 0, 0, 3500, 4500),        # the same x, so it's vertical
  y        = c(0.7, 0.7, 0.7, 0.7, 0.7),        # arrow's start
  yend     = c(0.5, 0.5, 0.5, 0.5, 0.5),         # arrow's end (below start => downward arrow)
  label    = c('Word\nInstruction', 'Word\nInstruction', "Word\nInstruction", "Response", "Response")
) %>%
  mutate(vjust = -.35 + .125*str_count(label, '\n'), soa = factor(tolower(soa), levels=c('early', 'late', 'catch')))



# Define labels and y-offsets
label_positions <- data.frame(
  comparison = c("Comparison"),
  y_label_offset = c(-0.64)  # Slightly above segments for clarity
)

g1 = ecomb_study_pd_sum_figure %>%
  mutate() %>% 
  ggplot(aes(x = bin, y = m, fill=condition, color = condition, linetype = condition, ymin=m-ci, ymax=m+ci)) +
  geom_line() +
  geom_ribbon(alpha=.5) +
  geom_segment(data = ecomb_plot_sig,
               aes(x = bin, xend = bin, y = y_offset, yend = y_offset - 0.02),
               inherit.aes = FALSE, size = 1) +
  # geom_text(data = label_positions,
  #           aes(x = min(massu_sig$bin) + 600,  # Align all labels at the same x position
  #               y = y_label_offset, 
  #               label = comparison),
  #           inherit.aes = FALSE, hjust = 1, size = 2.5) +
  labs(x = "Time (ms)", y = "z(Pupil)") +
  theme_classic() +
  geom_hline(
    yintercept = 0,
    colour='grey',
    alpha=.2
  ) + 
  facet_grid(soa~model) +
  theme(
    # Make facet strips black with white text
    strip.background = element_rect(fill = "black", color = NA),
    strip.text.x     = element_text(color = "white"),
    strip.text.y     = element_text(color = "white"),
    
    
    # Transparent panel and plot background
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    
    # Make axis & legend text visible on transparent background
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    #legend.text      = element_text(color = "white"),
    #legend.title     = element_text(color = "white")
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  geom_segment(
    data = arrow_data,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    # Add arrow at the "last" end of the segment (or "first", "both")
    arrow = arrow(
      length = unit(0.2, "cm"),   # arrow size
      ends = "last",             # which end to place arrow on
      type = "closed"            # closed arrow head
    ),
    color = "black",
    size  = 1
  ) +
  geom_text(
    data = arrow_data,
    aes(x = x, y = y, label = label, vjust=vjust),
    inherit.aes = FALSE,
    #vjust = -.2,       # Nudges text above the arrow's start
    hjust = 0.5,      # Center horizontally
    size  = 2.5,        # Font size
    color = "black",
    lineheight = 0.9
  ) +
  coord_cartesian(clip = "off", ylim=c(-.8,1.6)) +
  scale_color_manual(values = c("aloud" = "#F8766D", "control" = "#7CAE00", "silent" = "#00BFC4"))  # Customize colors

ggsave('figures/Figure 3 Pupil Data.pdf', g1, width=12, height=14, bg='transparent')
ggsave('figures/Figure 3 Pupil Data.tiff', g1, width=12, height=14, bg='transparent')

