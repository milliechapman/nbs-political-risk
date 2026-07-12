knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

# Install any missing packages before running:
# install.packages(c("tidyverse", "wbstats", "ggcorrplot", "ggrepel",
#                    "countrycode", "patchwork", "broom", "scales", "here"))

library(tidyverse)
library(wbstats)        # World Bank API
library(ggcorrplot)     # Correlation heatmaps
library(ggrepel)        # Non-overlapping labels on scatter plots
library(countrycode)    # Country code harmonisation
library(patchwork)      # Combine ggplots
library(broom)          # Tidy model output
library(scales)
library(here)

nbs_raw <- read_csv(here("data/all_pathways/all_pathways_tco2e_adm0.csv"))

# ── Pathway display labels ─────────────────────────────────────────────────────
pathway_labels <- c(
  # Individual pathways
  crp_cba_tco2e  = "Cropland Agroforestry",
  crp_rec_tco2e  = "Reduced Emissions (Crops)",
  crp_scc_tco2e  = "Soil Carbon (Crops)",
  for_afc_tco2e  = "Avoided Forest Conversion",
  for_csf_tco2e  = "Climate-Smart Forestry",
  for_ref_tco2e  = "Forest Reforestation",
  for_rwf_tco2e  = "Reduced Woodland Fire (Forest)",
  grs_agc_tco2e  = "Agroforestry (Grassland)",
  grs_asc_tco2e  = "Avoided Savanna Conversion",
  grs_grr_tco2e  = "Grassland Restoration",
  grs_reg_tco2e  = "Reforestation (Grassland)",
  grs_rwf_tco2e  = "Reduced Grassland Fire",
  grs_sba_tco2e  = "Savanna Burning Management",
  grs_scg_tco2e  = "Soil Carbon (Grassland)",
  grs_sfm_tco2e  = "Sustainable Forest Mgmt",
  wet_apc_tco2e  = "Avoided Peatland Conversion",
  wet_awc_tco2e  = "Avoided Wetland Conversion",
  wet_cwr_tco2e  = "Coastal Wetland Restoration",
  wet_ipm_tco2e  = "Improved Peatland Mgmt",
  wet_per_tco2e  = "Peatland Restoration",
  # Biome aggregates
  total_crop_tco2e    = "Total: Cropland",
  total_forest_tco2e  = "Total: Forest",
  total_grass_tco2e   = "Total: Grassland",
  total_wet_tco2e     = "Total: Wetland",
  total_peat_tco2e    = "Total: Peatland",
  # Activity-type aggregates
  total_protect_tco2e = "Total: Protection",
  total_manage_tco2e  = "Total: Management",
  total_restore_tco2e = "Total: Restoration",
  total_ncs_tco2e     = "Total NCS"
)

nbs_individual <- names(pathway_labels)[1:20]
nbs_biome      <- names(pathway_labels)[21:25]
nbs_activity   <- names(pathway_labels)[26:29]
nbs_grand      <- "total_ncs_tco2e"

cat("Countries in Naturebase:", nrow(nbs_raw), "\n")
glimpse(nbs_raw)

wgi_codes <- c(
  CC.EST = "Control of Corruption",
  GE.EST = "Government Effectiveness",
  PV.EST = "Political Stability & No Violence",
  RL.EST = "Rule of Law",
  RQ.EST = "Regulatory Quality",
  VA.EST  = "Voice & Accountability"
)

wgi_raw <- wb_data(
  indicator  = names(wgi_codes),
  start_date = 2018,
  end_date   = 2023,
  return_wide = TRUE
)

# Keep most recent non-NA value per country × indicator
wgi <- wgi_raw %>%
  pivot_longer(cols = all_of(names(wgi_codes)),
               names_to = "indicator", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(iso3c, indicator) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pivot_wider(names_from = indicator, values_from = value) %>%
  select(iso3c, country, date, all_of(names(wgi_codes)))

# Attach World Bank income groups and regions
wb_meta <- wb_countries() %>%
  select(iso3c, income_level, region) %>%
  filter(!income_level %in% c("", "Aggregates"))

wgi <- wgi %>%
  left_join(wb_meta, by = "iso3c") %>%
  mutate(
    income_level = factor(income_level,
      levels = c("Low income", "Lower middle income",
                 "Upper middle income", "High income"))
  )

cat("WGI records downloaded:", nrow(wgi), "\n")

df <- nbs_raw %>%
  rename(iso3c = adm0_id) %>%
  left_join(wgi, by = "iso3c") %>%
  filter(!is.na(CC.EST))  # require at least one WGI estimate

cat("Countries matched:", nrow(df), "\n")
cat("Countries without WGI (dropped):", nrow(nbs_raw) - nrow(df), "\n")

# Compute per-hectare NCS density for size-normalised comparisons
df <- df %>%
  mutate(
    across(
      all_of(c(nbs_individual, nbs_biome, nbs_activity, nbs_grand)),
      ~ . / adm0_ha,
      .names = "{.col}_dens"   # tCO2e / ha_country / yr
    )
  )

desc <- df %>%
  select(all_of(c(nbs_grand, nbs_biome, names(wgi_codes)))) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    n        = sum(!is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    median   = median(value, na.rm = TRUE),
    sd       = sd(value, na.rm = TRUE),
    min      = min(value, na.rm = TRUE),
    max      = max(value, na.rm = TRUE)
  ) %>%
  mutate(
    label = coalesce(pathway_labels[name], wgi_codes[name], name)
  ) %>%
  relocate(label, .after = name)

knitr::kable(desc, digits = 2, caption = "Descriptive statistics — NbS potential (tCO2e/yr) and WGI scores")

wgi_vars <- names(wgi_codes)

cor_biome_mat <- df %>%
  select(all_of(c(nbs_biome, nbs_grand, wgi_vars))) %>%
  mutate(across(all_of(c(nbs_biome, nbs_grand)), log1p)) %>%
  cor(method = "spearman", use = "pairwise.complete.obs")

cor_biome_sub <- cor_biome_mat[c(nbs_biome, nbs_grand), wgi_vars]
rownames(cor_biome_sub) <- pathway_labels[rownames(cor_biome_sub)]
colnames(cor_biome_sub) <- unname(wgi_codes)

knitr::kable(round(cor_biome_sub, 3),
             caption = "Spearman rho — biome-level NCS totals × WGI dimensions")

ggcorrplot(
  cor_biome_sub,
  method    = "square",
  lab       = TRUE,
  lab_size  = 3.5,
  colors    = c("#d73027", "white", "#4575b4"),
  title     = "NbS Biome Totals × WGI — Spearman rho",
  ggtheme   = theme_minimal(base_size = 12)
) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

cor_all_mat <- df %>%
  select(all_of(c(nbs_individual, wgi_vars))) %>%
  mutate(across(all_of(nbs_individual), log1p)) %>%
  cor(method = "spearman", use = "pairwise.complete.obs")

cor_all_sub <- cor_all_mat[nbs_individual, wgi_vars]
rownames(cor_all_sub) <- pathway_labels[rownames(cor_all_sub)]
colnames(cor_all_sub) <- unname(wgi_codes)

ggcorrplot(
  cor_all_sub,
  method    = "square",
  lab       = TRUE,
  lab_size  = 3,
  colors    = c("#d73027", "white", "#4575b4"),
  title     = "Individual NbS Pathways × WGI — Spearman rho",
  ggtheme   = theme_minimal(base_size = 10)
) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))



label_thresh <- quantile(df$total_ncs_tco2e, 0.9, na.rm = TRUE)
risk_thresh  <- quantile(df$PV.EST, 0.1, na.rm = TRUE)

df %>%
  filter(!is.na(PV.EST), total_ncs_tco2e > 0) %>%
  mutate(label = if_else(
    total_ncs_tco2e > label_thresh | PV.EST < risk_thresh,
    adm0_name, NA_character_
  )) %>%
  ggplot(aes(x = PV.EST, y = log10(total_ncs_tco2e))) +
  geom_point(aes(color = income_level), alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                  segment.color = "grey60") +
  scale_color_brewer(palette = "RdYlBu", name = "Income Group", na.value = "grey70") +
  labs(
    x     = "Political Stability & Absence of Violence (WGI, std. score)",
    y     = "Total NCS Potential (log₁₀ tCO₂e/yr)",
    title = "Total NCS Potential vs. Political Stability",
    caption = "Naturebase v3 (2025); WGI most-recent 2018–2023"
  ) +
  theme_minimal(base_size = 12)

df %>%
  filter(!is.na(CC.EST), total_forest_tco2e > 0) %>%
  mutate(label = if_else(
    total_forest_tco2e > quantile(total_forest_tco2e, 0.88, na.rm = TRUE),
    adm0_name, NA_character_
  )) %>%
  ggplot(aes(x = CC.EST, y = log10(total_forest_tco2e))) +
  geom_point(aes(color = income_level), alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 18,
                  segment.color = "grey60") +
  scale_color_brewer(palette = "RdYlBu", name = "Income Group", na.value = "grey70") +
  labs(
    x     = "Control of Corruption (WGI, std. score)",
    y     = "Forest NCS Potential (log₁₀ tCO₂e/yr)",
    title = "Forest NCS Potential vs. Control of Corruption"
  ) +
  theme_minimal(base_size = 12)

df %>%
  filter(!is.na(RL.EST), total_wet_tco2e > 0) %>%
  mutate(label = if_else(
    total_wet_tco2e > quantile(total_wet_tco2e, 0.88, na.rm = TRUE),
    adm0_name, NA_character_
  )) %>%
  ggplot(aes(x = RL.EST, y = log10(total_wet_tco2e))) +
  geom_point(aes(color = income_level), alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 18,
                  segment.color = "grey60") +
  scale_color_brewer(palette = "RdYlBu", name = "Income Group", na.value = "grey70") +
  labs(
    x     = "Rule of Law (WGI, std. score)",
    y     = "Wetland NCS Potential (log₁₀ tCO₂e/yr)",
    title = "Wetland NCS Potential vs. Rule of Law"
  ) +
  theme_minimal(base_size = 12)

income_cors <- df %>%
  filter(!is.na(income_level)) %>%
  group_by(income_level) %>%
  summarise(
    n = n(),
    across(
      all_of(wgi_vars),
      ~ cor(log1p(total_ncs_tco2e), .x,
            method = "spearman", use = "complete.obs"),
      .names = "{.col}"
    )
  ) %>%
  pivot_longer(all_of(wgi_vars), names_to = "wgi_var", values_to = "rho") %>%
  mutate(wgi_label = wgi_codes[wgi_var])

ggplot(income_cors,
       aes(x = rho, y = reorder(income_level, rho, FUN = mean),
           color = wgi_label)) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_brewer(palette = "Dark2", name = "WGI Dimension") +
  labs(
    x     = "Spearman rho (log1p Total NCS ~ WGI)",
    y     = NULL,
    title = "NCS–Governance Correlations by Income Group"
  ) +
  theme_minimal(base_size = 12)

df %>%
  filter(!is.na(income_level), total_ncs_tco2e > 0) %>%
  pivot_longer(cols = all_of(wgi_vars),
               names_to = "wgi_var", values_to = "wgi_val") %>%
  mutate(wgi_label = wgi_codes[wgi_var]) %>%
  ggplot(aes(x = wgi_val, y = log1p(total_ncs_tco2e), color = income_level)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(income_level ~ wgi_label,
             labeller = labeller(wgi_label = label_wrap_gen(12))) +
  scale_color_brewer(palette = "RdYlBu") +
  labs(
    x     = "WGI Score",
    y     = "Total NCS (log1p tCO₂e/yr)",
    title = "Total NCS ~ WGI, stratified by income group"
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(size = 7))

region_cors <- df %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(
    n = n(),
    across(
      all_of(wgi_vars),
      ~ cor(log1p(total_ncs_tco2e), .x,
            method = "spearman", use = "complete.obs"),
      .names = "{.col}"
    )
  ) %>%
  pivot_longer(all_of(wgi_vars), names_to = "wgi_var", values_to = "rho") %>%
  mutate(wgi_label = wgi_codes[wgi_var])

ggplot(region_cors,
       aes(x = rho, y = reorder(region, rho, FUN = mean),
           fill = wgi_label)) +
  geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2", name = "WGI Dimension") +
  labs(
    x     = "Spearman rho (log1p Total NCS ~ WGI)",
    y     = NULL,
    title = "Regional Variation in NCS–Governance Correlations"
  ) +
  theme_minimal(base_size = 11)

ols <- lm(
  log1p(total_ncs_tco2e) ~
    CC.EST + GE.EST + PV.EST + RL.EST + RQ.EST + VA.EST,
  data = df
)

summary(ols)

tidy(ols, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = wgi_codes[term]) %>%
  ggplot(aes(x = estimate,
             xmin = conf.low, xmax = conf.high,
             y = reorder(term, estimate))) +
  geom_pointrange(size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x        = "OLS Coefficient (95% CI)",
    y        = NULL,
    title    = "OLS: log1p(Total NCS) ~ WGI",
    subtitle = paste0("Adj. R² = ",
                      round(summary(ols)$adj.r.squared, 3),
                      " | N = ", nobs(ols))
  ) +
  theme_minimal(base_size = 12)

biome_models <- map(
  set_names(c(nbs_biome, nbs_grand)),
  ~ lm(
      reformulate(wgi_vars, response = paste0("log1p(", .x, ")")),
      data = df
    ) %>%
    tidy(conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      outcome   = pathway_labels[.x],
      wgi_label = wgi_codes[term]
    )
) %>%
  bind_rows()

ggplot(biome_models,
       aes(x = estimate, xmin = conf.low, xmax = conf.high,
           y = reorder(wgi_label, estimate),
           color = outcome)) +
  geom_pointrange(position = position_dodge(0.7), size = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_brewer(palette = "Set1", name = "NbS Outcome") +
  facet_wrap(~outcome, nrow = 2) +
  labs(
    x     = "OLS Coefficient (95% CI)",
    y     = NULL,
    title = "Per-biome OLS: log1p(NbS potential) ~ WGI"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))

export_cols <- c(
  "iso3c", "adm0_name", "adm0_ha",
  nbs_biome, nbs_activity, nbs_grand,
  wgi_vars, "income_level", "region"
)

df_export <- df %>%
  select(all_of(export_cols))

write_csv(df_export, here("data/nbs_wgi_merged.csv"))
cat("Exported", nrow(df_export), "countries →",
    here("data/nbs_wgi_merged.csv"), "\n")

sessionInfo()
