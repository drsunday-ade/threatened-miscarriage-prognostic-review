############################################################
# Section 4 – Ultrasound predictors and performance
# Systematic Review of First-Trimester Threatened Miscarriage
############################################################

## ---------------------------
## 1. Setup and read data
## ---------------------------

# install.packages(c("tidyverse", "knitr", "kableExtra", "forcats", "ggplot2"))

library(tidyverse)
library(knitr)
library(kableExtra)
library(forcats)
library(ggplot2)

extract_path <- "extracted_variables.csv"

dat_long <- read_csv(extract_path, show_col_types = FALSE)


## -----------------------------------------
## 2. Annotate each row with study_id
## -----------------------------------------

dat_annot <- dat_long %>%
  mutate(
    study_id_tmp = if_else(VariableName == "study_id", ExtractedValue, NA_character_)
  ) %>%
  tidyr::fill(study_id_tmp, .direction = "down") %>%
  rename(study_id = study_id_tmp)

unique(dat_annot$study_id)


## -----------------------------------------
## 3. Core metadata (labels, N, primary vs evidence-synthesis)
## -----------------------------------------

# 3.1 Study labels (author + year)
study_labels <- dat_annot %>%
  filter(VariableName %in% c("first_author", "publication_year")) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  group_by(study_id, VariableName) %>%
  summarise(ExtractedValue = ExtractedValue[1], .groups = "drop") %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue) %>%
  mutate(
    publication_year = suppressWarnings(as.integer(publication_year)),
    label = if_else(
      !is.na(publication_year),
      paste0(first_author, " (", publication_year, ")"),
      first_author
    )
  )

# 3.2 Minimal wide frame with N_total and study type
vars_for_N <- c(
  "sample_size_total",
  "n_tm",
  "n_total_recruited",
  "n_recruited",
  "n_completed",
  "n_total_women_meta",
  "n_miscarriage",
  "n_miscarriage_28w",
  "n_pregnancy_loss",
  "n_ongoing",
  "n_live_birth"
)

meta_wide <- dat_annot %>%
  filter(VariableName %in% vars_for_N) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  group_by(study_id, VariableName) %>%
  summarise(ExtractedValue = ExtractedValue[1], .groups = "drop") %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue) %>%
  mutate(
    across(
      c(sample_size_total, n_tm,
        n_total_recruited, n_recruited, n_completed,
        n_total_women_meta,
        n_miscarriage, n_miscarriage_28w,
        n_pregnancy_loss, n_ongoing, n_live_birth),
      ~ suppressWarnings(as.numeric(.x))
    ),
    N_total = coalesce(
      n_total_recruited,
      n_recruited,
      n_completed,
      sample_size_total,
      n_tm,
      n_total_women_meta
    )
  ) %>%
  left_join(study_labels %>% select(study_id, label, publication_year),
            by = "study_id")

# Evidence-synthesis / meta-analyses
evidence_ids <- c("Cao_2025", "Pillai_2016", "Pillai_2018", "Mishoe_2020")

meta_wide <- meta_wide %>%
  mutate(
    study_type = if_else(study_id %in% evidence_ids,
                         "Evidence synthesis / meta-analysis",
                         "Primary cohort / model")
  )

meta_primary  <- meta_wide %>% filter(study_type == "Primary cohort / model")
meta_evidence <- meta_wide %>% filter(study_type == "Evidence synthesis / meta-analysis")


############################################################
# PART A – ULTRASOUND FEATURE PRESENCE
############################################################

## 4.1 Identify ultrasound features per study (FHR, CRL, GS, hematoma, yolk sac)

us_features_raw <- dat_annot %>%
  mutate(
    has_FHR =
      (str_detect(Category,   regex("ultrasound", ignore_case = TRUE)) &
         str_detect(VariableName, regex("FHR|heart rate", ignore_case = TRUE))) |
      str_detect(VariableName, regex("FHR|heart rate", ignore_case = TRUE)) |
      str_detect(Notes,        regex("FHR|heart rate", ignore_case = TRUE)),
    
    has_CRL =
      (str_detect(Category,   regex("ultrasound", ignore_case = TRUE)) &
         str_detect(VariableName, regex("CRL|crown[- ]?rump", ignore_case = TRUE))) |
      str_detect(VariableName, regex("CRL|crown[- ]?rump", ignore_case = TRUE)) |
      str_detect(Notes,        regex("CRL|crown[- ]?rump", ignore_case = TRUE)),
    
    has_GS =
      (str_detect(Category,   regex("ultrasound", ignore_case = TRUE)) &
         str_detect(VariableName, regex("gestational sac|GSD|MSD|MASD|MGSD", ignore_case = TRUE))) |
      str_detect(VariableName, regex("gestational sac|GSD|MSD|MASD|MGSD", ignore_case = TRUE)) |
      str_detect(Notes,        regex("gestational sac|GSD|MSD|MASD|MGSD", ignore_case = TRUE)),
    
    has_hematoma =
      (str_detect(Category,   regex("ultrasound", ignore_case = TRUE)) &
         str_detect(VariableName, regex("hematoma", ignore_case = TRUE))) |
      str_detect(VariableName, regex("hematoma", ignore_case = TRUE)) |
      str_detect(Notes,        regex("hematoma", ignore_case = TRUE)),
    
    has_yolk =
      (str_detect(Category,   regex("ultrasound", ignore_case = TRUE)) &
         str_detect(VariableName, regex("yolk sac", ignore_case = TRUE))) |
      str_detect(VariableName, regex("yolk sac", ignore_case = TRUE)) |
      str_detect(Notes,        regex("yolk sac", ignore_case = TRUE))
  ) %>%
  group_by(study_id) %>%
  summarise(
    FHR              = any(has_FHR, na.rm = TRUE),
    CRL              = any(has_CRL, na.rm = TRUE),
    `Gestational sac`= any(has_GS, na.rm = TRUE),
    Hematoma         = any(has_hematoma, na.rm = TRUE),
    `Yolk sac`       = any(has_yolk, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  left_join(study_labels %>% select(study_id, label, publication_year),
            by = "study_id") %>%
  arrange(publication_year)

## 4.2 Table 4A – Ultrasound features assessed across studies

table4A <- us_features_raw %>%
  mutate(
    across(c(FHR, CRL, `Gestational sac`, Hematoma, `Yolk sac`),
           ~ if_else(.x, "Yes", "No"))
  ) %>%
  transmute(
    Study = label,
    FHR,
    CRL,
    `Gestational sac`,
    Hematoma,
    `Yolk sac`
  )

table4A_kable <- kbl(
  table4A,
  caption = "Table 4A. Ultrasound features evaluated across included studies.",
  align   = c("l", "c", "c", "c", "c", "c"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 9
  ) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE)

table4A_kable
# save_kable(table4A_kable, "table4A_ultrasound_features.html")


## 4.3 Figure 4C – Ultrasound feature heatmap

us_long <- us_features_raw %>%
  pivot_longer(
    cols      = c(FHR, CRL, `Gestational sac`, Hematoma, `Yolk sac`),
    names_to  = "Feature",
    values_to = "Has"
  ) %>%
  mutate(
    Has        = if_else(Has, "Yes", "No"),
    StudyLabel = fct_reorder(label, publication_year, .fun = min, .na_rm = TRUE)
  )

fig4C <- ggplot(us_long, aes(x = Feature, y = StudyLabel, fill = Has)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Yes" = "#1a9641", "No" = "#f0f0f0"),
    name   = "Feature evaluated"
  ) +
  labs(
    title = "Figure 4C. Ultrasound features evaluated across included studies",
    x     = "Ultrasound feature",
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0),
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "bottom"
  )

fig4C
# ggsave("figure4C_ultrasound_feature_heatmap.png", fig4C, width = 6.5, height = 4.5, dpi = 300)


############################################################
# PART B – QUANTITATIVE ULTRASOUND PERFORMANCE
############################################################

## 5.1 Pull structured ultrasound performance variables

us_vars_perf <- c(
  # FHR
  "fhr_cutoff_main",
  "fhr_unit",
  "fhr_sensitivity_main_cutoff",
  "fhr_specificity_main_cutoff",
  "auc_fhr",
  "fhr_summary_sensitivity_percent",
  "fhr_summary_specificity_percent",
  "fhr_total_women",
  
  # CRL / embryo growth
  "crl_cutoff_main",
  "crl_unit",
  "crl_sensitivity_main_cutoff",
  "crl_specificity_main_cutoff",
  "auc_crl",
  "crl_summary_sensitivity_percent",
  "crl_summary_specificity_percent",
  "crl_total_women",
  
  # Gestational sac metrics
  "gs_cutoff_main",
  "gs_unit",
  "gs_sensitivity_main_cutoff",
  "gs_specificity_main_cutoff",
  "auc_gs",
  "gs_summary_sensitivity_percent",
  "gs_summary_specificity_percent",
  "gs_total_women"
)

us_raw <- dat_annot %>%
  filter(VariableName %in% us_vars_perf) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  group_by(study_id, VariableName) %>%
  summarise(ExtractedValue = ExtractedValue[1], .groups = "drop") %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue) %>%
  left_join(meta_wide %>% select(study_id, label, N_total, study_type),
            by = "study_id")

# Ensure all potential columns exist (if not, create as NA)
for (col in us_vars_perf) {
  if (!col %in% names(us_raw)) {
    us_raw[[col]] <- NA_character_
  }
}

## 5.2 Derive unified metrics for FHR and CRL/GS

us_meta <- us_raw %>%
  mutate(
    # Classify source for tables
    source_class = case_when(
      study_id %in% c("Mishoe_2020", "Pillai_2018") ~
        "Ultrasound meta-analysis / pooled ED cohorts",
      study_id %in% evidence_ids ~ "Evidence synthesis / meta-analysis",
      TRUE ~ "Primary cohort / model"
    ),
    
    #### FHR ####
    N_fhr = case_when(
      !is.na(fhr_total_women) ~ suppressWarnings(as.numeric(fhr_total_women)),
      TRUE ~ N_total
    ),
    
    sens_fhr = coalesce(
      suppressWarnings(as.numeric(fhr_sensitivity_main_cutoff)),
      suppressWarnings(as.numeric(fhr_summary_sensitivity_percent))
    ),
    
    spec_fhr = coalesce(
      suppressWarnings(as.numeric(fhr_specificity_main_cutoff)),
      suppressWarnings(as.numeric(fhr_summary_specificity_percent))
    ),
    
    auc_fhr_num    = suppressWarnings(as.numeric(auc_fhr)),
    cutoff_fhr_num = suppressWarnings(as.numeric(fhr_cutoff_main)),
    
    #### CRL / GS ####
    N_crl = case_when(
      !is.na(crl_total_women) ~ suppressWarnings(as.numeric(crl_total_women)),
      !is.na(gs_total_women)  ~ suppressWarnings(as.numeric(gs_total_women)),
      TRUE ~ N_total
    ),
    
    sens_crl = coalesce(
      suppressWarnings(as.numeric(crl_sensitivity_main_cutoff)),
      suppressWarnings(as.numeric(crl_summary_sensitivity_percent)),
      suppressWarnings(as.numeric(gs_sensitivity_main_cutoff)),
      suppressWarnings(as.numeric(gs_summary_sensitivity_percent))
    ),
    
    spec_crl = coalesce(
      suppressWarnings(as.numeric(crl_specificity_main_cutoff)),
      suppressWarnings(as.numeric(crl_summary_specificity_percent)),
      suppressWarnings(as.numeric(gs_specificity_main_cutoff)),
      suppressWarnings(as.numeric(gs_summary_specificity_percent))
    ),
    
    auc_crl_num = coalesce(
      suppressWarnings(as.numeric(auc_crl)),
      suppressWarnings(as.numeric(auc_gs))
    ),
    
    cutoff_crl_num = coalesce(
      suppressWarnings(as.numeric(crl_cutoff_main)),
      suppressWarnings(as.numeric(gs_cutoff_main))
    )
  )
# NOTE: 'label' is already in us_raw via meta_wide join, so we do NOT re-join study_labels


## 5.3 Table 4B – FHR predictor performance (robust to empty data)

table4B_fhr <- us_meta %>%
  filter(!is.na(sens_fhr) | !is.na(spec_fhr)) %>%
  transmute(
    Feature          = "Fetal heart rate",
    Study            = label,
    `Study type`     = source_class,
    `N women`        = N_fhr,
    `Cut-off (bpm)`  = if_else(
      is.na(cutoff_fhr_num), "NR", as.character(cutoff_fhr_num)
    ),
    `Sensitivity (%)`= if_else(
      is.na(sens_fhr), NA_character_, sprintf("%.1f", sens_fhr)
    ),
    `Specificity (%)`= if_else(
      is.na(spec_fhr), NA_character_, sprintf("%.1f", spec_fhr)
    ),
    `AUC`            = if_else(
      is.na(auc_fhr_num), NA_character_, sprintf("%.2f", auc_fhr_num)
    )
  )

# If no FHR accuracy data, create a single informative placeholder row
if (nrow(table4B_fhr) == 0) {
  table4B_fhr <- tibble(
    Feature           = "Fetal heart rate",
    Study             = "No study reported extractable FHR accuracy metrics",
    `Study type`      = NA_character_,
    `N women`         = NA_real_,
    `Cut-off (bpm)`   = "NR",
    `Sensitivity (%)` = NA_character_,
    `Specificity (%)` = NA_character_,
    `AUC`             = NA_character_
  )
}

table4B_fhr <- table4B_fhr %>%
  arrange(`Study type`, Study)

table4B_fhr_kable <- kbl(
  table4B_fhr,
  caption = "Table 4B. Prognostic performance of fetal heart rate for miscarriage after threatened miscarriage.",
  align   = c("l", "l", "l", "c", "c", "c", "c", "c"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 9
  ) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE)

table4B_fhr_kable
# save_kable(table4B_fhr_kable, "table4B_fhr_performance.html")


## 5.4 Table 4C – CRL / gestational sac predictor performance (robust to empty data)

table4C_crl <- us_meta %>%
  filter(!is.na(sens_crl) | !is.na(spec_crl)) %>%
  transmute(
    Feature = "CRL / embryo size or sac–CRL discordance",
    Study   = label,
    `Study type` = source_class,
    `N women` = N_crl,
    `Cut-off (CRL or sac metric)` = if_else(
      is.na(cutoff_crl_num), "NR", as.character(cutoff_crl_num)
    ),
    `Sensitivity (%)` = if_else(
      is.na(sens_crl), NA_character_, sprintf("%.1f", sens_crl)
    ),
    `Specificity (%)` = if_else(
      is.na(spec_crl), NA_character_, sprintf("%.1f", spec_crl)
    ),
    `AUC` = if_else(
      is.na(auc_crl_num), NA_character_, sprintf("%.2f", auc_crl_num)
    )
  )

if (nrow(table4C_crl) == 0) {
  table4C_crl <- tibble(
    Feature = "CRL / embryo size or sac–CRL discordance",
    Study   = "No study reported extractable CRL / sac accuracy metrics",
    `Study type` = NA_character_,
    `N women` = NA_real_,
    `Cut-off (CRL or sac metric)` = "NR",
    `Sensitivity (%)` = NA_character_,
    `Specificity (%)` = NA_character_,
    `AUC` = NA_character_
  )
}

table4C_crl <- table4C_crl %>%
  arrange(`Study type`, Study)

table4C_crl_kable <- kbl(
  table4C_crl,
  caption = "Table 4C. Prognostic performance of CRL / gestational sac parameters for miscarriage after threatened miscarriage.",
  align   = c("l", "l", "l", "c", "l", "c", "c", "c"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 9
  ) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE)

table4C_crl_kable
# save_kable(table4C_crl_kable, "table4C_crl_gs_performance.html")


############################################################
# 6. FIGURE 4A – FHR sensitivity & specificity
############################################################

fhr_plot_data <- us_meta %>%
  filter(!is.na(sens_fhr) | !is.na(spec_fhr)) %>%
  transmute(
    Study = label,
    Type  = source_class,
    sens  = sens_fhr,
    spec  = spec_fhr
  ) %>%
  pivot_longer(cols = c(sens, spec),
               names_to = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(
    Metric = recode(Metric,
                    sens = "Sensitivity",
                    spec = "Specificity"),
    Study  = fct_reorder(Study, Value, .fun = mean, .na_rm = TRUE)
  )

if (nrow(fhr_plot_data) > 0) {
  fig4A <- ggplot(fhr_plot_data,
                  aes(x = Value, y = Study,
                      color = Metric, shape = Type)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "grey80") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    labs(
      title = "Figure 4A. Sensitivity and specificity of fetal heart rate cut-offs",
      x     = "Percentage (%)",
      y     = NULL,
      color = NULL,
      shape = "Study type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0),
      axis.text.y      = element_text(size = 9),
      axis.text.x      = element_text(size = 9),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
  
  fig4A
  # ggsave("figure4A_fhr_sens_spec.png", fig4A, width = 6.5, height = 4.5, dpi = 300)
} else {
  message("Figure 4A: No FHR accuracy data available to plot.")
}


############################################################
# 7. FIGURE 4B – CRL / GS sensitivity & specificity
############################################################

crl_plot_data <- us_meta %>%
  filter(!is.na(sens_crl) | !is.na(spec_crl)) %>%
  transmute(
    Study = label,
    Type  = source_class,
    sens  = sens_crl,
    spec  = spec_crl
  ) %>%
  pivot_longer(cols = c(sens, spec),
               names_to = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(
    Metric = recode(Metric,
                    sens = "Sensitivity",
                    spec = "Specificity"),
    Study  = fct_reorder(Study, Value, .fun = mean, .na_rm = TRUE)
  )

if (nrow(crl_plot_data) > 0) {
  fig4B <- ggplot(crl_plot_data,
                  aes(x = Value, y = Study,
                      color = Metric, shape = Type)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "grey80") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    labs(
      title = "Figure 4B. Sensitivity and specificity of CRL / gestational sac predictors",
      x     = "Percentage (%)",
      y     = NULL,
      color = NULL,
      shape = "Study type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0),
      axis.text.y      = element_text(size = 9),
      axis.text.x      = element_text(size = 9),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
  
  fig4B
  # ggsave("figure4B_crl_gs_sens_spec.png", fig4B, width = 6.5, height = 4.5, dpi = 300)
} else {
  message("Figure 4B: No CRL / GS accuracy data available to plot.")
}

############################################################
# End of Section 4 script
############################################################


############################################################
# OUTPUT FOLDER: SAVE ALL TABLES & FIGURES
############################################################

# Create output directory for Section 4
out_dir <- file.path("outputs", "section4")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Save tables as HTML
if (exists("table4A_kable")) {
  save_kable(table4A_kable,
             file.path(out_dir, "table4A_ultrasound_features.html"))
}
if (exists("table4B_fhr_kable")) {
  save_kable(table4B_fhr_kable,
             file.path(out_dir, "table4B_fhr_performance.html"))
}
if (exists("table4C_crl_kable")) {
  save_kable(table4C_crl_kable,
             file.path(out_dir, "table4C_crl_gs_performance.html"))
}

# Save figures as PNG (only if they were created)
if (exists("fig4C") && inherits(fig4C, "ggplot")) {
  ggsave(file.path(out_dir, "figure4C_ultrasound_feature_heatmap.png"),
         fig4C, width = 6.5, height = 4.5, dpi = 300)
}
if (exists("fig4A") && inherits(fig4A, "ggplot")) {
  ggsave(file.path(out_dir, "figure4A_fhr_sens_spec.png"),
         fig4A, width = 6.5, height = 4.5, dpi = 300)
}
if (exists("fig4B") && inherits(fig4B, "ggplot")) {
  ggsave(file.path(out_dir, "figure4B_crl_gs_sens_spec.png"),
         fig4B, width = 6.5, height = 4.5, dpi = 300)
}
