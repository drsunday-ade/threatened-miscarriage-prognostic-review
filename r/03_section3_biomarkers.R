############################################################
# Section 3 – Biochemical marker performance
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

# Define output folder for tables & figures
output_dir <- "section3_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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
## 3. Build core metadata (labels, N, type)
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

# 3.2 Minimal wide frame with N_total and type (primary vs evidence-synthesis)
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


############################################################
# 4. PROGESTERONE – PERFORMANCE SUMMARY
############################################################

## 4.1 Raw progesterone-related fields (cohort + meta-analysis)
prog_raw <- dat_annot %>%
  filter(VariableName %in% c(
    "progesterone_cutoff_main",
    "progesterone_sensitivity_main_cutoff",
    "progesterone_specificity_main_cutoff",
    "progesterone_unit",
    "auc_progesterone",
    "pilot_sensitivity_cutoff35",
    "pilot_specificity_cutoff35",
    "validation_sensitivity_cutoff35",
    "validation_specificity_cutoff35",
    "auc_pilot",
    "auc_validation",
    "progesterone_summary_sensitivity_percent",
    "progesterone_summary_sensitivity_95ci",
    "progesterone_summary_specificity_percent",
    "progesterone_summary_specificity_95ci",
    "progesterone_total_women"
  )) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  group_by(study_id, VariableName) %>%
  summarise(ExtractedValue = ExtractedValue[1], .groups = "drop") %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue) %>%
  left_join(meta_wide %>% select(study_id, label, N_total, study_type),
            by = "study_id")

## 4.2 Derive unified metrics (percent sensitivity/specificity, AUC, cutoff)
prog_meta <- prog_raw %>%
  mutate(
    marker = "Progesterone",
    
    # Study classification for the table
    source_class = case_when(
      study_id %in% c("AlMohamady_2016", "Ku_2015", "Ku_2025",
                      "Lek_2017", "Sammut_2025", "Shaamash_2020") ~
        "Primary cohort / model",
      study_id == "Pillai_2016" ~ "Meta-analysis of multiple cohorts",
      TRUE ~ study_type
    ),
    
    # N for the progesterone analyses
    N_marker = case_when(
      !is.na(N_total) ~ N_total,
      !is.na(progesterone_total_women) ~ suppressWarnings(as.numeric(progesterone_total_women)),
      TRUE ~ NA_real_
    ),
    
    # Sensitivity
    sens_prog_main = suppressWarnings(as.numeric(progesterone_sensitivity_main_cutoff)),
    sens_valid35   = suppressWarnings(as.numeric(validation_sensitivity_cutoff35) * 100),
    sens_pilot35   = suppressWarnings(as.numeric(pilot_sensitivity_cutoff35) * 100),
    sens_summary   = suppressWarnings(as.numeric(progesterone_summary_sensitivity_percent)),
    sens_percent   = coalesce(sens_prog_main, sens_valid35, sens_pilot35, sens_summary),
    
    # Specificity
    spec_prog_main = suppressWarnings(as.numeric(progesterone_specificity_main_cutoff)),
    spec_valid35   = suppressWarnings(as.numeric(validation_specificity_cutoff35) * 100),
    spec_pilot35   = suppressWarnings(as.numeric(pilot_specificity_cutoff35) * 100),
    spec_summary   = suppressWarnings(as.numeric(progesterone_summary_specificity_percent)),
    spec_percent   = coalesce(spec_prog_main, spec_valid35, spec_pilot35, spec_summary),
    
    # AUC
    auc_main      = suppressWarnings(as.numeric(auc_progesterone)),
    auc_pilot_num = suppressWarnings(as.numeric(auc_pilot)),
    auc_valid_num = suppressWarnings(as.numeric(auc_validation)),
    auc_value     = coalesce(auc_main, auc_valid_num, auc_pilot_num),
    
    # Cut-off (nmol/L or assay unit)
    cutoff_prog   = suppressWarnings(as.numeric(progesterone_cutoff_main)),
    cutoff_final  = case_when(
      !is.na(cutoff_prog) ~ cutoff_prog,
      study_id == "Lek_2017" ~ 35,  # validation threshold 35 nmol/L
      TRUE ~ NA_real_
    )
  )

## 4.3 Table 3A – Progesterone performance
table3A_prog <- prog_meta %>%
  filter(!is.na(sens_percent) | !is.na(spec_percent)) %>%
  transmute(
    Marker = marker,
    Study  = label,
    `Study type` = source_class,
    `N women` = N_marker,
    `Cut-off` = case_when(
      !is.na(cutoff_final) & !is.na(progesterone_unit) ~
        paste0(cutoff_final, " ", progesterone_unit),
      !is.na(cutoff_final) & is.na(progesterone_unit)  ~
        as.character(cutoff_final),
      study_id == "Pillai_2016" ~ "Various thresholds across studies",
      TRUE ~ "NR"
    ),
    `Sensitivity (%)` = if_else(is.na(sens_percent), NA_character_,
                                sprintf("%.1f", sens_percent)),
    `Specificity (%)` = if_else(is.na(spec_percent), NA_character_,
                                sprintf("%.1f", spec_percent)),
    `AUC` = if_else(is.na(auc_value), NA_character_,
                    sprintf("%.2f", auc_value))
  ) %>%
  arrange(Marker, `Study type`, Study)

table3A_kable <- kbl(
  table3A_prog,
  caption = "Table 3A. Prognostic performance of serum progesterone for miscarriage after threatened miscarriage.",
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

table3A_kable

# Save Table 3A into output folder
save_kable(
  table3A_kable,
  file = file.path(output_dir, "table3A_progesterone_performance.html")
)


############################################################
# 5. CA-125 – PERFORMANCE SUMMARY
############################################################

## 5.1 Raw CA-125 fields (cohort + meta-analyses)
ca_raw <- dat_annot %>%
  filter(VariableName %in% c(
    "ca125_cutoff_main",
    "ca125_sensitivity_main_cutoff",
    "ca125_specificity_main_cutoff",
    "ca125_unit",
    "auc_ca125",
    "ca125_summary_sensitivity_percent",
    "ca125_summary_sensitivity_95ci",
    "ca125_summary_specificity_percent",
    "ca125_summary_specificity_95ci",
    "ca125_total_women",
    "pooled_sensitivity",
    "pooled_sensitivity_95ci",
    "pooled_specificity",
    "pooled_specificity_95ci",
    "pooled_auc_sroc",
    "pooled_auc_95ci"
  )) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  group_by(study_id, VariableName) %>%
  summarise(ExtractedValue = ExtractedValue[1], .groups = "drop") %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue) %>%
  left_join(meta_wide %>% select(study_id, label, N_total, study_type),
            by = "study_id")

## 5.2 Derive unified metrics
ca_meta <- ca_raw %>%
  mutate(
    marker = "CA-125",
    
    source_class = case_when(
      study_id == "AlMohamady_2016" ~ "Primary cohort",
      study_id %in% c("Cao_2025", "Pillai_2016") ~
        "Meta-analysis of multiple cohorts",
      TRUE ~ study_type
    ),
    
    # N for CA-125 analyses
    N_marker = case_when(
      !is.na(N_total) ~ N_total,
      !is.na(ca125_total_women) ~ suppressWarnings(as.numeric(ca125_total_women)),
      TRUE ~ NA_real_
    ),
    
    # Sensitivity
    sens_main     = suppressWarnings(as.numeric(ca125_sensitivity_main_cutoff)),
    sens_summary  = suppressWarnings(as.numeric(ca125_summary_sensitivity_percent)),
    sens_pooled   = suppressWarnings(as.numeric(pooled_sensitivity) * 100),
    sens_percent  = coalesce(sens_main, sens_summary, sens_pooled),
    
    # Specificity
    spec_main     = suppressWarnings(as.numeric(ca125_specificity_main_cutoff)),
    spec_summary  = suppressWarnings(as.numeric(ca125_summary_specificity_percent)),
    spec_pooled   = suppressWarnings(as.numeric(pooled_specificity) * 100),
    spec_percent  = coalesce(spec_main, spec_summary, spec_pooled),
    
    # AUC
    auc_ca_num    = suppressWarnings(as.numeric(auc_ca125)),
    auc_pooled    = suppressWarnings(as.numeric(pooled_auc_sroc)),
    auc_value     = coalesce(auc_ca_num, auc_pooled),
    
    # Cut-off
    cutoff_ca     = suppressWarnings(as.numeric(ca125_cutoff_main))
  )

## 5.3 Table 3B – CA-125 performance
table3B_ca <- ca_meta %>%
  filter(!is.na(sens_percent) | !is.na(spec_percent)) %>%
  transmute(
    Marker = marker,
    Study  = label,
    `Study type` = source_class,
    `N women` = N_marker,
    `Cut-off` = case_when(
      !is.na(cutoff_ca) & !is.na(ca125_unit) ~
        paste0(cutoff_ca, " ", ca125_unit),
      !is.na(cutoff_ca) & is.na(ca125_unit) ~
        as.character(cutoff_ca),
      study_id %in% c("Pillai_2016", "Cao_2025") ~
        "Various thresholds across studies",
      TRUE ~ "NR"
    ),
    `Sensitivity (%)` = if_else(is.na(sens_percent), NA_character_,
                                sprintf("%.1f", sens_percent)),
    `Specificity (%)` = if_else(is.na(spec_percent), NA_character_,
                                sprintf("%.1f", spec_percent)),
    `AUC` = if_else(is.na(auc_value), NA_character_,
                    sprintf("%.2f", auc_value))
  ) %>%
  arrange(Marker, `Study type`, Study)

table3B_kable <- kbl(
  table3B_ca,
  caption = "Table 3B. Prognostic performance of serum CA-125 for miscarriage after threatened miscarriage.",
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

table3B_kable

# Save Table 3B into output folder
save_kable(
  table3B_kable,
  file = file.path(output_dir, "table3B_ca125_performance.html")
)


############################################################
# 6. FIGURE 3A – Progesterone sensitivity & specificity
############################################################

prog_plot_data <- prog_meta %>%
  filter(!is.na(sens_percent) | !is.na(spec_percent)) %>%
  transmute(
    Study = label,
    Type  = source_class,
    sens  = sens_percent,
    spec  = spec_percent
  ) %>%
  pivot_longer(cols = c(sens, spec),
               names_to = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(
    Metric = recode(Metric,
                    sens = "Sensitivity",
                    spec = "Specificity"),
    Study = fct_reorder(Study, Value, .fun = mean, .na_rm = TRUE)
  )

fig3A <- ggplot(prog_plot_data,
                aes(x = Value, y = Study,
                    color = Metric, shape = Type)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey80") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "Figure 3A. Sensitivity and specificity of serum progesterone cut-offs",
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

fig3A

# Save Figure 3A into output folder
ggsave(
  filename = file.path(output_dir, "figure3A_progesterone_sens_spec.png"),
  plot     = fig3A,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)


############################################################
# 7. FIGURE 3B – CA-125 sensitivity & specificity
############################################################

ca_plot_data <- ca_meta %>%
  filter(!is.na(sens_percent) | !is.na(spec_percent)) %>%
  transmute(
    Study = label,
    Type  = source_class,
    sens  = sens_percent,
    spec  = spec_percent
  ) %>%
  pivot_longer(
    cols      = c(sens, spec),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(
    Metric = recode(Metric,
                    sens = "Sensitivity",
                    spec = "Specificity"),
    Study = fct_reorder(Study, Value, .fun = mean, .na_rm = TRUE)
  )

fig3B <- ggplot(ca_plot_data,
                aes(x = Value, y = Study,
                    color = Metric, shape = Type)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey80") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "Figure 3B. Sensitivity and specificity of serum CA-125",
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

fig3B

# Save Figure 3B into output folder
ggsave(
  filename = file.path(output_dir, "figure3B_ca125_sens_spec.png"),
  plot     = fig3B,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)

############################################################
# End of Section 3 script
############################################################
