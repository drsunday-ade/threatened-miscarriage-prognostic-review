############################################################
# Section 5 – Multivariable and Machine-Learning Models
# Systematic Review of First-Trimester Threatened Miscarriage
############################################################

## ---------------------------
## 0. Output folder
## ---------------------------

output_dir <- "section5_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

## -----------------------------------------
## 3. Study labels and N per study
## -----------------------------------------

# 3.1 First author + year
study_labels <- dat_annot %>%
  filter(VariableName %in% c("first_author", "publication_year")) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  group_by(study_id, VariableName) %>%
  summarise(ExtractedValue = ExtractedValue[1], .groups = "drop") %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue) %>%
  mutate(
    publication_year = suppressWarnings(as.integer(publication_year)),
    StudyLabel = if_else(
      !is.na(publication_year),
      paste0(first_author, " (", publication_year, ")"),
      first_author
    )
  )

# 3.2 N_total per study
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
  left_join(
    study_labels %>% select(study_id, StudyLabel, publication_year),
    by = "study_id"
  )

############################################################
# 4. Build model-level dataset (Ku 2015, Ku 2025, Sammut 2025)
############################################################

## 4.1 Ku 2015 – clinical + biochemical + US logistic model

ku2015_details <- dat_annot %>%
  filter(
    study_id == "Ku_2015",
    VariableName %in% c(
      "model_name",
      "model_type",
      "model_purpose",
      "predictors_in_model",
      "internal_validation_method",
      "external_validation"
    )
  ) %>%
  select(VariableName, ExtractedValue) %>%
  pivot_wider(names_from = VariableName, values_from = ExtractedValue)

ku2015_auc <- dat_annot %>%
  filter(study_id == "Ku_2015", VariableName == "auc_roc") %>%
  summarise(
    auc = suppressWarnings(as.numeric(ExtractedValue[1]))
  ) %>%
  pull(auc)

ku2015_model <- ku2015_details %>%
  mutate(
    study_id    = "Ku_2015",
    model_id    = if_else(!is.na(model_name), model_name, "Ku_2015_risk_model"),
    model_class = "Multivariable logistic regression",
    auc         = ku2015_auc,
    auc_ci_text = NA_character_
  ) %>%
  select(
    study_id,
    model_id,
    model_class,
    model_type,
    model_purpose,
    predictors_in_model,
    internal_validation_method,
    external_validation,
    auc,
    auc_ci_text
  )

## 4.2 Ku 2025 – holistic models (Models 1–5, with AUCs)

ku2025_models <- purrr::map_dfr(1:5, function(k) {
  desc_row <- dat_annot %>%
    filter(study_id == "Ku_2025",
           VariableName == paste0("model", k, "_description")) %>%
    slice(1)
  
  auc_row <- dat_annot %>%
    filter(study_id == "Ku_2025",
           VariableName == paste0("model", k, "_auroc")) %>%
    slice(1)
  
  if (nrow(desc_row) == 0 & nrow(auc_row) == 0) {
    return(NULL)
  }
  
  desc_val <- if (nrow(desc_row) > 0) desc_row$ExtractedValue[1] else NA_character_
  auc_char <- if (nrow(auc_row) > 0) auc_row$ExtractedValue[1] else NA_character_
  auc_num  <- suppressWarnings(as.numeric(auc_char))
  auc_ci   <- if (nrow(auc_row) > 0) auc_row$Notes[1] else NA_character_
  
  tibble(
    study_id    = "Ku_2025",
    model_id    = paste0("Model ", k),
    model_class = "Multivariable logistic regression",
    model_type  = "Logistic regression",
    model_purpose = "Clinical + ultrasound + biochemical prediction model",
    predictors_in_model = desc_val,
    internal_validation_method = NA_character_,
    external_validation        = "No (single cohort)",
    auc         = auc_num,
    auc_ci_text = auc_ci
  )
})

## 4.3 Sammut 2025 – MLR vs Random Forest

# MLR (multivariable logistic regression)
sammut_mlr_auc <- dat_annot %>%
  filter(study_id == "Sammut_2025", VariableName == "mlr_auc") %>%
  summarise(
    auc = suppressWarnings(as.numeric(ExtractedValue[1]))
  ) %>%
  pull(auc)

sammut_mlr_pred <- dat_annot %>%
  filter(study_id == "Sammut_2025", VariableName == "mlr_final_predictors") %>%
  summarise(
    predictors = ExtractedValue[1]
  ) %>%
  pull(predictors)

# Random forest (testing AUC)
sammut_rf_auc <- dat_annot %>%
  filter(study_id == "Sammut_2025", VariableName == "rf_auc_testing") %>%
  summarise(
    auc = suppressWarnings(as.numeric(ExtractedValue[1]))
  ) %>%
  pull(auc)

sammut_models <- tibble(
  study_id    = c("Sammut_2025",        "Sammut_2025"),
  model_id    = c("MLR risk model",     "Random forest model"),
  model_class = c("Multivariable logistic regression",
                  "Machine learning (random forest)"),
  model_type  = c("Logistic regression", "Random forest"),
  model_purpose = c(
    "Multivariable clinical + ultrasound + biochemical model",
    "Non-linear machine-learning risk model"
  ),
  predictors_in_model = c(sammut_mlr_pred, sammut_mlr_pred),
  internal_validation_method = c("Train/test split", "Train/test split"),
  external_validation        = c("No (single centre)", "No (single centre)"),
  auc         = c(sammut_mlr_auc, sammut_rf_auc),
  auc_ci_text = c(NA_character_, NA_character_)
)

## 4.4 Combine all models and join N_total + StudyLabel

models_all <- bind_rows(
  ku2015_model,
  ku2025_models,
  sammut_models
) %>%
  left_join(
    meta_wide %>% select(study_id, N_total, StudyLabel, publication_year),
    by = "study_id"
  ) %>%
  mutate(
    ModelLabel = paste0(StudyLabel, " – ", model_id),
    auc_label = case_when(
      !is.na(auc) & !is.na(auc_ci_text) ~
        paste0(sprintf("%.2f", auc), " (", auc_ci_text, ")"),
      !is.na(auc) ~ sprintf("%.2f", auc),
      TRUE        ~ "NR"
    )
  )

############################################################
# PART A – Table 5A: Multivariable & ML models
############################################################

table5A <- models_all %>%
  transmute(
    Study              = StudyLabel,
    `Model`            = model_id,
    `Model class`      = model_class,
    `N (derivation sample)` =
      if_else(is.na(N_total), "NR", as.character(N_total)),
    `AUC / c-statistic` = auc_label,
    `Predictors included` = predictors_in_model,
    `Internal validation` =
      coalesce(internal_validation_method, "Not reported"),
    `External validation` =
      coalesce(external_validation, "Not reported")
  ) %>%
  arrange(Study, `Model`)

table5A_kable <- kbl(
  table5A,
  caption = "Table 5A. Multivariable and machine-learning models predicting miscarriage after threatened miscarriage.",
  align   = c("l", "l", "l", "c", "c", "l", "l", "l"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 9
  ) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE)

table5A_kable
# Original (left untouched):
# save_kable(table5A_kable, "table5A_multivariable_models.html")

# NEW: save to output folder
save_kable(
  table5A_kable,
  file.path(output_dir, "table5A_multivariable_models.html")
)

############################################################
# PART B – Figure 5A: AUC dot plot for all models
############################################################

models_plot <- models_all %>%
  filter(!is.na(auc)) %>%
  mutate(
    model_class = fct_relevel(
      model_class,
      "Multivariable logistic regression",
      "Machine learning (random forest)"
    ),
    ModelLabel_f = fct_reorder(
      ModelLabel,
      auc,
      .fun = function(x) ifelse(is.na(x), 0, x)
    )
  )

fig5A <- ggplot(models_plot,
                aes(x = auc, y = ModelLabel_f, colour = model_class)) +
  geom_vline(xintercept = 0.5,
             linetype  = "dashed",
             colour    = "grey70",
             linewidth = 0.5) +
  geom_point(size = 3) +
  scale_x_continuous(
    limits = c(0.5, 1.0),
    breaks = seq(0.5, 1.0, by = 0.1)
  ) +
  labs(
    title  = "Figure 5A. Discriminative performance (AUC) of multivariable models",
    x      = "Area under the ROC curve (AUC)",
    y      = NULL,
    colour = "Model class"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", hjust = 0),
    axis.text.y    = element_text(size = 9),
    legend.position = "bottom"
  )

fig5A
# Original:
# ggsave("figure5A_multivariable_model_auc.png", fig5A,
#        width = 7, height = 4.5, dpi = 300)

# NEW: save to output folder
ggsave(
  filename = file.path(output_dir, "figure5A_multivariable_model_auc.png"),
  plot     = fig5A,
  width    = 7,
  height   = 4.5,
  dpi      = 300
)

############################################################
# PART C – Figure 5B: Logistic vs Random Forest (Sammut 2025)
############################################################

head2head <- models_all %>%
  filter(
    study_id == "Sammut_2025",
    !is.na(auc)
  ) %>%
  mutate(
    model_class = fct_relevel(
      model_class,
      "Multivariable logistic regression",
      "Machine learning (random forest)"
    )
  )

fig5B <- ggplot(head2head,
                aes(x = auc,
                    y = StudyLabel,
                    colour = model_class)) +
  geom_vline(xintercept = 0.5,
             linetype  = "dashed",
             colour    = "grey80",
             linewidth = 0.4) +
  geom_line(aes(group = StudyLabel),
            colour    = "grey70",
            linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(
    limits = c(0.5, 1.0),
    breaks = seq(0.5, 1.0, by = 0.1)
  ) +
  labs(
    title  = "Figure 5B. Logistic regression vs random forest (Sammut et al. 2025)",
    x      = "Area under the ROC curve (AUC)",
    y      = NULL,
    colour = "Model class"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", hjust = 0),
    axis.text.y    = element_text(size = 9),
    legend.position = "bottom"
  )

fig5B
# Original:
# ggsave("figure5B_logistic_vs_rf.png", fig5B,
#        width = 6.5, height = 3.5, dpi = 300)

# NEW: save to output folder
ggsave(
  filename = file.path(output_dir, "figure5B_logistic_vs_rf.png"),
  plot     = fig5B,
  width    = 6.5,
  height   = 3.5,
  dpi      = 300
)

############################################################
# PART D – Figure 5C: Heatmap of predictor domains per model
############################################################

domains_long <- models_all %>%
  mutate(
    predictors_text = tolower(coalesce(predictors_in_model, "")),
    has_clinical =
      str_detect(predictors_text, "age")     |
      str_detect(predictors_text, "bmi")     |
      str_detect(predictors_text, "parity")  |
      str_detect(predictors_text, "bleeding")|
      str_detect(predictors_text, "nausea")  |
      str_detect(predictors_text, "history"),
    has_ultrasound =
      str_detect(predictors_text, "crl")            |
      str_detect(predictors_text, "crown")          |
      str_detect(predictors_text, "gestational sac")|
      str_detect(predictors_text, "gsd")            |
      str_detect(predictors_text, "msd")            |
      str_detect(predictors_text, "fhr")            |
      str_detect(predictors_text, "heart rate")     |
      str_detect(predictors_text, "yolk")           |
      str_detect(predictors_text, "hematoma"),
    has_biomarker =
      str_detect(predictors_text, "progesterone") |
      str_detect(predictors_text, "ca-125")      |
      str_detect(predictors_text, "ca125")       |
      str_detect(predictors_text, "hcg")         |
      str_detect(predictors_text, "β-hcg")       |
      str_detect(predictors_text, "pibf")        |
      str_detect(predictors_text, "plgf")        |
      str_detect(predictors_text, "papp-a")
  ) %>%
  arrange(publication_year, StudyLabel, model_id) %>%
  mutate(
    ModelLabel   = paste0(StudyLabel, " – ", model_id),
    ModelLabel_f = factor(ModelLabel, levels = unique(ModelLabel))
  ) %>%
  select(ModelLabel_f, has_clinical, has_ultrasound, has_biomarker) %>%
  pivot_longer(
    cols      = starts_with("has_"),
    names_to  = "Domain",
    values_to = "Included"
  ) %>%
  mutate(
    Domain = recode(
      Domain,
      "has_clinical"   = "Clinical predictors",
      "has_ultrasound" = "Ultrasound predictors",
      "has_biomarker"  = "Biochemical predictors"
    ),
    Included = if_else(Included, "Yes", "No")
  )

fig5C <- ggplot(domains_long,
                aes(x = Domain, y = ModelLabel_f, fill = Included)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Yes" = "#2166ac", "No" = "#f0f0f0"),
    name   = "Predictor domain"
  ) +
  labs(
    title = "Figure 5C. Predictor domains in multivariable models",
    x     = "Predictor domain",
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title     = element_text(face = "bold", hjust = 0),
    axis.text.y    = element_text(size = 8),
    axis.text.x    = element_text(size = 9),
    panel.grid     = element_blank(),
    legend.position = "bottom"
  )

fig5C
# Original:
# ggsave("figure5C_predictor_domain_heatmap.png", fig5C,
#        width = 6.5, height = 4.5, dpi = 300)

# NEW: save to output folder
ggsave(
  filename = file.path(output_dir, "figure5C_predictor_domain_heatmap.png"),
  plot     = fig5C,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)
