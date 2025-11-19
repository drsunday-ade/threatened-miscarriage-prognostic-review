############################################################
# Section 6 – Risk of Bias and Publication Bias
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

# ---- NEW: create output folder for Section 6 ----
out_dir <- "section6_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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
## 3. Study labels and N (for context)
## -----------------------------------------

# 3.1 First author + year -> StudyLabel
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

# 3.2 N_total per study (same logic as previous sections)
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
# 4. Risk of bias – tools and domain-level concerns
############################################################

## 4.1 Pull risk-of-bias tool, domains, and overall judgement

rob_tools <- dat_annot %>%
  filter(VariableName == "risk_of_bias_tool") %>%
  group_by(study_id) %>%
  summarise(
    risk_of_bias_tool = ExtractedValue[1],
    rob_tool_notes    = Notes[1],
    .groups = "drop"
  )

rob_domains_free_text <- dat_annot %>%
  filter(VariableName == "risk_of_bias_domains") %>%
  group_by(study_id) %>%
  summarise(
    risk_of_bias_domains = ExtractedValue[1],
    .groups = "drop"
  )

rob_overall <- dat_annot %>%
  filter(VariableName == "risk_of_bias_overall") %>%
  group_by(study_id) %>%
  summarise(
    risk_of_bias_overall = ExtractedValue[1],
    rob_overall_notes    = Notes[1],
    .groups = "drop"
  )

## 4.2 Domain-specific textual concerns (patient selection, index test, etc.)

rob_domains_detail <- dat_annot %>%
  filter(VariableName %in% c(
    "common_sources_of_bias",
    "common_bias_patient_selection",
    "common_bias_index_test",
    "common_bias_reference_standard",
    "common_bias_flow_timing"
  )) %>%
  select(study_id, VariableName, ExtractedValue) %>%
  pivot_wider(
    names_from  = VariableName,
    values_from = ExtractedValue
  )

## 4.3 Combine into a single risk-of-bias table per study

rob_study <- study_labels %>%
  select(study_id, StudyLabel, publication_year) %>%
  left_join(meta_wide %>% select(study_id, N_total), by = "study_id") %>%
  left_join(rob_tools,              by = "study_id") %>%
  left_join(rob_domains_free_text,  by = "study_id") %>%
  left_join(rob_overall,            by = "study_id") %>%
  left_join(rob_domains_detail,     by = "study_id") %>%
  arrange(publication_year)

############################################################
# PART A – Table 6A: Risk-of-bias summary by study
############################################################

table6A <- rob_study %>%
  transmute(
    Study = StudyLabel,
    `N (analysed)` =
      if_else(is.na(N_total), "Not reported", as.character(N_total)),
    `Risk-of-bias tool` =
      coalesce(risk_of_bias_tool, "Not formally specified"),
    `Domains assessed` =
      coalesce(risk_of_bias_domains, "Standard diagnostic domains"),
    `Overall risk of bias` =
      coalesce(risk_of_bias_overall, "Not graded study-level"),
    `Sources of bias (overall)` =
      coalesce(common_sources_of_bias, "Not specifically summarised"),
    `Patient selection` =
      coalesce(common_bias_patient_selection, "No major concern highlighted"),
    `Index test / predictors` =
      coalesce(common_bias_index_test, "No major concern highlighted"),
    `Reference standard / outcome` =
      coalesce(common_bias_reference_standard, "No major concern highlighted"),
    `Flow & timing` =
      coalesce(common_bias_flow_timing, "No major concern highlighted")
  )

table6A_kable <- kbl(
  table6A,
  caption = "Table 6A. Risk-of-bias tools and domain-level concerns across included studies.",
  align   = c("l", "c", "l", "l", "l", "l", "l", "l", "l", "l"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 8
  )

table6A_kable
# save_kable(table6A_kable, "table6A_risk_of_bias_summary.html")

############################################################
# PART B – Figure 6A: Risk-of-bias tools used across the evidence base
############################################################

rob_tool_counts <- rob_study %>%
  filter(!is.na(risk_of_bias_tool)) %>%
  count(risk_of_bias_tool) %>%
  arrange(desc(n)) %>%
  mutate(
    risk_of_bias_tool = fct_reorder(risk_of_bias_tool, n)
  )

fig6A <- ggplot(rob_tool_counts,
                aes(x = n, y = risk_of_bias_tool)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n),
            hjust = -0.1,
            size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Figure 6A. Risk-of-bias tools used in included studies",
    x     = "Number of studies",
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0),
    axis.text.y = element_text(size = 9)
  )

fig6A
# ggsave("figure6A_rob_tools_barplot.png", fig6A,
#        width = 6, height = 3.5, dpi = 300)

############################################################
# PART C – Figure 6B: Heatmap of bias domains with concerns
############################################################

# Create flags: did that study have a specific bias note in each domain?
rob_domains_binary <- rob_study %>%
  mutate(
    has_patient_selection_issue   = !is.na(common_bias_patient_selection),
    has_index_test_issue          = !is.na(common_bias_index_test),
    has_reference_standard_issue  = !is.na(common_bias_reference_standard),
    has_flow_timing_issue         = !is.na(common_bias_flow_timing)
  ) %>%
  select(
    StudyLabel,
    has_patient_selection_issue,
    has_index_test_issue,
    has_reference_standard_issue,
    has_flow_timing_issue
  )

rob_domains_long <- rob_domains_binary %>%
  pivot_longer(
    cols      = starts_with("has_"),
    names_to  = "Domain",
    values_to = "Issue"
  ) %>%
  mutate(
    Domain = recode(
      Domain,
      "has_patient_selection_issue"  = "Patient selection",
      "has_index_test_issue"         = "Index test / predictors",
      "has_reference_standard_issue" = "Reference standard / outcome",
      "has_flow_timing_issue"        = "Flow & timing"
    ),
    Issue = if_else(Issue, "Issue described", "No specific issue")
  )

# Keep order of studies as in rob_study (publication year ascending)
rob_domains_long <- rob_domains_long %>%
  mutate(
    StudyLabel = factor(
      StudyLabel,
      levels = rob_study$StudyLabel
    ),
    Domain = factor(
      Domain,
      levels = c(
        "Patient selection",
        "Index test / predictors",
        "Reference standard / outcome",
        "Flow & timing"
      )
    )
  )

fig6B <- ggplot(rob_domains_long,
                aes(x = Domain, y = StudyLabel, fill = Issue)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Issue described" = "#d73027",
               "No specific issue" = "#f0f0f0"),
    name   = "Domain-level concern"
  ) +
  labs(
    title = "Figure 6B. Bias domains with explicitly described concerns",
    x     = "Risk-of-bias domain",
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

fig6B
# ggsave("figure6B_bias_domain_heatmap.png", fig6B,
#        width = 6.5, height = 4.5, dpi = 300)

############################################################
# PART D – Publication bias assessment (Cao 2025 meta-analysis)
############################################################

pub_bias <- dat_annot %>%
  filter(VariableName == "publication_bias_assessment") %>%
  group_by(study_id) %>%
  summarise(
    publication_bias_assessment = ExtractedValue[1],
    pub_bias_notes              = Notes[1],
    .groups = "drop"
  ) %>%
  left_join(
    study_labels %>% select(study_id, StudyLabel),
    by = "study_id"
  )

table6B <- pub_bias %>%
  transmute(
    Study = StudyLabel,
    `Publication bias assessment` =
      coalesce(publication_bias_assessment, "Not formally assessed"),
    `Details` = coalesce(pub_bias_notes, "")
  )

table6B_kable <- kbl(
  table6B,
  caption = "Table 6B. Publication bias assessment in biomarker or ultrasound meta-analyses.",
  align   = c("l", "l", "l"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 9
  )

table6B_kable
# save_kable(table6B_kable, "table6B_publication_bias.html")

############################################################
# NEW: Save all figures and tables into the output folder
############################################################

# Save tables as HTML
save_kable(
  table6A_kable,
  file = file.path(out_dir, "table6A_risk_of_bias_summary.html")
)

save_kable(
  table6B_kable,
  file = file.path(out_dir, "table6B_publication_bias.html")
)

# Save figures as PNG (same dimensions as your commented examples)
ggsave(
  filename = file.path(out_dir, "figure6A_rob_tools_barplot.png"),
  plot     = fig6A,
  width    = 6,
  height   = 3.5,
  dpi      = 300
)

ggsave(
  filename = file.path(out_dir, "figure6B_bias_domain_heatmap.png"),
  plot     = fig6B,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)
