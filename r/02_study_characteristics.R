############################################################
# EXTRA: Create output folder and save Section 2 tables/figures
############################################################

# Create a dedicated output folder for Section 2
output_dir <- "section2_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## ---- Save tables ----

# Table 2A – primary cohort design & population
save_kable(
  table2A_kable,
  file = file.path(output_dir, "table2A_primary_design_setting.html")
)

# Table 2B – primary cohort N & outcomes
save_kable(
  table2B_kable,
  file = file.path(output_dir, "table2B_primary_sample_outcomes.html")
)

# Table 2C – evidence-synthesis/meta-analysis characteristics
save_kable(
  table2C_kable,
  file = file.path(output_dir, "table2C_evidence_synthesis.html")
)

## ---- Save figures ----

# Figure 2A – sample size per primary study
ggsave(
  filename = file.path(output_dir, "figure2A_primary_sample_size.png"),
  plot     = fig2A,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# Figure 2B – miscarriage rate per primary study
ggsave(
  filename = file.path(output_dir, "figure2B_primary_miscarriage_rates.png"),
  plot     = fig2B,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# Figure 2C – predictor domain heatmap
ggsave(
  filename = file.path(output_dir, "figure2C_domain_heatmap.png"),
  plot     = fig2C,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)
