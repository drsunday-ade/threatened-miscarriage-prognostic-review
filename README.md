# threatened-miscarriage-prognostic-review
Reproducible R code, anonymised study-level data, and outputs for a systematic review and meta-analysis of biochemical, ultrasound, and multivariable predictors of miscarriage after first-trimester threatened miscarriage.
# Threatened Miscarriage Prognostic Review

This repository contains the reproducible code, data extraction template, and analysis outputs for the systematic review and meta-analysis:

> **Do Progestational, CA-125, Ultrasound and Multivariable Models Predict Miscarriage after First-Trimester Threatened Miscarriage?**  
> Sunday Adetunji et al., 2025.

All code is written in **R**, organised by manuscript section, and generates the figures and tables used in the paper.

---

## Repository structure

```text
data/
  raw/
    Extraction_Template.csv       # empty extraction template
    extracted_variables.csv       # anonymised extracted dataset
  processed/                      # derived analysis datasets (if used)

R/
  01_study_selection.R            # PRISMA counts, study flow
  02_study_characteristics.R      # tables of study design, setting, outcomes
  03_section3_biomarkers.R        # biochemical predictors (progesterone, CA-125)
  04_section4_ultrasound.R        # ultrasound predictors
  05_section5_multivariable_models.R  # multivariable & ML models
  06_section6_risk_of_bias.R      # risk-of-bias and publication-bias analyses
  data.Rproj                      # RStudio project file

outputs/
  figures/
    prisma_flow.svg
    figure2A_primary_sample_size.png
    figure2B_primary_miscarriage_rates.png
    figure2C_domain_heatmap.png
    figure3A_progesterone_sens_spec.png
    figure3B_ca125_sens_spec.png
    figure4C_ultrasound_feature_heatmap.png
    figure5A_multivariable_model_auc.png
    figure5B_logistic_vs_rf.png
    figure5C_predictor_domain_heatmap.png
    figure6A_rob_tools_barplot.png
    figure6B_bias_domain_heatmap.png
  tables/
    table1_included_studies.html
    tableS1_exclusion_reasons.html
    table2A_primary_design_setting.html
    table2B_primary_sample_outcomes.html
    table2C_evidence_synthesis.html
    table3A_progesterone_performance.html
    table3B_ca125_performance.html
    table4A_ultrasound_features.html
    table4B_fhr_performance.html
    table4C_crl_gs_performance.html
    table5A_multivariable_models.html
    table6A_risk_of_bias_summary.html
    table6B_publication_bias.html
  section_reports/
    03_section3_biomarkers.html
    Section-5---Multivariable-and-Machine-Learning-Models.html
    Section-6.html

manuscript/
  final_manuscript.docx           # optional

