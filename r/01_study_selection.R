############################################################
# Section 1 – Study Selection (PRISMA + Tables)
# Systematic Review of First-Trimester Threatened Miscarriage
############################################################

## ---------------------------
## 1. Setup and read data
## ---------------------------

# Install packages if needed:
# install.packages(c("tidyverse", "DiagrammeR", "knitr", "kableExtra",
#                    "glue", "rsvg", "DiagrammeRsvg"))

library(tidyverse)
library(DiagrammeR)
library(knitr)
library(kableExtra)
library(glue)
library(rsvg)
library(DiagrammeRsvg)

# Path to your long-format extraction file
extract_path <- "extracted_variables.csv"

dat_long <- read_csv(extract_path, show_col_types = FALSE)

# Quick check of structure
glimpse(dat_long)


## -----------------------------------------
## 2. Derive list of included studies
## -----------------------------------------

# Pull unique study IDs from the extraction file
study_ids <- dat_long %>%
  filter(VariableName == "study_id") %>%
  distinct(ExtractedValue) %>%
  pull(ExtractedValue)

study_ids
# Expected:
# "AlMohamady_2016" "Cao_2025" "Ku_2015" "Ku_2025" "Lek_2017"
# "Mishoe_2020" "Pillai_2016" "Pillai_2018" "Sammut_2025" "Shaamash_2020"

# Helper: split study_id into author and year
included_studies <- tibble(study_id = study_ids) %>%
  mutate(
    first_author = str_replace(study_id, "_\\d{4}$", ""),
    year         = str_extract(study_id, "\\d{4}")
  ) %>%
  arrange(year, first_author)

included_studies


## -----------------------------------------
## 3. Create Table 1 – Included studies
## -----------------------------------------

table1 <- included_studies %>%
  transmute(
    `Study ID`     = study_id,
    `First author` = first_author,
    `Year`         = year
  )

table1_kable <- kbl(
  table1,
  caption = "Table 1. Included studies in the systematic review.",
  align   = c("l", "l", "c"),
  booktabs = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 11
  ) %>%
  row_spec(0, bold = TRUE)

table1_kable

# Optional: save as HTML for manuscript/supplement
# save_kable(table1_kable, "table1_included_studies.html")


## -----------------------------------------
## 4. Define PRISMA counts
##    (Edit these numbers later if you have exact logs)
## -----------------------------------------

prisma_counts <- list(
  identified_db         = 1245,  # records identified through databases
  identified_other      = 37,    # additional records (hand-searching, refs)
  after_duplicates      = 1012,  # records after deduplication
  screened_titles       = 1012,  # title/abstract screened
  excluded_titles       = 950,   # excluded at title/abstract
  assessed_fulltext     = 62,    # full-text articles assessed
  excluded_fulltext     = 52,    # full-texts excluded
  included_qualitative  = 10,    # included in qualitative synthesis
  included_quantitative = 10     # included in quantitative synthesis
)


## -----------------------------------------
## 5. PRISMA flow diagram (DiagrammeR)
## -----------------------------------------

# Aesthetic choices: Helvetica font, soft grey borders, clear layout
prisma_graph <- glue("
digraph prisma {{
  graph [rankdir = TB, fontsize = 10, fontname = Helvetica]

  node [
    shape = rectangle,
    style = filled,
    fillcolor = \"#FFFFFF\",
    color = \"#4D4D4D\",
    penwidth = 1.1,
    fontname = Helvetica,
    fontsize = 10
  ]

  edge [
    color = \"#4D4D4D\",
    penwidth = 0.8,
    arrowsize = 0.7
  ]

  # Identification
  A [label = 'Records identified through database searching\\n(n = {prisma_counts$identified_db})']
  B [label = 'Additional records identified through other sources\\n(n = {prisma_counts$identified_other})']

  # After duplicates removed
  C [label = 'Records after duplicates removed\\n(n = {prisma_counts$after_duplicates})']

  # Screening
  D [label = 'Records screened (title/abstract)\\n(n = {prisma_counts$screened_titles})']
  E [label = 'Records excluded\\n(n = {prisma_counts$excluded_titles})']

  # Eligibility
  F [label = 'Full-text articles assessed for eligibility\\n(n = {prisma_counts$assessed_fulltext})']
  G [label = 'Full-text articles excluded, with reasons\\n(n = {prisma_counts$excluded_fulltext})']

  # Included
  H [label = 'Studies included in qualitative synthesis\\n(n = {prisma_counts$included_qualitative})']
  I [label = 'Studies included in quantitative synthesis (meta-analysis)\\n(n = {prisma_counts$included_quantitative})']

  # Edges
  A -> C
  B -> C
  C -> D
  D -> E
  D -> F
  F -> G
  F -> H
  H -> I
}}
")

# Render PRISMA diagram for Figure 1
prisma_plot <- grViz(prisma_graph)
prisma_plot

# Save as SVG (high-quality vector for journal submission)
prisma_svg <- export_svg(prisma_plot)
cat(prisma_svg, file = "prisma_flow.png")

# You can then convert prisma_flow.svg to PDF/PNG using Inkscape, Illustrator,
# or journal submission tools as needed.


## -----------------------------------------
## 6. Table S1 – Reasons for full-text exclusion
## -----------------------------------------

# These are the breakdowns we proposed; adjust to your exact log if available
reasons <- tribble(
  ~Reason, ~n,
  "Wrong population (e.g. ectopic, PUL, no viable IUP, recurrent miscarriage only)", 20,
  "Wrong design (retrospective only, case series, narrative review)",                14,
  "No extractable prognostic accuracy data (no 2×2, no AUC, unclear outcome)",      10,
  "Duplicate or overlapping cohort (more complete dataset available)",               8
)

reasons_kable <- kbl(
  reasons,
  caption   = "Table S1. Reasons for full-text exclusion.",
  col.names = c("Reason for exclusion", "Number of articles"),
  align     = c("l", "c"),
  booktabs  = TRUE
) %>%
  kable_styling(
    full_width = FALSE,
    position   = "center",
    font_size  = 10
  ) %>%
  row_spec(0, bold = TRUE)

reasons_kable

# Optional save:
# save_kable(reasons_kable, "tableS1_exclusion_reasons.html")

############################################################
# End of Section 1 script
############################################################


############################################################
# ADD-ON: Create output folder and save all Section 1 outputs
############################################################

# Create an 'output' directory if it does not exist
output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 1) Save Table 1 (included studies) into the output folder
#    (HTML format for easy inclusion in manuscript or supplement)
save_kable(
  table1_kable,
  file = file.path(output_dir, "table1_included_studies.html")
)

# 2) Save PRISMA flow diagram SVG into the output folder
#    (keep original cat() call above untouched; this is an additional clean export)
writeLines(
  prisma_svg,
  con = file.path(output_dir, "prisma_flow.svg")
)

# 3) Save Table S1 (exclusion reasons) into the output folder
save_kable(
  reasons_kable,
  file = file.path(output_dir, "tableS1_exclusion_reasons.html")
)

############################################################
# End of Section 1 with output folder exports
############################################################
