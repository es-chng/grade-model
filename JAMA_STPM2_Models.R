# Load Libraries
library(tidyverse) # Includes dplyr, stringr, ggplot2, forcats
library(broom)
library(gtsummary)
library(survival)
library(survminer)
library(rstpm2)
library(flextable) # For table export

# --- 1. Data Loading ---
D <- read.csv("GRADE/criteria 1.csv", header = TRUE, stringsAsFactors = FALSE)
n_initial_load <- nrow(D)
cat("Initial data loaded: ", n_initial_load, " rows\n")

# --- 2. Feature Engineering and Preprocessing ---
radio_values_yes <- c(
  "Beam radiation",
  "Radioactive implants (includes brachytherapy) (1988+)",
  "Combination of beam with implants or isotopes",
  "Radiation, NOS method or source not specified",
  "Radioisotopes (1988+)"
)
D <- D %>%
  mutate(
    # Age
    Age = as.numeric(str_sub(Age.recode.with.single.ages.and.85., 1, 2)),
    
    # Race
    race_numeric = case_when(
      Race.recode..W..B..AI..API. == "White" ~ 1,
      Race.recode..W..B..AI..API. == "Black" ~ 2,
      Race.recode..W..B..AI..API. == "Unknown" ~ NA_real_,
      TRUE ~ 3 # Others
    ),
    Race = factor(race_numeric, levels = 1:3, labels = c("White", "Black", "Others")),
    
    # Tumor Grade
    # Clean original grade columns
    Grade.Recode..thru.2017. = na_if(Grade.Recode..thru.2017., "Unknown"),
    Grade.Pathological..2018.. = na_if(Grade.Pathological..2018.., "Blank(s)"),
    # Combine grade information
    grade_numeric = case_when(
      Grade.Recode..thru.2017. == "Well differentiated; Grade I" & is.na(Grade.Pathological..2018..) ~ 1,
      Grade.Recode..thru.2017. == "Moderately differentiated; Grade II" & is.na(Grade.Pathological..2018..) ~ 2,
      Grade.Recode..thru.2017. == "Poorly differentiated; Grade III" & is.na(Grade.Pathological..2018..) ~ 3,
      is.na(Grade.Recode..thru.2017.) & Grade.Pathological..2018.. == "1" ~ 1, # Assuming Pathological grade is character
      is.na(Grade.Recode..thru.2017.) & Grade.Pathological..2018.. == "2" ~ 2,
      is.na(Grade.Recode..thru.2017.) & Grade.Pathological..2018.. == "3" ~ 3,
      TRUE ~ NA_real_
    ),
    `Tumor grade` = factor(grade_numeric, levels = 1:3, labels = c("G1", "G2", "G3")),
    
    # ER, PR, HER2 Status
    ER_numeric = case_when(
      ER.Status.Recode.Breast.Cancer..1990.. == "Positive" ~ 1,
      ER.Status.Recode.Breast.Cancer..1990.. == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    ER = factor(ER_numeric, levels = 0:1, labels = c("Negative", "Positive")),
    
    PR_numeric = case_when(
      PR.Status.Recode.Breast.Cancer..1990.. == "Positive" ~ 1,
      PR.Status.Recode.Breast.Cancer..1990.. == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    PR = factor(PR_numeric, levels = 0:1, labels = c("Negative", "Positive")),
    
    HER2_numeric = case_when(
      Derived.HER2.Recode..2010.. == "Positive" ~ 1,
      Derived.HER2.Recode..2010.. == "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    HER2 = factor(HER2_numeric, levels = 0:1, labels = c("Negative", "Positive")),
    
    # HR/HER2 Status Derivation
    HRHER2_combined = case_when( # Renamed to avoid conflict with D_for_survival$HRHER2
      (ER == "Positive" | PR == "Positive") & HER2 == "Negative" ~ "HR+/HER2-",
      (ER == "Positive" | PR == "Positive") & HER2 == "Positive" ~ "HR+/HER2+",
      (ER == "Negative" & PR == "Negative") & HER2 == "Positive" ~ "HR-/HER2+",
      (ER == "Negative" & PR == "Negative") & HER2 == "Negative" ~ "HR-/HER2-",
      TRUE ~ NA_character_
    ),
    `HR/HER2 status` = factor(HRHER2_combined, levels = c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "HR-/HER2-")),
    
    # Cancer Stage (AJCC, SEER, EOD)
    ajcc_cleaned = na_if(Derived.AJCC.Stage.Group..7th.ed..2010.2015., "Blank(s)"),
    ajcc_cleaned = na_if(ajcc_cleaned, "UNK Stage"),
    ajcc_cleaned = na_if(ajcc_cleaned, "0"),
    ajcc_numeric = recode(ajcc_cleaned,
                          "IA" = "1", "IB" = "1",
                          "IIA" = "2", "IIB" = "2",
                          "IIIA" = "3", "IIIB" = "3", "IIIC" = "3", "IIINOS" = "3",
                          "IV" = "4",
                          .default = NA_character_),
    
    seer_cleaned = str_sub(Derived.SEER.Cmb.Stg.Grp..2016.2017., 1, 1),
    seer_numeric = na_if(seer_cleaned, "B"),
    seer_numeric = na_if(seer_numeric, "9"),
    seer_numeric = na_if(seer_numeric, "N"),
    seer_numeric = na_if(seer_numeric, "0"),
    
    eod_cleaned = str_sub(Derived.EOD.2018.Stage.Group..2018.., 1, 1),
    eod_numeric = na_if(eod_cleaned, "B"),
    eod_numeric = na_if(eod_numeric, "9"),
    eod_numeric = na_if(eod_numeric, "0"),
    
    stage_unified = coalesce(ajcc_numeric, seer_numeric, eod_numeric),
    Stage = factor(stage_unified, levels = c("1", "2", "3", "4"), labels = c("I", "II", "III", "IV")),
    
    # Mastectomy Type
    `Mastectomy type` = case_when(
      RX.Summ..Surg.Prim.Site..1998.. %in% c(30:76) ~ "Total mastectomy",
      RX.Summ..Surg.Prim.Site..1998.. %in% c(20:24) ~ "Partial mastectomy",
      TRUE ~ NA_character_
    ),
    `Mastectomy type` = factor(`Mastectomy type`, levels = c("Partial mastectomy", "Total mastectomy")),
    
    # Radiotherapy and Chemotherapy
    radiotherapy_numeric = ifelse(Radiation.recode %in% radio_values_yes, 2, 1), # 2:Yes, 1:No/Unknown
    Radiotherapy = factor(radiotherapy_numeric, levels = 1:2, labels = c("No/unknown", "Yes")),
    
    chemotherapy_numeric = ifelse(Chemotherapy.recode..yes..no.unk. == "Yes", 2, 1), # 2:Yes, 1:No/Unknown
    Chemotherapy = factor(chemotherapy_numeric, levels = 1:2, labels = c("No/unknown", "Yes")),
    
    # Survival
    BCSS = ifelse(COD.to.site.recode.ICD.O.3.2023.Revision.Expanded..1999.. == "Breast", 1, 0), # 1:Death due to Breast Cancer, 0:Otherwise
    Time = Survival.months
  )

# Number of rows after feature engineering (should be the same as initial load)
n_after_feature_engineering <- nrow(D)
cat("Number of rows after feature engineering: ", n_after_feature_engineering, " rows (no rows removed in this step)\n")


# --- 3. Data Filtering (Detailed Breakdown - Simplified) ---
cat("\n3. Detailed Data Filtering for NAs and Time != 0:\n")
n_before_step3_filter <- nrow(D)
cat("   - Starting rows before any filter in Step 3: ", n_before_step3_filter, " rows\n")

filter_conditions <- list(
  "`!is.na(Tumor grade)`" = quote(!is.na(`Tumor grade`)),
  "`!is.na(Stage)`" = quote(!is.na(Stage)),
  "`!is.na(ER)`" = quote(!is.na(ER)),
  "`!is.na(PR)`" = quote(!is.na(PR)),
  "`!is.na(HER2)`" = quote(!is.na(HER2)),
  "`!is.na(Mastectomy type)`" = quote(!is.na(`Mastectomy type`)),
  "`Time != 0`" = quote(Time != 0),
  "`!is.na(Time)`" = quote(!is.na(Time))
)

D_filtered <- D

for (i in seq_along(filter_conditions)) {
  filter_name <- names(filter_conditions)[i]
  condition_expr <- filter_conditions[[i]]
  
  rows_before_this_filter <- nrow(D_filtered)
  D_filtered <- D_filtered %>% filter(!!condition_expr)
  rows_after_this_filter <- nrow(D_filtered)
  removed_by_this_filter <- rows_before_this_filter - rows_after_this_filter
  
  cat("     - After filter ", filter_name, ": ", rows_after_this_filter, " rows (Removed by this filter: ", removed_by_this_filter, ")\n")
}

cat("   - Total rows remaining after all Step 3 filters: ", nrow(D_filtered), " rows\n")
cat("   - Total rows removed in Step 3 (cumulative): ", n_before_step3_filter - nrow(D_filtered), " rows\n")


# --- 4. Select Columns for Descriptive Table ---
DT_analysis <- D_filtered %>%
  select(
    Age, Race, `Tumor grade`, ER, PR, HER2, `HR/HER2 status`, Stage,
    `Mastectomy type`, Radiotherapy, Chemotherapy
  )
# No rows are removed in this step, only columns are selected.

# --- 5. Table Generation and Export ---
# Ensure `HR/HER2 status` is a factor for the 'by' argument in tbl_summary
DT_analysis$`HR/HER2 status` <- factor(DT_analysis$`HR/HER2 status`)

table_export <- DT_analysis %>%
  tbl_summary(
    by = `HR/HER2 status`,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "oneway.test",
      all_categorical() ~ "chisq.test.no.correct" # This will apply to all *other* categorical vars
    ),
    # This is the key: exclude ER, PR, and HER2 by their column names
    include = -c(ER, PR, HER2),
    test.args = list(all_continuous() ~ list(var.equal = TRUE)),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**HR/HER2 Status**") %>%
  bold_labels() %>%
  as_flex_table()

# print(table_export) # View table in console
save_as_docx(table_export, path = "table_export_simplified.docx") # Uncomment to save as Word

# --- 6. Prepare Data for Survival Analysis ---
# Select relevant columns and rename for clarity
D_for_survival <- D_filtered %>%
  select(
    Age,
    Grade = `Tumor grade`,
    HRHER2 = `HR/HER2 status`, # This uses the factor created from HRHER2_combined
    Stage,
    Time, # Survival time in months
    BCSS  # Event status (Breast Cancer Specific Survival: 1=event, 0=censored)
  )

# Ensure factor levels are correctly ordered/defined for modeling
D_for_survival <- D_for_survival %>%
  mutate(
    HRHER2 = factor(HRHER2, levels = c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "HR-/HER2-")),
    Grade = factor(Grade, levels = c("G1", "G2", "G3")),
    Stage = factor(Stage, levels = c("I", "II", "III", "IV"))
  )

# Check for NAs in the final survival dataset (should ideally be 0 for modeling variables)
cat("\nNA counts in D_for_survival (columns for modeling):\n")
sapply(D_for_survival %>% select(Age, Grade, HRHER2, Stage, Time, BCSS), function(x) sum(is.na(x))) %>% print()

cat("\nStructure of D_for_survival:\n")
str(D_for_survival)

# End of script. D_for_survival is ready for survival modeling.
# table_export contains the descriptive statistics table.


# --- Cox Proportional Hazards Model Fitting and Assumption Check ---
# Step 1: Fit Cox model
cox_main <- coxph(Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage, data = D_for_survival)
cox.zph(cox_main)

# --- Flexible Parametric Survival Models (stpm2) ---
# Load the rstpm2 library for fitting flexible parametric survival models
library(rstpm2)

# Define common TVC degrees of freedom to avoid repetition
# This applies to all models where TVC is (Age = X, Grade = X, Stage = X, HRHER2 = X)
tvc_df_2 <- list(Age = 2, Grade = 2, Stage = 2, HRHER2 = 2)


# --- Model 1: Baseline model with main effects and fixed baseline hazard flexibility ---
model_stpm2_1 <- stpm2(
  Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage,
  data = D_for_survival,
  df = 4,
  tvc = tvc_df_2
)
summary(model_stpm2_1)

# --- Model 2: Adding interaction between Grade and HRHER2 ---
model_stpm2_2<- stpm2(
  Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage + Grade*HRHER2,
  data = D_for_survival,
  df = 4,
  tvc = tvc_df_2
)
summary(model_stpm2_2)

# --- Model 3: Adding more interaction terms, fixed baseline hazard flexibility ---
model_stpm2_3 <- stpm2(
  Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage + Grade*HRHER2 + Grade*Stage + HRHER2*Stage,
  data = D_for_survival,
  df = 4,
  tvc = tvc_df_2
)
summary(model_stpm2_3)

# --- Model Comparison ---
AIC(model_stpm2_1, model_stpm2_2, model_stpm2_3) 
BIC(model_stpm2_1, model_stpm2_2, model_stpm2_3)

# Saving the models
save(model_stpm2_1, model_stpm2_2, model_stpm2_3, file = "stpm2_models.RData")
# Loading the models
load("stpm2_models.RData")