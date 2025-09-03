# Load Libraries
library(tidyverse) 
library(broom)
library(gtsummary)
library(survival)
library(survminer)
library(rstpm2)
library(flextable)

# --- 1. Data Loading ---
D <- read.csv("criteria 1.csv", header = TRUE, stringsAsFactors = FALSE)
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
  dplyr::select(
    Age, Race, `Tumor grade`, ER, PR, HER2, `HR/HER2 status`, Stage,
    `Mastectomy type`, Radiotherapy, Chemotherapy, `Breast cancer deaths` = BCSS
  )
# No rows are removed in this step, only columns are selected.
head(D_filtered)
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

# print(table_export) # View table in console. TABLE 1 for manuscript.
save_as_docx(table_export, path = "table_export_simplified.docx") # Uncomment to save as Word

# --- 6. Prepare Data for Survival Analysis ---
# Select relevant columns and rename for clarity
D_for_survival <- D_filtered %>%
  dplyr::select(
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

get_reverse_KM_stats <- function(time, event) {
  fit <- survfit(Surv(time, event) ~ 1)  # overall fit
  surv_summary <- summary(fit)
  
  # Median follow-up (time where survival = 0.5)
  median_followup <- summary(fit)$table["median"]
  median_LCL <- summary(fit)$table["0.95LCL"]
  median_UCL <- summary(fit)$table["0.95UCL"]
  
  # To get IQR from KM curve, approximate by times at survival ~0.75 and 0.25
  # Find closest times where survival crosses 0.75 and 0.25
  times <- surv_summary$time
  surv_prob <- surv_summary$surv
  
  q75_time <- min(times[surv_prob <= 0.75])
  q25_time <- min(times[surv_prob <= 0.25])
  
  return(tibble(
    median = median_followup,
    median_LCL = median_LCL,
    median_UCL = median_UCL,
    Q3 = q75_time,
    Q1 = q25_time,
    IQR = q75_time - q25_time
  ))
}
overall_stats <- get_reverse_KM_stats(
  time = D_for_survival$Time,
  event = D_for_survival$BCSS == 0
)
print(overall_stats)

# Check for NAs in the final survival dataset (should ideally be 0 for modeling variables)
cat("\nNA counts in D_for_survival (columns for modeling):\n")
sapply(D_for_survival %>% dplyr::select(Age, Grade, HRHER2, Stage, Time, BCSS), function(x) sum(is.na(x))) %>% print()

cat("\nStructure of D_for_survival:\n")
str(D_for_survival)

# End of script. D_for_survival is ready for survival modeling.
# table_export contains the descriptive statistics table.

# Check excluded cohort
# --- 2. Data Filtering ---
filter_conditions <- list(
  quote(!is.na(`Tumor grade`)),
  quote(!is.na(Stage)),
  quote(!is.na(ER)),
  quote(!is.na(PR)),
  quote(!is.na(HER2)),
  quote(!is.na(`Mastectomy type`)),
  quote(Time != 0),
  quote(!is.na(Time))
)

D_final_filtered <- D
for (condition_expr in filter_conditions) {
  D_final_filtered <- D_final_filtered %>% filter(!!condition_expr)
}

cohort_filtered_data <- D_final_filtered
cohort_filtered_out_data <- anti_join(D, cohort_filtered_data, by = names(D))

# --- 3. Combine Cohorts + Variables ---
combined_cohorts_for_gtsummary <- bind_rows(
  cohort_filtered_data %>% mutate(cohort_type = "Filtered Data"),
  cohort_filtered_out_data %>% mutate(cohort_type = "Filtered Out Data")
) %>%
  mutate(cohort_type = factor(cohort_type, levels = c("Filtered Data", "Filtered Out Data"))) %>%
  mutate(
    Time_FollowUp_Value = if_else(Time > 0, Time, NA_real_),
    
    Age_Group = cut(Age,
                    breaks = c(0, 40, 50, 60, 70, 80, Inf),
                    labels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+"),
                    right = FALSE,
                    include.lowest = TRUE),
    Age_Group = factor(Age_Group, levels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+")),
    
    Year_Diagnosis_Group = case_when(
      Year.of.diagnosis < 2016 ~ "Before 2016",
      Year.of.diagnosis >= 2016 ~ "2016 and After",
      TRUE ~ NA_character_
    ),
    Year_Diagnosis_Group = factor(Year_Diagnosis_Group, levels = c("Before 2016", "2016 and After"))
  )

# --- 4. gtsummary Table ---
comparison_table <- combined_cohorts_for_gtsummary %>%
  dplyr::select(
    cohort_type,
    Age_Group,
    Race,
    Year_Diagnosis_Group,
    `Tumor grade`,
    Stage,
    ER,
    PR,
    HER2,
    `Mastectomy type`,
    Time_FollowUp_Value,
    `Breast cancer deaths` = BCSS
  ) %>%
  tbl_summary(
    by = cohort_type,
    percent = "column",
    statistic = list(
      Time_FollowUp_Value ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany",
    label = list(
      Age_Group ~ "Age group (Years)",
      Race ~ "Race/Ethnicity",
      Year_Diagnosis_Group ~ "Year of diagnosis",
      `Tumor grade` ~ "Tumor grade",
      Stage ~ "Stage",
      ER ~ "ER",
      PR ~ "PR",
      HER2 ~ "HER2",
      `Mastectomy type` ~ "Mastectomy Type",
      Time_FollowUp_Value ~ "Time: Follow-up > 0 (Months)",
      `Breast cancer deaths` ~ "Breast cancer deaths"
    )
  ) %>%
  add_p() %>%
  add_n() %>%
  modify_header(label = "**Characteristic**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Cohort**") %>%
  modify_caption("**Table 1. Comparison of Filtered Data vs. Filtered Out Data**")

# --- 5. Custom: Change "Missing" → "Time = 0" for Time only ---
comparison_table <- comparison_table %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        label = if_else(
          variable == "Time_FollowUp_Value" & row_type == "missing",
          "Time = 0",
          label
        )
      )
  ) %>%
  bold_labels() %>%
  as_flex_table()

# --- Display Table ---
comparison_table
save_as_docx(comparison_table, path = "comparison_table.docx")


# --- Cox Proportional Hazards Model Fitting and Assumption Check ---
# Step 1: Fit Cox model
cox_main <- coxph(Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage, data = D_for_survival)
cox.zph(cox_main)
summary(cox_main)


# --- Flexible Parametric Survival Models (stpm2) ---

# Define common TVC degrees of freedom to avoid repetition
# This applies to all models where TVC is (Age = X, Grade = X, Stage = X, HRHER2 = X)
tvc_df_2 <- list(Age = 5, Grade = 2,  HRHER2 = 3, Stage = 3)

# --- Model 1: Baseline model with main effects and fixed baseline hazard flexibility ---
model_stpm2_I <- stpm2(
  Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage,
  data = D_for_survival,
  df = 4,
  tvc = tvc_df_2
)
summary(model_stpm2_I)

# --- Model 2: Adding interaction between Grade and HRHER2 ---
model_stpm2_II<- stpm2(
  Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage + Grade*HRHER2,
  data = D_for_survival,
  df = 4,
  tvc = tvc_df_2
)
summary(model_stpm2_II)

# --- Model 3: Adding more interaction terms, fixed baseline hazard flexibility ---
model_stpm2_III <- stpm2(
  Surv(Time, BCSS) ~ Age + Grade + HRHER2 + Stage + Grade*HRHER2 + Grade*Stage + HRHER2*Stage,
  data = D_for_survival,
  df = 4,
  tvc = tvc_df_2
)
summary(model_stpm2_III)

# --- Model Comparison ---
AIC(model_stpm2_I, model_stpm2_II, model_stpm2_III) 
BIC(model_stpm2_I, model_stpm2_II, model_stpm2_III)

# Saving the models
save(model_stpm2_I, model_stpm2_II, model_stpm2_III, file = "stpm2_modelsIII.RData")
# Loading the models
load("stpm2_modelsIII.RData")

final_model <- model_stpm2_III

time_point <- 60
predicted_hazards <- predict(final_model,
                             newdata = D_for_survival,
                             type = "hazard",
                             time = time_point)

# Check for negative hazard values
negative_hazards <- predicted_hazards[predicted_hazards < 0]

if (length(negative_hazards) == 0) {
  cat("✅ No negative hazard values found at time =", time_point, "\n")
} else {
  cat("⚠️ WARNING: Negative hazard values detected at time =", time_point, "!\n")
  
  neg_indices <- which(predicted_hazards < 0)
  
  cat("\nDetails of negative hazard predictions:\n")
  for (idx in neg_indices) {
    cat(sprintf("Individual (row in D_for_survival): %d, Predicted Hazard: %.6f\n",
                idx, predicted_hazards[idx]))
    
    # Optionally print covariates for that individual:
    # print(D_for_survival[idx, ])
  }
}


# Define function to extract survival probabilities from stpm2 models
library(broom)
library(dplyr)
library(pec)

# Define simplified function to extract survival probabilities from stpm2 models
predictSurvProb.stpm2 <- function(object, newdata, times, ...) {
  # Basic checks
  stopifnot(is.data.frame(newdata), nrow(newdata) > 0, length(times) > 0)
  
  # Prepare expanded newdata and corresponding times
  all_newdata_rows <- do.call(rbind, replicate(length(times), newdata, simplify = FALSE))
  all_times <- rep(times, each = nrow(newdata))
  all_newdata2 <- data.frame(time = all_times)
  
  # Predict survival probabilities directly
  predictions_vector <- predict(object,
                                newdata = all_newdata_rows,
                                type = "surv",
                                newdata2 = all_newdata2)
  
  # Reshape vector into matrix (subjects x times)
  p_matrix <- matrix(predictions_vector,
                     nrow = nrow(newdata),
                     ncol = length(times),
                     byrow = FALSE)
  
  # 5. Add robust dimension check (inspired by predictSurvProb.coxph)
  if (NROW(p_matrix) != NROW(newdata) || NCOL(p_matrix) != length(times)) {
    stop(paste("\nInternal Error: Prediction matrix has wrong dimensions after reshaping.\nRequested newdata x times: ",
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ",
               NROW(p_matrix), " x ", NCOL(p_matrix), "\n\n", sep = ""))
  }
  
  # Return the matrix of survival probabilities
  return(p_matrix)
}

# Compute prediction error curves (Brier Score) for model_stpm2_3
pec_res1 <- pec(
  object = list(fpm = final_model),  # Specify the survival model for evaluation
  formula = Surv(Time, BCSS) ~ 1,  # Survival formula defining the event and time variable
  data = D_for_survival,  # Dataset for performance evaluation
  times = seq(12, 96, by = 12),  # Define evaluation time points
  cens.model = "marginal",  # Handle censoring using Kaplan-Meier estimation
  exact = FALSE,
  splitMethod = "none"  # Compute apparent performance without cross-validation
)

# Extract and summarize prediction error results
crps(pec_res1)[1:2]  # Compute Continuous Ranked Probability Score (related to IBS)
summary(pec_res1, times = seq(12, 96, by = 12))  # Summarize Brier scores at specific intervals

# Load additional packages for AUC calculations
library(parallel)
library(risksetROC)
library(rstpm2)

# Define evaluation time points
time_points <- seq(12, 96, by = 12)

# Set up parallel processing for AUC estimation
n_chunks <- 4  # Define number of parallel computation chunks
chunk_size <- ceiling(nrow(D_for_survival) / n_chunks)
indices <- split(1:nrow(D_for_survival), ceiling(seq_along(1:nrow(D_for_survival)) / chunk_size))

# Extract survival times and event status
Stime <- D_for_survival$Time
status <- D_for_survival$BCSS

run_roc_on_chunk <- function(idx, predict_time) {
  survival_probabilities <- predict(
    final_model,
    newdata = D_for_survival[idx, ],
    type = "surv",
    time = predict_time
  )
  
  marker <- 1 - survival_probabilities
  
  # Compute time-dependent AUC using the Schoenfeld method
  risksetROC(
    Stime = Stime[idx],
    status = status[idx],
    marker = marker,
    method = "Schoenfeld",
    predict.time = predict_time,
    span = 1
  )
}

# Initialize storage for AUC results
auc_summary <- data.frame(
  time = time_points,
  AUC = NA_real_,
  CI_lower = NA_real_,
  CI_upper = NA_real_
)

# Compute AUC for each evaluation time point
for (t in time_points) {
  cat("Processing time =", t, "\n")
  
  # Set up parallel computation cluster
  cl <- makeCluster(n_chunks)
  
  # Export necessary data and functions to worker nodes
  clusterExport(cl, varlist = c("D_for_survival", "final_model", "Stime", "status", "run_roc_on_chunk"))
  clusterEvalQ(cl, {
    library(risksetROC)
    library(rstpm2)
  })
  
  # Perform parallel AUC estimation for survival prediction models
  results <- parLapply(cl, indices, run_roc_on_chunk, predict_time = t)
  stopCluster(cl)
  
  # Extract estimated AUC values
  auc_values <- sapply(results, function(res) res$AUC)
  mean_auc <- mean(auc_values)
  
  # Bootstrap confidence intervals for the AUC estimate
  set.seed(123)
  n_boot <- 1000
  boot_means <- replicate(n_boot, {
    sample_auc <- sample(auc_values, replace = TRUE)
    mean(sample_auc)
  })
  
  ci_lower <- quantile(boot_means, 0.025)
  ci_upper <- quantile(boot_means, 0.975)
  
  # Store computed AUC and confidence interval values
  auc_summary[auc_summary$time == t, "AUC"] <- mean_auc
  auc_summary[auc_summary$time == t, "CI_lower"] <- ci_lower
  auc_summary[auc_summary$time == t, "CI_upper"] <- ci_upper
}

# Display final AUC results across all evaluation time points
print(auc_summary)


# --- Extract model coefficients, calculate HRs and CIs, and combine results ---

# Extract coefficients for each model, calculate HRs and CIs, and add a 'model' column
tidy_m1 <- tidy(model_stpm2_I, conf.int = TRUE) %>%
  mutate(HR = exp(estimate), CI_low = exp(conf.low), CI_high = exp(conf.high), model = "Model 1")

tidy_m2 <- tidy(model_stpm2_II, conf.int = TRUE) %>%
  mutate(HR = exp(estimate), CI_low = exp(conf.low), CI_high = exp(conf.high), model = "Model 2")

tidy_m3 <- tidy(model_stpm2_III, conf.int = TRUE) %>%
  mutate(HR = exp(estimate), CI_low = exp(conf.low), CI_high = exp(conf.high), model = "Model 3")


# Combine processed results from all models into a single data frame.
tidy_all <- bind_rows(
  tidy_m1 %>% dplyr::select(term, HR, CI_low, CI_high, p.value, model),
  tidy_m2 %>% dplyr::select(term, HR, CI_low, CI_high, p.value, model),
  tidy_m3 %>% dplyr::select(term, HR, CI_low, CI_high, p.value, model)
)
df <- tidy_all

# Save the combined, unprocessed data (optional, based on your needs)
write.csv(df, "df.csv", row.names = FALSE)

# Filter out intercept and spline terms
df_filtered <- df %>%
  filter(term != "(Intercept)", !grepl("nsx", term, fixed = TRUE))

# Save the filtered data
write.csv(df_filtered, "df_filtered.csv", row.names = FALSE)

  
  ## Prepare Data for Plotting
  # Identify unique terms present in each model to categorize them later.
terms_m1 <- df_filtered %>% filter(model == "Model 1") %>% pull(term) %>% unique()
terms_m2 <- df_filtered %>% filter(model == "Model 2") %>% pull(term) %>% unique()
terms_m3 <- df_filtered %>% filter(model == "Model 3") %>% pull(term) %>% unique()


# Determine common and unique terms across models for grouping.
common_all    <- intersect(intersect(terms_m1, terms_m2), terms_m3)
shared_m2_m3  <- setdiff(intersect(terms_m2, terms_m3), terms_m1)
unique_m3     <- setdiff(terms_m3, union(terms_m1, terms_m2))

# Assign terms to broader categories like "Main Effects" or "Interaction Terms".
df_labeled <- df_filtered %>%
  mutate(group = case_when(
    term %in% common_all ~ "Main Effects",
    grepl(":", term) ~ "Interaction Terms",
    TRUE ~ "Other Terms" # Fallback, though ideally all terms should be categorized
  ))


# Create a complete grid of all terms and models to ensure all possible
# combinations are represented, even if a term is absent in a specific model.
all_terms <- df_labeled %>%
  dplyr::select(term, group) %>%
  distinct()

full_grid <- expand.grid(term = unique(all_terms$term),
                         model = c("Model 1", "Model 2", "Model 3"),
                         stringsAsFactors = FALSE) %>%
  left_join(all_terms, by = "term") %>%
  left_join(df_labeled, by = c("term", "model", "group"))


# Format terms for display on the Y-axis, grouping them by category.
full_grid <- full_grid %>%
  mutate(
    group = case_when(
      term %in% common_all ~ "Main Effects",
      grepl(":", term) ~ "Interaction Terms",
      TRUE ~ "Other Terms"
    ),
    group = factor(group, levels = c("Main Effects", "Interaction Terms"))
  ) %>%
  mutate(
    interaction_type_for_sort = case_when(
      grepl("^GradeG[0-9]:HRHER2", term) ~ "A_Grade_HRHER2_Interactions",
      grepl("^GradeG[0-9]:Stage", term) ~ "B_Grade_Stage_Interactions",
      grepl("^HRHER2.+:Stage", term) ~ "C_HRHER2_Stage_Interactions",
      TRUE ~ "D_Other_Interactions"
    ),
    grade_numeric_for_sort = case_when(
      grepl("GradeG1", term) ~ 1,
      grepl("GradeG2", term) ~ 2,
      grepl("GradeG3", term) ~ 3,
      TRUE ~ NA_real_
    ),
    hrher2_order_for_sort = case_when(
      grepl("HR\\+/HER2-", term) ~ 1,
      grepl("HR\\+/HER2\\+", term) ~ 2,
      grepl("HR-/HER2\\+", term) ~ 3,
      grepl("HR-/HER2-", term) ~ 4,
      TRUE ~ NA_real_
    ),
    stage_numeric_for_sort = case_when(
      grepl("StageI", term) ~ 1,
      grepl("StageII", term) ~ 2,
      grepl("StageIII", term) ~ 3,
      grepl("StageIV", term) ~ 4,
      TRUE ~ NA_real_
    ),
    main_sort_flag = case_when(
      term == "Age" ~ 0,
      term == "GradeG2" ~ 1,
      term == "GradeG3" ~ 2,
      TRUE ~ 99
    )
  ) %>%
  arrange(
    group,
    main_sort_flag,
    interaction_type_for_sort,
    hrher2_order_for_sort,
    grade_numeric_for_sort,
    stage_numeric_for_sort,
    term
  ) %>%
  mutate(
    term_clean = term %>%
      gsub("HRHER2", "", .) %>%
      gsub("GradeG", "G", .)
  ) %>%
  mutate(
    term_display = ifelse(
      is.na(group),
      NA_character_,
      paste(group, term_clean, sep = ": ")
    )
  ) %>%
  mutate(term_display = fct_rev(fct_inorder(term_display))) %>%
  dplyr::select(-c(interaction_type_for_sort, grade_numeric_for_sort,
            hrher2_order_for_sort, stage_numeric_for_sort,
            main_sort_flag, term_clean))



# Assign significance labels and colors based on p-value for plotting.
full_grid <- full_grid %>%
  mutate(sig_label = ifelse(!is.na(p.value) & p.value < 0.05, "Significant (p < 0.05)", "Not Significant"),
         sig_color = ifelse(sig_label == "Significant (p < 0.05)", "red", "black"))

F1<-ggplot(full_grid, aes(x = HR, y = term_display)) +
  
  # Add a vertical reference line at HR = 1 (no effect).
  geom_vline(xintercept = 1, linetype = "solid", color = "grey30", linewidth = 0.6) +
  
  # Add points for Hazard Ratios and horizontal error bars for confidence intervals.
  geom_point(aes(color = sig_label), na.rm = TRUE, size = 2) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high, color = sig_label),
    height = 0.2,
    na.rm = TRUE,
    linewidth = 0.8
  ) +
  
  # Define colors for significant and non-significant terms.
  scale_color_manual(
    name = "Significance",
    values = c("Significant (p < 0.05)" = "red", "Not Significant" = "black")
  ) +
  
  # Use a logarithmic scale for the X-axis for better visualization of HRs.
  scale_x_log10() +
  
  # Create separate panels (facets) for each model to compare HRs across models.
  facet_grid(. ~ model, labeller = label_value) +
  
  # Apply a clean theme and add plot labels.
  theme_bw() +
  labs(
    x = "Hazard Ratio (log scale)",
    y = "Term",
    title = "Hazard Ratio by Different Models"
  ) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )
F1

ggsave("Figure1.jpg", plot = F1, width = 8, height = 6, units = "in", dpi = 300)

# --- Global Plotting Variable Definitions ---
# Extract all unique levels for key categorical variables from the D_for_survival dataset.
hrher2_levels <- levels(D_for_survival$HRHER2)
stage_levels <- levels(D_for_survival$Stage)
grade_levels <- levels(D_for_survival$Grade)

# Define a sequence of time points for predictions to create smooth curves.
time_seq <- seq(min(D_for_survival$Time), 100, length.out = 100)

# Initialize lists and data frames to store plots and extracted data across all iterations.
hrher2_plot_list <- list() # Stores combined plots for each HR/HER2 status group.
all_surv_data <- data.frame() # Accumulates all survival prediction data for export.
all_hr_data <- data.frame() # Accumulates all hazard ratio data for export.

# --- Loop Through HR/HER2 Status (Outer Loop) ---
# This loop iterates through each HR/HER2 status (e.g., HR+/HER2-), creating a set of plots for each.
for (current_hrher2 in hrher2_levels) {
  
  # --- Survival Curve Generation (Left Panels) ---
  # This nested loop generates survival curves for each cancer stage,
  # keeping HR/HER2 status constant for the current outer loop iteration.
  surv_plots <- list() # Stores survival plots for the current HR/HER2 status.
  
  for (current_stage in stage_levels) {
    surv_all <- data.frame() # Accumulates survival data for different grades within the current stage.
    
    # Innermost loop: Generate survival curves for each tumor grade.
    for (current_grade in grade_levels) {
      # Create a 'newdata' data frame for prediction, holding other variables constant
      # (e.g., Age fixed at median, HR/HER2 and Stage at current loop values).
      newdata_single <- data.frame(
        Grade = factor(current_grade, levels = grade_levels),
        Age = median(D_for_survival$Age, na.rm = TRUE), # Use D_for_survival
        HRHER2 = factor(current_hrher2, levels = hrher2_levels),
        Stage = factor(current_stage, levels = stage_levels),
        Time = time_seq # Predict across the defined time sequence.
      )
      
      # Predict survival probability and its confidence intervals from model_stpm2_3.
      pred <- predict(final_model, newdata = newdata_single, type = "surv", se.fit = TRUE) # Model updated
      newdata_single$survival <- pred$Estimate
      newdata_single$lower <- pred$lower
      newdata_single$upper <- pred$upper
      newdata_single$Grade <- current_grade # Re-attach Grade for plotting aesthetics.
      
      surv_all <- rbind(surv_all, newdata_single) # Add predictions for this grade to current stage's data.
    }
    
    # Store metadata for the accumulated survival data and append to the global survival data frame.
    newdata_single$HRHER2 <- current_hrher2
    newdata_single$Stage <- current_stage
    all_surv_data <- rbind(all_surv_data, surv_all)
    
    # Create the title for the current survival panel.
    panel_title <- paste(current_hrher2, ":", "Stage", current_stage)
    
    # Generate the ggplot for survival curves.
    p_surv <- ggplot(surv_all, aes(x = Time, y = survival * 100, color = Grade)) +
      geom_line(linewidth = 0.5) + # Plot the survival curve.
      geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100, fill = Grade), alpha = 0.2, color = NA) + # Add confidence intervals.
      labs(
        title = panel_title, # Panel-specific title.
        x = NULL, # X-axis label will be set at the very bottom of the combined plot.
        y = "Survival (%)", # Y-axis label.
        color = "Grade", # Legend title.
        fill = "Grade" # Fill legend title.
      ) +
      ylim(0, 100) + # Set Y-axis limit for survival percentage.
      scale_color_manual(values = c("G1" = "blue", "G2" = "orange", "G3" = "red")) + # Custom colors for grades.
      scale_fill_manual(values = c("G1" = "blue", "G2" = "orange", "G3" = "red")) + # Custom fill colors for grades.
      theme_bw() + # Black and white theme.
      theme(legend.position = "none", 
            plot.title = element_text(size = 7, hjust = 0),
            axis.title.x = element_text(size = 6), # Set font size for x-axis title
            axis.title.y = element_text(size = 6), # Set font size for y-axis title
            axis.text.x = element_text(size = 5),  # Set font size for x-axis text (numbers)
            axis.text.y = element_text(size = 5)   # Set font size for y-axis text (numbers)
      ) # Hide individual legends, style title.
    
    surv_plots[[current_stage]] <- p_surv # Store the generated survival plot for this stage.
  }
  
  
  # --- Hazard Ratio Curve Generation (Right Panels) ---
  # This nested loop generates time-dependent Hazard Ratio curves for different grade comparisons,
  # keeping HR/HER2 status constant for the current outer loop iteration.
  hr_plots <- list() # Stores HR plots for the current HR/HER2 status.
  
  for (current_stage in stage_levels) {
    # Define a reference profile (Grade G1) for HR calculation.
    newdata_ref <- data.frame(
      Grade = factor("G1", levels = grade_levels),
      Age = median(D_for_survival$Age, na.rm = TRUE), # Use D_for_survival
      HRHER2 = factor(current_hrher2, levels = hrher2_levels),
      Stage = factor(current_stage, levels = stage_levels),
      Time = time_seq
    )
    
    # Predict Hazard Ratios (HR) for G2 vs G1, G3 vs G1, and G3 vs G2.
    # The 'exposed' argument allows dynamic comparison.
    hr_g2 <- predict(final_model, newdata = newdata_ref, type = "hr", # Model updated
                     exposed = function(data) transform(data, Grade = factor("G2", levels = grade_levels)),
                     se.fit = TRUE)
    hr_g3 <- predict(final_model, newdata = newdata_ref, type = "hr", # Model updated
                     exposed = function(data) transform(data, Grade = factor("G3", levels = grade_levels)),
                     se.fit = TRUE)
    
    # Need to correctly define the reference for G3 vs G2
    newdata_ref_g2 <- transform(newdata_ref, Grade = factor("G2", levels = grade_levels)) # Reference for G2
    hr_g3g2 <- predict(final_model, newdata = newdata_ref_g2, type = "hr", # Model updated
                       exposed = function(data) transform(data, Grade = factor("G3", levels = grade_levels)),
                       se.fit = TRUE)
    
    
    # Combine HR data for all comparisons into a single data frame for plotting.
    hr_data_sub <- rbind(
      data.frame(Time = time_seq, HR = hr_g2$Estimate, Upper = hr_g2$upper, Lower = hr_g2$lower, Comparison = "G2 vs G1"),
      data.frame(Time = time_seq, HR = hr_g3$Estimate, Upper = hr_g3$upper, Lower = hr_g3$lower, Comparison = "G3 vs G1"),
      data.frame(Time = time_seq, HR = hr_g3g2$Estimate, Upper = hr_g3g2$upper, Lower = hr_g3g2$lower, Comparison = "G3 vs G2")
    )
    
    # Generate the ggplot for Hazard Ratio curves.
    p_hr <- ggplot(hr_data_sub, aes(x = Time, y = HR, color = Comparison, fill = Comparison)) +
      geom_line(linewidth = 0.5) + # Plot the HR curve.
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + # Add confidence intervals.
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + # Reference line at HR=1 (no effect).
      labs(
        x = "Time (Months)", # X-axis label.
        y = "Hazard Ratio", # Y-axis label.
        color = "Comparison", # Legend title.
        fill = "Comparison" # Fill legend title.
      ) +
      scale_color_manual(values = c("G2 vs G1" = "blue", "G3 vs G1" = "red", "G3 vs G2" = "orange")) + # Custom colors for comparisons.
      scale_fill_manual(values = c("G2 vs G1" = "blue", "G3 vs G1" = "red", "G3 vs G2" = "orange")) + # Custom fill colors.
      theme_bw() + # Black and white theme.
      theme(
        legend.position = "none", # Hide individual legends.
        axis.title.x = element_text(size = 6), # Set font size for x-axis title
        axis.title.y = element_text(size = 6), # Set font size for y-axis title
        axis.text.x = element_text(size = 5),  # Set font size for x-axis text (numbers)
        axis.text.y = element_text(size = 5)   # Set font size for y-axis text (numbers)
      )
    hr_plots[[current_stage]] <- p_hr # Store the generated HR plot for this stage.
    
    # Store metadata for the accumulated HR data and append to the global HR data frame.
    hr_data_sub$HRHER2 <- current_hrher2
    hr_data_sub$Stage <- current_stage
    all_hr_data <- rbind(all_hr_data, hr_data_sub)
  }
  
  # --- Combine Survival and HR Plots for Each Stage (using patchwork) ---
  # This loop combines the survival plot (left) and HR plot (right) for each stage.
  stage_combined_plots <- list()
  for (current_stage in stage_levels) {
    stage_combined_plots[[current_stage]] <- patchwork::wrap_plots(surv_plots[[current_stage]], hr_plots[[current_stage]], ncol = 2) # Added patchwork:: prefix
  }
  # Combine all stage-specific combined plots into a single row for the current HR/HER2 group.
  hrher2_plot_list[[current_hrher2]] <- patchwork::wrap_plots(stage_combined_plots, nrow = 1) # Added patchwork:: prefix
}

# --- Post-Processing for Legends (Important for Readability) ---
# This section ensures that legends are only displayed in the first plot of each row,
# preventing clutter in the final combined layout.
for (i in seq_along(hrher2_plot_list)) {
  # Add legend to the first survival plot (first plot of the first combined stage panel) in each HRHER2 row.
  hrher2_plot_list[[i]][[1]][[1]] <- hrher2_plot_list[[i]][[1]][[1]] +
    theme(legend.position = "bottom")
  
  # Add legend to the first hazard ratio plot (second plot of the first combined stage panel) in each HRHER2 row.
  hrher2_plot_list[[i]][[1]][[2]] <- hrher2_plot_list[[i]][[1]][[2]] +
    theme(legend.position = "bottom")
}

# Stack all HR/HER2 groups vertically into the final grand plot.
final_plot <- patchwork::wrap_plots(hrher2_plot_list, ncol = 1) + # Added patchwork:: prefix
  patchwork::plot_layout(guides = "collect") & # Added patchwork:: prefix
  theme(legend.position = "bottom") # Position the collected legend at the bottom of the final plot

ggsave("Figure2.jpg", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)

write.csv(all_surv_data, "all_surv_data.csv", row.names = FALSE)
write.csv(all_hr_data, "all_hr_data.csv", row.names = FALSE)
survival_data_at_60_months <- all_surv_data[all_surv_data$Time == 60, ]

# Compute survival difference between G3 and G1
survival_data_at_60_months %>%
  filter(Grade %in% c("G1", "G3")) %>%
  group_by(HRHER2, Stage) %>%
  summarise(
    drop = survival[Grade == "G1"] - survival[Grade == "G3"]
  ) %>%
  arrange(desc(drop))

###Maximun HR
# Initialize an empty list to store the maximum HR rows
max_hr_rows_list <- list()
counter <- 1 # Counter for the list index

# Get the unique values from the HRHER2 column (e.g., HR+/HER2-, HR-/HER2+, etc.).
unique_subtypes <- unique(all_hr_data$HRHER2)

# Outer loop: Iterate through each unique HR/HER2 subtype.
for (subtype in unique_subtypes) {
  # Filter the complete hazard ratio data for the current HR/HER2 subtype.
  subset_subtype_data <- all_hr_data[all_hr_data$HRHER2 == subtype, ]
  
  # Check if there's any data for the current subtype before proceeding.
  if (nrow(subset_subtype_data) > 0) {
    # Get the unique stages present within this specific HR/HER2 subtype.
    unique_stages <- unique(subset_subtype_data$Stage)
    
    # Inner loop: Iterate through each unique cancer stage within the current HR/HER2 subtype.
    for (stage in unique_stages) {
      # Filter the subtype data further for the current stage.
      subset_stage_data <- subset_subtype_data[subset_subtype_data$Stage == stage, ]
      
      # Check if there's any data for the current stage.
      if (nrow(subset_stage_data) > 0) {
        # Find the row index where the HR value is at its maximum within this specific stage subgroup.
        max_hr_index <- which.max(subset_stage_data$HR)
        
        # Extract the entire row corresponding to the maximum HR.
        max_hr_row <- subset_stage_data[max_hr_index, ]
        
        # Store the relevant columns of this row in our list
        max_hr_rows_list[[counter]] <- data.frame(
          HRHER2 = max_hr_row$HRHER2,
          Stage = max_hr_row$Stage,
          Comparison = max_hr_row$Comparison,
          Time_at_Max_HR = max_hr_row$Time,
          Max_HR_Estimate = max_hr_row$HR,
          Max_HR_Lower_CI = max_hr_row$Lower,
          Max_HR_Upper_CI = max_hr_row$Upper
        )
        counter <- counter + 1 # Increment the counter for the next entry
      }
    }
  }
}

# Combine all stored rows into a single dataframe
maximum_hr_table <- do.call(rbind, max_hr_rows_list)

# Print the resulting table
print(maximum_hr_table)


library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)
library(scales) # Required for scales::number_format
hrher2_levels <- levels(D_for_survival$HRHER2)
stage_levels <- levels(D_for_survival$Stage)
grade_levels <- levels(D_for_survival$Grade)


# RMST plots
# Tau values to evaluate
tau_values <- c(24, 60, 96)

# Initialize empty data frames to store all results
all_rmst_values <- data.frame()
all_rmst_differences <- data.frame()

for (tau in tau_values) {
  cat("Processing tau =", tau, "\n")
  
  # Reuse your existing code here, just adding `tau` to each result set
  rmst_values <- data.frame()
  rmst_results <- list()
  
  for (current_hrher2 in hrher2_levels) {
    for (current_stage in stage_levels) {
      for (current_grade in grade_levels) {
        newdata <- data.frame(
          Grade = factor(current_grade, levels = grade_levels),
          Age = median(D_for_survival$Age, na.rm = TRUE),
          HRHER2 = factor(current_hrher2, levels = hrher2_levels),
          Stage = factor(current_stage, levels = stage_levels),
          Time = tau
        )
        
        pred_rmst <- predict(final_model, newdata = newdata, type = "rmst", se.fit = TRUE, t = tau, type.relsurv = NULL)
        
        rmst_values <- rbind(
          rmst_values,
          data.frame(
            HRHER2 = current_hrher2,
            Stage = current_stage,
            Grade = current_grade,
            RMST = pred_rmst$Estimate,
            Lower = pred_rmst$lower,
            Upper = pred_rmst$upper,
            tau = tau
          )
        )
      }
    }
  }
  
  # Save all RMST values
  all_rmst_values <- rbind(all_rmst_values, rmst_values)
  
  # RMST difference logic
  for (hr in hrher2_levels) {
    for (stg in stage_levels) {
      age_med <- median(D_for_survival$Age, na.rm = TRUE)
      subgroup_results <- list()
      
      comparisons <- list(
        list(ref = "G1", exp = "G3"),
        list(ref = "G1", exp = "G2"),
        list(ref = "G2", exp = "G3")
      )
      
      for (comp in comparisons) {
        ref_grade <- comp$ref
        exp_grade <- comp$exp
        
        newdata_ref <- data.frame(
          Grade = factor(ref_grade, levels = grade_levels),
          Age = age_med,
          HRHER2 = factor(hr, levels = hrher2_levels),
          Stage = factor(stg, levels = stage_levels),
          Time = tau
        )
        
        rmst_diff <- tryCatch({
          predict(
            final_model,
            newdata = newdata_ref,
            type = "rmstdiff",
            exposed = function(data) transform(data, Grade = factor(exp_grade, levels = grade_levels)),
            t = tau,
            se.fit = TRUE, type.relsurv = NULL
          )
        }, error = function(e) return(NULL))
        
        if (!is.null(rmst_diff)) {
          subgroup_results[[length(subgroup_results) + 1]] <- data.frame(
            HRHER2 = hr,
            Stage = stg,
            Comparison = paste0(exp_grade, " vs ", ref_grade),
            Estimate = -rmst_diff$Estimate,
            Lower = -rmst_diff$lower,
            Upper = -rmst_diff$upper,
            tau = tau
          )
        }
      }
      
      if (length(subgroup_results) > 0) {
        subgroup_df <- do.call(rbind, subgroup_results)
        subgroup_df$p.value <- pnorm(abs(subgroup_df$Estimate) / ((subgroup_df$Upper - subgroup_df$Lower) / (2 * qnorm(0.975))))
        subgroup_df$adjusted.p.value.BH <- p.adjust(subgroup_df$p.value, method = "BH")
        subgroup_df$adjusted.p.value.Bonferroni <- p.adjust(subgroup_df$p.value, method = "bonferroni")
        subgroup_df$significant.Bonferroni <- subgroup_df$adjusted.p.value.Bonferroni < 0.05
        rmst_results[[length(rmst_results) + 1]] <- subgroup_df
      }
    }
  }
  
  rmst_summary_adjusted_subgroup <- do.call(rbind, rmst_results)
  
  # Save all RMST differences
  all_rmst_differences <- rbind(all_rmst_differences, rmst_summary_adjusted_subgroup)
}

# View or save
head(all_rmst_values)
head(all_rmst_differences)
write.csv(all_rmst_differences, "all_rmst_differences.csv", row.names = FALSE)
write.csv(all_rmst_values, "all_rmst_values.csv", row.names = FALSE)

# Order by descending Estimate (from most negative to least)
df_ranked <- all_rmst_differences %>%
  filter(Comparison == "G3 vs G1", tau == 60) %>%
  group_by(Stage) %>%
  arrange(Estimate, .by_group = TRUE) %>%  # ascending: most negative to least
  mutate(Rank = row_number()) %>%
  ungroup()

# View the result
print(df_ranked)

df_ranked96 <- all_rmst_differences %>%
  filter(Comparison == "G3 vs G1", tau == 96) %>%
  group_by(Stage) %>%
  arrange(Estimate, .by_group = TRUE) %>%  # ascending: most negative to least
  mutate(Rank = row_number()) %>%
  ungroup()

print(df_ranked96)


# Define plot_rmst_stage function
plot_rmst_stage <- function(stage_value, hrher2_value, data_df, diff_limits) {
  is_stage1 <- stage_value == "I"
  
  rmst_stage <- data_df %>%
    filter(Stage == stage_value, HRHER2 == hrher2_value) %>%
    droplevels()
  
  if (nrow(rmst_stage) == 0) {
    warning(paste("No data for Stage:", stage_value, "and HRHER2:", hrher2_value))
    return(NULL) # Return NULL if no data
  }
  
  # Define the desired order for Comparison levels
  comparison_order <- c("G3 vs G1", "G2 vs G1", "G3 vs G2")
  
  rmst_stage <- rmst_stage %>%
    mutate(
      Comparison = factor(Comparison, levels = comparison_order), # Set factor levels here
      Grade1 = str_extract(Comparison, "^G[123]"),
      Grade2 = str_extract(Comparison, "(?<=vs )G[123]")
    )
  
  segments_data <- rmst_stage %>%
    rowwise() %>%
    mutate(
      RMST1 = get(paste0("RMST_", Grade1)),
      RMST2 = get(paste0("RMST_", Grade2))
    ) %>%
    ungroup()
  
  points_data <- segments_data %>%
    dplyr::select(Comparison, tau, Grade1, Grade2, RMST1, RMST2) %>%
    pivot_longer(cols = c(RMST1, RMST2), names_to = "RMST_order", values_to = "RMST") %>%
    mutate(
      Grade = ifelse(RMST_order == "RMST1", Grade1, Grade2)
    ) %>%
    filter(!is.na(RMST)) %>%
    dplyr::select(Comparison, tau, Grade, RMST)
  
  tau_levels_ordered <- sort(unique(rmst_stage$tau))
  
  tau_bg_df_for_rect <- data.frame(
    tau = factor(tau_levels_ordered), # Convert tau to a factor here
    bg_color_value = rep(c("white", "gray95"), length.out = length(tau_levels_ordered))
  )
  
  fill_colors_map <- setNames(tau_bg_df_for_rect$bg_color_value, tau_bg_df_for_rect$tau)
  
  p_rmst <- ggplot(data = rmst_stage) +
    geom_rect(
      data = tau_bg_df_for_rect,
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = tau),
      inherit.aes = FALSE,
      alpha = 0.4 # Increased alpha for better visibility of gray
    ) +
    scale_fill_manual(values = fill_colors_map, guide = "none") +
    
    facet_wrap(~ tau, ncol = 1, strip.position = "right") +
    theme_minimal() +
    theme(strip.text.y = element_blank()) +
    
    geom_segment(data = segments_data,
                 aes(x = RMST1, xend = RMST2, y = Comparison, yend = Comparison),
                 color = "gray50", linewidth = 1) +
    
    geom_point(data = points_data,
               aes(x = RMST, y = Comparison, color = Grade),
               size = 1) +
    
    scale_color_manual(values = c("G1" = "green", "G2" = "orange", "G3" = "red")) +
    
    labs(
      x = "RMST (months)",
      y = if (is_stage1) "Comparison" else NULL,
      color = "Grade",
      title = paste(hrher2_value, ": Stage", stage_value)
    ) +
    scale_x_continuous(
      limits = c(15, 96),
      breaks = seq(15, 95, by = 10),
      labels = scales::number_format(accuracy = 1), # Ensures no decimals
      expand = c(0, 0)
    ) +
    theme(axis.title.x = element_text(size = 6),
          axis.text.x = element_text(size = 6)) +
    theme_minimal(base_size = 6) +
    theme(
      axis.text.y = if (is_stage1) element_text(size = 6) else element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = if (is_stage1) element_text(size = 6) else element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray", linetype = "dashed", size = 0.4), # Added x-grid line
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      legend.position = "none", # Set legend position to none for individual plots
      plot.margin = margin(5, 5, 5, 5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           fill = "none")
  
  return(p_rmst)
}

# Define plot_rmst_diff_bar function
plot_rmst_diff_bar <- function(stage_value, hrher2_value, data_df, diff_limits) {
  # Define the desired order for Comparison levels
  comparison_order <- c("G3 vs G1", "G2 vs G1", "G3 vs G2")
  
  rmst_stage <- data_df %>%
    filter(Stage == stage_value, HRHER2 == hrher2_value) %>%
    droplevels() %>%
    mutate(
      Comparison = factor(Comparison, levels = comparison_order), # Set factor levels here
      tau_label = paste0("τ:", tau, "m"),
      Significant = ifelse(p_bonf < 0.05, "*", "")
    )
  
  if (nrow(rmst_stage) == 0) {
    warning(paste("No diff data for Stage:", stage_value, "and HRHER2:", hrher2_value))
    return(NULL) # Return NULL if no data
  }
  
  tau_levels <- sort(unique(rmst_stage$tau_label))
  tau_bg_df <- data.frame(
    tau_label = tau_levels,
    xmin = diff_limits[[stage_value]][1],
    xmax = diff_limits[[stage_value]][2],
    ymin = 0.5,
    ymax = length(levels(factor(rmst_stage$Comparison))) + 0.5,
    bg_color = rep(c("white", "gray95"), length.out = length(tau_levels))
  )
  
  p_diff <- ggplot() +
    geom_rect(data = tau_bg_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = bg_color),
              inherit.aes = FALSE, color = NA) +
    scale_fill_identity() +
    
    geom_col(data = rmst_stage, aes(x = RMST_Diff, y = Comparison), fill = "#e31a1c", width = 0.5) +
    
    geom_text(data = rmst_stage, aes(label = Significant, x = RMST_Diff, y = Comparison),
              hjust = -0.2, color = "black", size = 3) +
    
    facet_wrap(~ tau_label, ncol = 1, strip.position = "right") +
    
    labs(x = "RMST Difference (months)", y = "Comparison") +
    scale_x_continuous(
      limits = diff_limits[[stage_value]],
      expand = c(0, 0)
    ) +
    
    theme_minimal(base_size = 11) +
    theme(
      strip.text.y = if (stage_value == "IV")
        element_text(size = 5)
      else
        element_blank(),
      
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      panel.grid.major.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      plot.margin = margin(5, 10, 5, 5),
      legend.position = "none" # Set legend position to none for individual diff plots
    )
  
  return(p_diff)
}
# Assuming 'all_rmst_values' and 'all_rmst_differences' are already loaded
# and contain data for various HRHER2 levels.

# Step 1: Combine all RMST and Difference data into a single 'final_df'
# This 'final_df' will contain data for ALL HRHER2 levels
hr_rmst_all <- all_rmst_values %>%
  dplyr::select(HRHER2, Stage, Grade, RMST, RMST_Lower = Lower, RMST_Upper = Upper, tau = tau)

hr_diff_all <- all_rmst_differences %>%
  mutate(Grade = sub(" .*", "", Comparison)) %>%
  dplyr::select(HRHER2, Stage, Grade, Comparison,
         RMST_Diff = Estimate,
         Diff_Lower = Lower,
         Diff_Upper = Upper,
         p_value = p.value,
         p_bonf = adjusted.p.value.Bonferroni,
         tau = tau)

all_grades <- c("G1", "G2", "G3")
final_rows_global <- list()

for (i in 1:nrow(hr_diff_all)) {
  row <- hr_diff_all[i, ]
  
  grades_in_comp <- unlist(strsplit(as.character(row$Comparison), " vs "))
  
  rmst_vals <- setNames(rep(NA, length(all_grades)), paste0("RMST_", all_grades))
  rmst_lower <- setNames(rep(NA, length(all_grades)), paste0("RMST_", all_grades, "_Lower"))
  rmst_upper <- setNames(rep(NA, length(all_grades)), paste0("RMST_", all_grades, "_Upper"))
  
  for (g in grades_in_comp) {
    rmst_row <- subset(hr_rmst_all, HRHER2 == row$HRHER2 & Stage == row$Stage & Grade == g & tau == row$tau)
    if (nrow(rmst_row) == 1) {
      rmst_vals[paste0("RMST_", g)] <- rmst_row$RMST
      rmst_lower[paste0("RMST_", g, "_Lower")] <- rmst_row$RMST_Lower
      rmst_upper[paste0("RMST_", g, "_Upper")] <- rmst_row$RMST_Upper
    }
  }
  
  combined_row <- data.frame(
    HRHER2 = row$HRHER2,
    Stage = row$Stage,
    tau = row$tau,
    Comparison = row$Comparison,
    RMST_Diff = row$RMST_Diff,
    Diff_Lower = row$Diff_Lower,
    Diff_Upper = row$Diff_Upper, 
    p_value = row$p_value,
    p_bonf = row$p_bonf
  )
  
  combined_row <- cbind(combined_row, 
                        as.data.frame(as.list(rmst_vals)), 
                        as.data.frame(as.list(rmst_lower)), 
                        as.data.frame(as.list(rmst_upper)))
  
  final_rows_global[[length(final_rows_global) + 1]] <- combined_row
}

# This 'final_df' now contains data for all HRHER2 levels
final_df_global <- do.call(rbind, final_rows_global)

write.csv(final_df_global, "final_df_global.csv", row.names = FALSE)

# Calculate stage_x_diff_limits once for all stages across ALL HRHER2 levels
stage_x_diff_limits_global <- list()
for (s in unique(final_df_global$Stage)) {
  subset_diffs <- final_df_global %>% 
    dplyr::filter(Stage == s)
  
  if (nrow(subset_diffs) > 0) {
    min_val <- min(subset_diffs$RMST_Diff, na.rm = TRUE)
    max_val <- max(subset_diffs$RMST_Diff, na.rm = TRUE)
    padding <- 0.01 * (max_val - min_val)
    lower_bound <- floor(min_val - padding)
    upper_bound <- ceiling(max_val + padding)
    stage_x_diff_limits_global[[s]] <- c(lower_bound, upper_bound)
  } else {
    stage_x_diff_limits_global[[s]] <- c(-10, 10) 
  }
}


# Get all unique HRHER2 levels to loop through
all_hrher2_levels <- unique(final_df_global$HRHER2)

# List to store combined plots for each HRHER2 level
all_combined_plots <- list()

# Loop through each HRHER2 level
for (current_hrher2_value in all_hrher2_levels) {
  
  # Filter the global final_df for the current HRHER2 level for plotting
  current_hrher2_data_for_plotting <- final_df_global %>%
    filter(HRHER2 == current_hrher2_value)
  
  # Generate plots for all stages for the current HRHER2 level
  # Pass the filtered data and the GLOBAL diff_limits
  plot_I <- (plot_rmst_stage("I", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global) +
               plot_rmst_diff_bar("I", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global)) +
    plot_layout(widths = c(2, 1)) # Set 2:1 ratio for stage and diff plots
  plot_II <- (plot_rmst_stage("II", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global) +
                plot_rmst_diff_bar("II", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global)) +
    plot_layout(widths = c(2, 1)) # Set 2:1 ratio for stage and diff plots
  plot_III <- (plot_rmst_stage("III", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global) +
                 plot_rmst_diff_bar("III", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global)) +
    plot_layout(widths = c(2, 1)) # Set 2:1 ratio for stage and diff plots
  plot_IV <- (plot_rmst_stage("IV", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global) +
                plot_rmst_diff_bar("IV", current_hrher2_value, current_hrher2_data_for_plotting, stage_x_diff_limits_global)) +
    plot_layout(widths = c(2, 1)) # Set 2:1 ratio for stage and diff plots
  
  # Combine the stage-specific plots horizontally for the current HRHER2 level
  # Removed guides = "collect" and theme(legend.position = "bottom") from here
  combined_plot_for_hrher2 <- (plot_I | plot_II | plot_III | plot_IV) +
    plot_layout(ncol = 4, widths = c(1, 1, 1, 1)) & # Removed guides = "collect"
    theme(plot.title = element_text(hjust = 0.5, size = 8)) # Removed legend.position
  
  # Store the combined plot in the list, named by the HRHER2 level
  all_combined_plots[[current_hrher2_value]] <- combined_plot_for_hrher2
}

# Combine all generated plots vertically
final_stacked_plot <- wrap_plots(all_combined_plots, ncol = 1) + # Stack all plots vertically
  plot_layout(guides = "collect") & # Collect all legends from all plots
  theme(legend.position = "bottom") # Place the single collected legend at the bottom

# Display the final stacked plot
final_stacked_plot
ggsave("Figure3.jpg", plot = final_stacked_plot, width = 8, height = 6, units = "in", dpi = 300)

