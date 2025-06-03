# Loading the models
load("stpm2_models.RData")
library(broom)
library(dplyr)

##### AUC and IBS calculation
library(riskRegression) # Load the riskRegression package for survival model validation.

# Enables riskRegression to get risk probabilities from stpm2 models.
predictRisk.stpm2 <- function(object, newdata, times, ...) {
  sapply(times, function(t) {
    1 - predict(object, newdata = newdata, type = "surv", time = t)
  })
}

library(pec) # Load the pec package for prediction error curves and C-index.

# Enables pec to get survival probabilities from stpm2 models.
predictSurvProb.stpm2 <- function(object, newdata, times, ...) {
  all_newdata_rows <- do.call(rbind, replicate(length(times), newdata, simplify = FALSE))
  all_times <- rep(times, each = nrow(newdata))
  all_newdata2 <- data.frame(time = all_times)
  predictions <- predict(object, newdata = all_newdata_rows, type = "surv", newdata2 = all_newdata2)
  matrix(predictions, nrow = nrow(newdata), ncol = length(times), byrow = TRUE)
}

# Define time points for evaluation (in months).
times <- c(12, 24, 36, 48, 60, 72, 84, 96)

# Compute time-dependent AUC for model_stpm2_3.
score_result <- Score(
  object = list(stpm2_model = model_stpm2_3), # Using model_stpm2_3 as the specified model
  formula = Surv(Time, BCSS) ~ 1, # Specifies the survival outcome.
  data = D_for_survival, # The dataset for evaluation.
  times = times, # Evaluation time points.
  metrics = "AUC", # Requesting the AUC metric.
  summary = "risk",
  contrasts = FALSE # Disable pairwise model comparisons.
)

print(score_result) # Display the calculated AUC values.

# Compute prediction error curves (Brier Score) for model_stpm2_3.
pec_res1 <- pec(
  object = list(fpm = model_stpm2_3), # Using model_stpm2_3 as the specified model
  formula = Surv(Time, BCSS) ~ 1, # Specifies the survival outcome.
  data = D_for_survival, # The dataset for evaluation.
  times = seq(12, 96, length.out = 12), # Evaluation time points for the curve.
  cens.model = "marginal", # Uses Kaplan-Meier to handle censoring.
  exact = FALSE,
  splitMethod = "none" # Computes apparent performance (no cross-validation).
)

# Extract and summarize the prediction error results.
crps(pec_res1)[1:2] # Show the Continuous Ranked Probability Score (related to IBS).
summary(pec_res1, times = seq(12, 96, by = 12)) # Summarize Brier scores at specific intervals.



# --- Extract model coefficients, calculate HRs and CIs, and combine results ---

# Extract coefficients for each model, calculate HRs and CIs, and add a 'model' column
tidy_m1 <- tidy(model_stpm2_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate), CI_low = exp(conf.low), CI_high = exp(conf.high), model = "Model 1")

tidy_m2 <- tidy(model_stpm2_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate), CI_low = exp(conf.low), CI_high = exp(conf.high), model = "Model 2")

tidy_m3 <- tidy(model_stpm2_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate), CI_low = exp(conf.low), CI_high = exp(conf.high), model = "Model 3")


# Combine processed results from all models into a single data frame.
tidy_all <- bind_rows(
  tidy_m1 %>% select(term, HR, CI_low, CI_high, p.value, model),
  tidy_m2 %>% select(term, HR, CI_low, CI_high, p.value, model),
  tidy_m3 %>% select(term, HR, CI_low, CI_high, p.value, model)
)
df <- tidy_all

# Save the combined, unprocessed data (optional, based on your needs)
write.csv(df, "df.csv", row.names = FALSE)

# Filter out intercept and spline terms
df_filtered <- df %>%
  filter(term != "(Intercept)", !grepl("nsx", term, fixed = TRUE))

# Save the filtered data
write.csv(df_filtered, "df_filtered.csv", row.names = FALSE)


---
  
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
  select(term, group) %>%
  distinct()

full_grid <- expand.grid(term = unique(all_terms$term),
                         model = c("Model 1", "Model 2", "Model 3"),
                         stringsAsFactors = FALSE) %>%
  left_join(all_terms, by = "term") %>%
  left_join(df_labeled, by = c("term", "model", "group"))


# Format terms for display on the Y-axis, grouping them by category.
full_grid <- full_grid %>%
  mutate(
    group = factor(group, levels = c("Main Effects", "Interaction Terms", "Other Terms")),
    term_display = paste(group, term, sep = ": ")
  ) %>%
  arrange(group, term) %>%
  mutate(term_display = fct_rev(fct_inorder(term_display)))


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

ggsave("Figure1_JAMA.png", plot = F1, width = 8, height = 6, units = "in", dpi = 300)

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
      pred <- predict(model_stpm2_3, newdata = newdata_single, type = "surv", se.fit = TRUE) # Model updated
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
    hr_g2 <- predict(model_stpm2_3, newdata = newdata_ref, type = "hr", # Model updated
                     exposed = function(data) transform(data, Grade = factor("G2", levels = grade_levels)),
                     se.fit = TRUE)
    hr_g3 <- predict(model_stpm2_3, newdata = newdata_ref, type = "hr", # Model updated
                     exposed = function(data) transform(data, Grade = factor("G3", levels = grade_levels)),
                     se.fit = TRUE)
    
    # Need to correctly define the reference for G3 vs G2
    newdata_ref_g2 <- transform(newdata_ref, Grade = factor("G2", levels = grade_levels)) # Reference for G2
    hr_g3g2 <- predict(model_stpm2_3, newdata = newdata_ref_g2, type = "hr", # Model updated
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

ggsave("Figure2_JAMA.png", plot = final_plot, width = 8, height = 6, units = "in", dpi = 300)

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


# --- 1. Define Evaluation Parameters and Reference Grid ---
# Define the specific time point (in months) at which to calculate Hazard Ratios.
time_point <- 60  # 60 months (5 years)

# Get the median age from your dataset to use as a fixed covariate value for HR calculations.
median_age <- median(D_for_survival$Age, na.rm = TRUE)

# Create a comprehensive grid covering all combinations of HRHER2, Grade, and Stage.
base_grid <- expand.grid(
  HRHER2 = factor(c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "HR-/HER2-"),
                  levels = levels(D_for_survival$HRHER2)),
  Grade = factor(c("G1", "G2", "G3"), levels = levels(D_for_survival$Grade)),
  Stage = factor(c("I", "II", "III", "IV"), levels = levels(D_for_survival$Stage)),
  Age = median_age
)

# Create a unique identifier for each stratum.
base_grid$stratum_id <- with(base_grid, paste(HRHER2, Grade, Stage, sep = "_"))
# Identify the reference stratum: HR+/HER2-, Grade G1, Stage I.
base_grid$reference <- base_grid$stratum_id == "HR+/HER2-_G1_I"

# --- 2. Calculate Hazard Ratios at 60 Months ---
# Define the reference profile at 60 months.
newdata_ref <- base_grid[base_grid$reference, , drop = FALSE]
newdata_ref$Time <- time_point

# Calculate HRs for each subgroup.
hr_list <- lapply(1:nrow(base_grid), function(i) {
  row_i <- base_grid[i, ]
  
  if (row_i$reference) {
    return(data.frame(HR = 1, lower = 1, upper = 1))
  }
  
  # Define the 'exposed' function for predict.stpm2.
  exposed_fn <- function(data) {
    data$HRHER2 <- factor(row_i$HRHER2, levels = levels(D_for_survival$HRHER2))
    data$Grade  <- factor(row_i$Grade, levels = levels(D_for_survival$Grade))
    data$Stage  <- factor(row_i$Stage, levels = levels(D_for_survival$Stage))
    data$Age    <- row_i$Age
    data$Time   <- time_point
    return(data)
  }
  
  # Predict the hazard ratio.
  pred <- predict(model_stpm2_3, newdata = newdata_ref, type = "hr",
                  exposed = exposed_fn, se.fit = TRUE)
  
  return(data.frame(HR = pred$Estimate, lower = pred$lower, upper = pred$upper))
})

# Combine the HR results.
hr_df <- do.call(rbind, hr_list)
final_hr_table <- cbind(base_grid, Time = time_point, hr_df)

# Order for clarity.
final_hr_table <- final_hr_table[order(final_hr_table$HRHER2,
                                       final_hr_table$Stage,
                                       final_hr_table$Grade), ]

# Write to CSV.
write.csv(final_hr_table, "final_hr_table_60months.csv", row.names = FALSE)

# --- 1. Prepare data with HR and CI labels ---
heatmap_data_60m <- final_hr_table %>%
  mutate(
    HRHER2_display = factor(HRHER2, levels = c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "HR-/HER2-")),
    Grade_display = factor(Grade, levels = c("G1", "G2", "G3")),
    Stage_display = factor(Stage, levels = c("I", "II", "III", "IV")),
    HR_label = sprintf("%.2f", HR),
    CI_label = paste0("(", sprintf("%.2f", lower), "–", sprintf("%.2f", upper), ")")
  )

# --- 2. Plot heatmap of 60-month hazard ratios ---
F3<-ggplot(heatmap_data_60m, aes(x = Stage_display, y = HRHER2_display, fill = HR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  
  # Add HR value
  geom_text(aes(label = HR_label), color = "black", size = 2.2, vjust = -0.6) +
  
  # Add 95% CI below HR
  geom_text(aes(label = CI_label), color = "black", size = 1.8, vjust = 1.4) +
  
  facet_wrap(~ Grade_display, ncol = 3) +
  
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 1,
    limits = c(0, max(heatmap_data_60m$HR, na.rm = TRUE) * 1.05),
    name = "Hazard Ratio"
  ) +
  
  labs(
    title = "60-Month Hazard Ratios by Subtype, Stage, and Grade",
    subtitle = "Reference: HR+/HER2−, Stage I, Grade 1",
    x = "Stage",
    y = "HR/HER2 Subtype"
  ) +
  
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey95", color = NA),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
ggsave("Figure3_JAMA.png", plot = F3, width = 8, height = 6, units = "in", dpi = 300)

