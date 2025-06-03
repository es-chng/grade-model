library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)
library(scales) # Required for scales::number_format
hrher2_levels <- levels(D_for_survival$HRHER2)
stage_levels <- levels(D_for_survival$Stage)
grade_levels <- levels(D_for_survival$Grade)

###### Figure 4

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
        
        pred_rmst <- predict(model_stpm2_3, newdata = newdata, type = "rmst", se.fit = TRUE, t = tau, type.relsurv = NULL)
        
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
            model_stpm2_3,
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
    select(Comparison, tau, Grade1, Grade2, RMST1, RMST2) %>%
    pivot_longer(cols = c(RMST1, RMST2), names_to = "RMST_order", values_to = "RMST") %>%
    mutate(
      Grade = ifelse(RMST_order == "RMST1", Grade1, Grade2)
    ) %>%
    filter(!is.na(RMST)) %>%
    select(Comparison, tau, Grade, RMST)
  
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
  select(HRHER2, Stage, Grade, RMST, RMST_Lower = Lower, RMST_Upper = Upper, tau = tau)

hr_diff_all <- all_rmst_differences %>%
  mutate(Grade = sub(" .*", "", Comparison)) %>%
  select(HRHER2, Stage, Grade, Comparison,
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
ggsave("Figure4_JAMA.png", plot = final_stacked_plot, width = 8, height = 6, units = "in", dpi = 300)
