# Longitudinal Trajectory Analysis Using Latent Class Mixed Models
# -----------------------------------------------------------------
# Authors: Le Wang, Yifei Lin, Jin Huang
# Description: This script performs trajectory analysis for a single biomarker 
#              (e.g., red blood cell count) using LCMM (latent class mixed models).

# Load Required Packages
library(data.table)
library(openxlsx)
library(dplyr)
library(lcmm)
library(tidyr)
library(ggplot2)

#-----------------------------------------------
# Load and Prepare Data
#-----------------------------------------------

data <- fread("data/data_male/data_female.csv", encoding = "UTF-8")
data$dah <- as.numeric(data$dah)

# Select biomarker and time variables
df <- data[, .(dah, value = `Red Blood Cell Count`, time)]  # Replace with your biomarker if needed
df <- df[!value %in% c("", "/", "-----")]
df$value <- as.numeric(df$value)
df <- df[!is.na(value)]

# Retain subjects with at least 2 observations
subject_counts <- df[, .N, by = dah]
subjects_to_keep <- subject_counts[N >= 2]$dah
rundata <- df[dah %in% subjects_to_keep]

# Log transformation
rundata$value <- log(rundata$value)

#-----------------------------------------------
# Fit Baseline (1-class) Models
#-----------------------------------------------

m1     <- hlme(value ~ time,                          random = ~ 1 + time, subject = "dah", data = rundata)
m1.1   <- hlme(value ~ time + I(time^2),              random = ~ 1 + time, subject = "dah", data = rundata)
m1.1.1 <- hlme(value ~ time + I(time^2) + I(time^3),  random = ~ 1 + time, subject = "dah", data = rundata)

#-----------------------------------------------
# Fit Multi-class Models (2 to 5 classes)
#-----------------------------------------------

model_list <- list(
  # Linear
  m2 = hlme(value ~ time, random = ~ 1 + time, subject = "dah", ng = 2, mixture = ~ time, B = m1, data = rundata),
  m3 = hlme(value ~ time, random = ~ 1 + time, subject = "dah", ng = 3, mixture = ~ time, B = m1, data = rundata),
  m4 = hlme(value ~ time, random = ~ 1 + time, subject = "dah", ng = 4, mixture = ~ time, B = m1, data = rundata),
  m5 = hlme(value ~ time, random = ~ 1 + time, subject = "dah", ng = 5, mixture = ~ time, B = m1, data = rundata),
  
  # Quadratic
  m2.1 = hlme(value ~ time + I(time^2), random = ~ 1 + time, subject = "dah", ng = 2, mixture = ~ time + I(time^2), B = m1.1, data = rundata),
  m3.1 = hlme(value ~ time + I(time^2), random = ~ 1 + time, subject = "dah", ng = 3, mixture = ~ time + I(time^2), B = m1.1, data = rundata),
  m4.1 = hlme(value ~ time + I(time^2), random = ~ 1 + time, subject = "dah", ng = 4, mixture = ~ time + I(time^2), B = m1.1, data = rundata),
  m5.1 = hlme(value ~ time + I(time^2), random = ~ 1 + time, subject = "dah", ng = 5, mixture = ~ time + I(time^2), B = m1.1, data = rundata),
  
  # Cubic
  m2.1.1 = hlme(value ~ time + I(time^2) + I(time^3), random = ~ 1 + time, subject = "dah", ng = 2, mixture = ~ time + I(time^2) + I(time^3), B = m1.1.1, data = rundata),
  m3.1.1 = hlme(value ~ time + I(time^2) + I(time^3), random = ~ 1 + time, subject = "dah", ng = 3, mixture = ~ time + I(time^2) + I(time^3), B = m1.1.1, data = rundata),
  m4.1.1 = hlme(value ~ time + I(time^2) + I(time^3), random = ~ 1 + time, subject = "dah", ng = 4, mixture = ~ time + I(time^2) + I(time^3), B = m1.1.1, data = rundata),
  m5.1.1 = hlme(value ~ time + I(time^2) + I(time^3), random = ~ 1 + time, subject = "dah", ng = 5, mixture = ~ time + I(time^2) + I(time^3), B = m1.1.1, data = rundata)
)

#-----------------------------------------------
# Summarize Model Fit Statistics
#-----------------------------------------------

outcome <- data.frame(summarytable(
  m1, m1.1, m1.1.1,
  model_list$m2, model_list$m2.1, model_list$m2.1.1,
  model_list$m3, model_list$m3.1, model_list$m3.1.1,
  model_list$m4, model_list$m4.1, model_list$m4.1.1,
  model_list$m5, model_list$m5.1, model_list$m5.1.1,
  which = c("G", "loglik", "conv", "npm", "AIC", "BIC", "entropy", "%class")
))

outcome$`mean posterior probabilities` <- ""
outcome$`posterior probabilities > 0.7 (%)` <- ""
outcome$model <- row.names(outcome)

# Posterior probabilities summary (example for 2-class to 5-class models)
for (ng in 2:5) {
  for (suffix in c("", ".1", ".1.1")) {
    mod_name <- paste0("m", ng, suffix)
    mod <- model_list[[mod_name]]
    if (!is.null(mod)) {
      post1 <- data.frame(postprob(mod)[2])
      post2 <- data.frame(postprob(mod)[3])
      mean_probs <- paste0(post1[, 1:ng], collapse = "/")
      class70 <- paste0(post2[, 1:ng], collapse = "/")
      row_index <- which(outcome$G == ng & startsWith(outcome$model, mod_name))
      if (length(row_index) > 0) {
        outcome$`mean posterior probabilities`[row_index] <- mean_probs
        outcome$`posterior probabilities > 0.7 (%)`[row_index] <- class70
      }
    }
  }
}

# Save result table
write.xlsx(outcome, file = "output/trajectory_summary_RBC.xlsx", rowNames = FALSE)

#-----------------------------------------------
# Save Class Assignments and Trajectory Plots
#-----------------------------------------------

save_path_pprob <- "output/trajectory_class/"
save_path_plot  <- "output/trajectory_plot/"
dir.create(save_path_pprob, showWarnings = FALSE, recursive = TRUE)
dir.create(save_path_plot, showWarnings = FALSE, recursive = TRUE)

for (k in seq_along(model_list)) {
  mod <- model_list[[k]]
  mod_name <- names(model_list)[k]
  
  # Save class assignments
  if (!is.null(mod$pprob)) {
    class_data <- mod$pprob[, c("dah", "class")]
    fwrite(class_data, file = file.path(save_path_pprob, paste0("RBC_", mod_name, "_class.txt")), sep = "\t")
  }
  
  # Generate and save predicted trajectory plot
  time_seq <- data.frame(time = seq(min(rundata$time), max(rundata$time), length.out = 100))
  pred_data <- predictY(mod, newdata = time_seq, draw = FALSE)
  
  pred_df <- as.data.frame(pred_data$pred)
  pred_df$time <- time_seq$time
  pred_long <- pivot_longer(pred_df, cols = -time, names_to = "trajectory", values_to = "value")
  
  p <- ggplot(pred_long, aes(x = time, y = value, color = trajectory)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("#f16c23", "#2b6a99", "#1b7c3d", "#7f318d", "gray")) +
    labs(
      x = "Time to Diagnosis (Years)",
      y = "Log-Transformed Biomarker Value",
      title = paste0("RBC Trajectory - ", mod_name)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      axis.line = element_line(colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")
    )
  
  ggsave(filename = file.path(save_path_plot, paste0("RBC_", mod_name, ".jpeg")),
         plot = p, width = 8, height = 6, dpi = 300)
  
  print(paste0("Finished: RBC - ", mod_name))
}
