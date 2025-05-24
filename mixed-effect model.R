# Mixed Effects Logistic Regression for Biomarker and Depression Status
# ----------------------------------------------------------------------
# Author: Wang Le
# Description: This script performs mixed effects logistic regression to evaluate 
#              the association between blood cell indicators and depression status,
#              accounting for repeated measures per subject.

# Load Required Packages
library(data.table)
library(dplyr)
library(openxlsx)
library(lme4)
library(broom.mixed)

#-----------------------------------------------
# Load Data
#-----------------------------------------------

# Load main dataset and biomarker list
data <- fread("data.csv", encoding = "UTF-8")
indicator <- read.xlsx("indicator_name.xlsx", sheet = 1)

# Initialize result storage
out <- data.frame(
  lag = "", indicator = "", OR = "", CI95 = "", P = "",
  Estimate = "", SE = "", Z = "", CI_L = "", CI_U = "",
  model = "", converge = ""
)

#-----------------------------------------------
# Run Mixed Effects Models (unadjusted + adjusted)
#-----------------------------------------------

set.seed(128)

for (i in indicator$indicator) {
  
  d <- data %>%
    select(dah, status, gender, age, smoke, drink, hypertension, diabetes, bmi, value = all_of(i)) %>%
    filter(!value %in% c("", "/", ".", "-----", ",")) %>%
    mutate(value = as.numeric(value)) %>%
    filter(!is.na(value)) %>%
    mutate(value = scale(value))
  
  # Fit unadjusted model
  fit <- glmer(status ~ value + (1 | dah), data = d, family = binomial,
               control = glmerControl(optimizer = "bobyqa"), nAGQ = 20)
  
  estimate <- data.frame(summary(fit)[["coefficients"]])
  estimate$lag <- row.names(estimate)
  estimate$OR <- exp(estimate$Estimate)
  estimate$CI_L <- exp(estimate$Estimate - 1.96 * estimate$Std..Error)
  estimate$CI_U <- exp(estimate$Estimate + 1.96 * estimate$Std..Error)
  estimate$CI95 <- paste0(round(estimate$OR, 3), " (", round(estimate$CI_L, 3), ", ", round(estimate$CI_U, 3), ")")
  estimate$indicator <- i
  estimate$P <- estimate$`Pr(>|z|)`
  estimate <- estimate[, c("lag", "indicator", "OR", "CI95", "P", "Estimate", "Std..Error", "z value", "CI_L", "CI_U")]
  colnames(estimate)[6:8] <- c("Estimate", "SE", "Z")
  estimate$model <- "unadjusted"
  estimate$converge <- ifelse(length(fit@optinfo$conv$lme4) > 0, "no", "yes")
  out <- rbind(out, estimate)
  print(paste0(i, " - unadjusted model complete"))
  
  # Adjusted Model 1: + gender + age
  fit1 <- glmer(status ~ value + gender + age + (1 | dah), data = d, family = binomial,
                control = glmerControl(optimizer = "bobyqa"), nAGQ = 20)
  
  estimate <- data.frame(summary(fit1)[["coefficients"]])
  estimate$lag <- row.names(estimate)
  estimate$OR <- exp(estimate$Estimate)
  estimate$CI_L <- exp(estimate$Estimate - 1.96 * estimate$Std..Error)
  estimate$CI_U <- exp(estimate$Estimate + 1.96 * estimate$Std..Error)
  estimate$CI95 <- paste0(round(estimate$OR, 3), " (", round(estimate$CI_L, 3), ", ", round(estimate$CI_U, 3), ")")
  estimate$indicator <- i
  estimate$P <- estimate$`Pr(>|z|)`
  estimate <- estimate[, c("lag", "indicator", "OR", "CI95", "P", "Estimate", "Std..Error", "z value", "CI_L", "CI_U")]
  colnames(estimate)[6:8] <- c("Estimate", "SE", "Z")
  estimate$model <- "adjusted model 1"
  estimate$converge <- ifelse(length(fit1@optinfo$conv$lme4) > 0, "no", "yes")
  out <- rbind(out, estimate)
  print(paste0(i, " - model 1 complete"))
  
  # Adjusted Model 2: + smoke + drink
  fit2 <- glmer(status ~ value + gender + age + smoke + drink + (1 | dah), data = d, family = binomial,
                control = glmerControl(optimizer = "bobyqa"), nAGQ = 20)
  
  estimate <- data.frame(summary(fit2)[["coefficients"]])
  estimate$lag <- row.names(estimate)
  estimate$OR <- exp(estimate$Estimate)
  estimate$CI_L <- exp(estimate$Estimate - 1.96 * estimate$Std..Error)
  estimate$CI_U <- exp(estimate$Estimate + 1.96 * estimate$Std..Error)
  estimate$CI95 <- paste0(round(estimate$OR, 3), " (", round(estimate$CI_L, 3), ", ", round(estimate$CI_U, 3), ")")
  estimate$indicator <- i
  estimate$P <- estimate$`Pr(>|z|)`
  estimate <- estimate[, c("lag", "indicator", "OR", "CI95", "P", "Estimate", "Std..Error", "z value", "CI_L", "CI_U")]
  colnames(estimate)[6:8] <- c("Estimate", "SE", "Z")
  estimate$model <- "adjusted model 2"
  estimate$converge <- ifelse(length(fit2@optinfo$conv$lme4) > 0, "no", "yes")
  out <- rbind(out, estimate)
  print(paste0(i, " - model 2 complete"))
  
  # Adjusted Model 3: + bmi + comorbidities
  fit3 <- glmer(status ~ value + gender + age + bmi + smoke + drink + diabetes + hypertension + (1 | dah), 
                data = d, family = binomial,
                control = glmerControl(optimizer = "bobyqa"), nAGQ = 20)
  
  estimate <- data.frame(summary(fit3)[["coefficients"]])
  estimate$lag <- row.names(estimate)
  estimate$OR <- exp(estimate$Estimate)
  estimate$CI_L <- exp(estimate$Estimate - 1.96 * estimate$Std..Error)
  estimate$CI_U <- exp(estimate$Estimate + 1.96 * estimate$Std..Error)
  estimate$CI95 <- paste0(round(estimate$OR, 3), " (", round(estimate$CI_L, 3), ", ", round(estimate$CI_U, 3), ")")
  estimate$indicator <- i
  estimate$P <- estimate$`Pr(>|z|)`
  estimate <- estimate[, c("lag", "indicator", "OR", "CI95", "P", "Estimate", "Std..Error", "z value", "CI_L", "CI_U")]
  colnames(estimate)[6:8] <- c("Estimate", "SE", "Z")
  estimate$model <- "adjusted model 3"
  estimate$converge <- ifelse(length(fit3@optinfo$conv$lme4) > 0, "no", "yes")
  out <- rbind(out, estimate)
  
  print(paste0(i, " - model 3 complete"))
}

#-----------------------------------------------
# Save Results
#-----------------------------------------------
write.xlsx(out, file = "output/mixed_model_results_ratio.xlsx", rowNames = FALSE)
