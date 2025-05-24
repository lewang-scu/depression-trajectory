# Poisson Regression with Robust Standard Errors
# ---------------------------------------------
# Author: Wang Le (modified for robust SE)
# Description: This script performs Poisson regression using robust SE
#              on time-dependent data to analyze the association between
#              blood cell indicators and depressive outcomes.

# Load Required Packages
library(data.table)
library(dplyr)
library(openxlsx)
library(lubridate)
library(sandwich)
library(lmtest)

# Load dataset
data <- fread("data/data_1.csv", encoding = "UTF-8")

# Load indicator list
indicator <- read.xlsx("data/indicator_name.xlsx", sheet = 1)

# Compute time between visits
data_transformed <- data %>%
  arrange(dah, jcsj) %>%
  group_by(dah) %>%
  mutate(
    status_next = lead(status),
    jcsj_next = lead(jcsj),
    days_between_visits = as.numeric(difftime(jcsj_next, jcsj, units = "days"))
  ) %>%
  ungroup() %>%
  filter(!is.na(status_next))

#-----------------------------------------------
# Example: Poisson Model with Robust SE
#-----------------------------------------------
d <- data_transformed %>%
  select(dah, value = `Red Blood Cell Count`, status_next, days_between_visits) %>%
  filter(!value %in% c("", "/", "-----")) %>%
  mutate(value = scale(as.numeric(value))) %>%
  filter(!is.na(value))

fit <- glm(status_next ~ value + offset(log(days_between_visits)), family = poisson, data = d)
robust_result <- coeftest(fit, vcovCL(fit, cluster = ~dah))

est <- robust_result[, 1]
se <- robust_result[, 2]
CI_L <- est - 1.96 * se
CI_U <- est + 1.96 * se
OR <- exp(est)

estimate <- data.frame(
  estimate = est,
  OR = OR,
  CI_L = exp(CI_L),
  CI_U = exp(CI_U),
  P = robust_result[, 4],
  CI95 = paste0(round(OR, 3), " (", round(exp(CI_L), 3), ", ", round(exp(CI_U), 3), ")"),
  indicator = "Red Blood Cell Count",
  flag = rownames(robust_result)
)
estimate <- estimate[c("indicator", "flag", "CI95", "P", "estimate", "OR", "CI_L", "CI_U")]
out <- estimate[1, ]  # keep value row only

#-----------------------------------------------
# Unadjusted Poisson Models with Robust SE
#-----------------------------------------------
for (i in indicator$indicator) {
  d <- data_transformed %>%
    select(dah, value = all_of(i), status_next, days_between_visits) %>%
    filter(!value %in% c("", "/", "-----")) %>%
    mutate(value = scale(as.numeric(value))) %>%
    filter(!is.na(value))
  
  fit <- glm(status_next ~ value + offset(log(days_between_visits)), family = poisson, data = d)
  robust_result <- coeftest(fit, vcovCL(fit, cluster = ~dah))
  
  est <- robust_result[, 1]
  se <- robust_result[, 2]
  CI_L <- est - 1.96 * se
  CI_U <- est + 1.96 * se
  OR <- exp(est)
  
  estimate <- data.frame(
    estimate = est,
    OR = OR,
    CI_L = exp(CI_L),
    CI_U = exp(CI_U),
    P = robust_result[, 4],
    CI95 = paste0(round(OR, 3), " (", round(exp(CI_L), 3), ", ", round(exp(CI_U), 3), ")"),
    indicator = i,
    flag = rownames(robust_result)
  )
  estimate <- estimate[c("indicator", "flag", "CI95", "P", "estimate", "OR", "CI_L", "CI_U")]
  out <- rbind(out, subset(estimate, flag == "value"))
  print(paste0(i, " finished"))
}

out$FDR <- p.adjust(out$P, method = "fdr")
write.xlsx(out, file = "output/1_poisson.xlsx", rowNames = FALSE)

#-----------------------------------------------
# Adjusted Poisson Models (Model 1–3) with Robust SE
#-----------------------------------------------
out <- data.frame()

for (i in indicator$indicator) {
  d <- data_transformed %>%
    select(
      dah, value = all_of(i), days_between_visits, status_next,
      age, gender, bmi, smoke, drink, hypertension, diabetes
    ) %>%
    filter(!value %in% c("", "/", "-----")) %>%
    mutate(value = scale(as.numeric(value))) %>%
    filter(!is.na(value))
  
  model_list <- list(
    "adjusted model 1" = status_next ~ value + offset(log(days_between_visits)) + gender + age,
    "adjusted model 2" = status_next ~ value + offset(log(days_between_visits)) + gender + age + drink + smoke,
    "adjusted model 3" = status_next ~ value + offset(log(days_between_visits)) + gender + age + bmi + smoke + drink + diabetes + hypertension
  )
  
  for (model_name in names(model_list)) {
    fit <- glm(model_list[[model_name]], family = poisson, data = d)
    robust_result <- coeftest(fit, vcovCL(fit, cluster = ~dah))
    
    est <- robust_result[, 1]
    se <- robust_result[, 2]
    CI_L <- est - 1.96 * se
    CI_U <- est + 1.96 * se
    OR <- exp(est)
    
    estimate <- data.frame(
      estimate = est,
      OR = OR,
      CI_L = exp(CI_L),
      CI_U = exp(CI_U),
      P = robust_result[, 4],
      CI95 = paste0(round(OR, 3), " (", round(exp(CI_L), 3), ", ", round(exp(CI_U), 3), ")"),
      indicator = i,
      flag = rownames(robust_result),
      model = model_name
    )
    out <- rbind(out, subset(estimate, flag == "value"))
  }
  print(paste0(i, " finished"))
}

out$FDR <- p.adjust(out$P, method = "fdr")
write.xlsx(out, file = "output/1_poisson_multi.xlsx", rowNames = FALSE)
