# Logistic and Ordinal Regression Analysis Script
# -----------------------------------------------
# Authors: Le Wang, Yifei Lin, Jin Huang
# Description: This script performs time-dependent logistic and ordinal logistic regression 
#              on blood cell indicators in patients with depressive disorder.

# Load Required Packages
library(data.table)
library(dplyr)
library(openxlsx)
library(glm2)
library(lubridate)
library(MASS)

#logistic regression####
# Load Data
data <- fread("data/data_1.csv", encoding = "UTF-8")
indicator <- read.xlsx("indicator_name.xlsx")

# Calculate Time Between Visits
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
# Logistic Regression: Unadjusted Model Example
#-----------------------------------------------
# Replace "Red Blood Cell Count" with the actual column name from your dataset
d <- data_transformed %>%
  select(
    dah,
    value = `Red Blood Cell Count`,  # <--- Replace this column name if necessary
    status_next,
    days_between_visits,
    gender,
    age,
    bmi,
    smoke,
    drink,
    diabetes,
    hypertension
  ) %>%
  filter(!value %in% c("", "/", "-----")) %>%
  mutate(value = scale(as.numeric(value))) %>%
  filter(!is.na(value))

fit <- glm(status_next ~ value + days_between_visits, family = binomial(link = "logit"), data = d)
summary(fit)

# Extract results
e <- coefficients(summary(fit))
confitOR <- exp(confint(fit))
p <- summary(fit)$coefficients[, 4]

estimate <- data.frame(
  estimate = e[, 1],
  OR = exp(e[, 1]),
  CI_L = confitOR[, 1],
  CI_U = confitOR[, 2],
  P = p
)
estimate$CI95 <- paste0(round(estimate$OR, 3), " (", round(estimate$CI_L, 3), ", ", round(estimate$CI_U, 3), ")")
estimate$indicator <- "Red Blood Cell Count"
estimate$flag <- row.names(estimate)
estimate <- estimate[c("indicator", "flag", "CI95", "P", "estimate", "OR", "CI_L", "CI_U")]

# Save results
write.xlsx(estimate, "output/1_logi_multi_time.xlsx", rowNames = FALSE)

#-----------------------------------------------
# Logistic Regression: Adjusted Models
#-----------------------------------------------

fit1 <- glm(status_next ~ value + days_between_visits + gender + age,
            family = binomial(link = "logit"), data = d)

fit2 <- glm(status_next ~ value + days_between_visits + gender + age + drink + smoke,
            family = binomial(link = "logit"), data = d)

fit3 <- glm(status_next ~ value + days_between_visits + gender + age + bmi + smoke + drink + diabetes + hypertension,
            family = binomial(link = "logit"), data = d)

#-----------------------------------------------
# Ordinal Logistic Regression Example
#-----------------------------------------------

fit_ordinal <- polr(class_next ~ value + days_between_visits + Covariates,
                    method = "logistic", data = d, Hess = TRUE)


