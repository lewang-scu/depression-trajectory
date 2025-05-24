# Poisson Dispersion Analysis (Overall + Gender Stratified)
# ----------------------------------------------------------
# Author: Wang Le
# Description: Calculates Pearson-based dispersion statistics to assess overdispersion 
#              in Poisson models for blood cell indicators associated with depression.

# Load Required Packages
library(data.table)
library(dplyr)
library(openxlsx)
library(lubridate)
library(tidyr)

#-----------------------------------------------
# Load and Preprocess Data
#-----------------------------------------------
data <- fread("data.csv", encoding = "UTF-8")
indicator <- read.xlsx("indicator_name.xlsx", sheet = 1)

# Calculate time between visits
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
# Dispersion Calculation Function
#-----------------------------------------------
calc_dispersion <- function(fit) {
  pearson_chisq <- sum(residuals(fit, type = "pearson")^2)
  df <- df.residual(fit)
  pearson_chisq / df
}

#-----------------------------------------------
# Overall Dispersion Analysis
#-----------------------------------------------
dispersion_overall <- data.frame()

for (i in indicator$indicator) {
  d <- data_transformed %>%
    select(dah, value = all_of(i), status_next, days_between_visits, gender, age, bmi, smoke, drink, hypertension, diabetes) %>%
    filter(!value %in% c("", "/", "-----")) %>%
    mutate(value = scale(as.numeric(value))) %>%
    filter(!is.na(value))
  
  # Unadjusted
  fit0 <- glm(status_next ~ value + offset(log(days_between_visits)), family = poisson, data = d)
  dispersion_overall <- rbind(dispersion_overall, data.frame(indicator = i, model = "unadjusted model", dispersion = calc_dispersion(fit0)))
  
  # Adjusted model 1
  fit1 <- glm(status_next ~ value + offset(log(days_between_visits)) + gender + age, family = poisson, data = d)
  dispersion_overall <- rbind(dispersion_overall, data.frame(indicator = i, model = "adjusted model 1", dispersion = calc_dispersion(fit1)))
  
  # Adjusted model 2
  fit2 <- glm(status_next ~ value + offset(log(days_between_visits)) + gender + age + drink + smoke, family = poisson, data = d)
  dispersion_overall <- rbind(dispersion_overall, data.frame(indicator = i, model = "adjusted model 2", dispersion = calc_dispersion(fit2)))
  
  # Adjusted model 3
  fit3 <- glm(status_next ~ value + offset(log(days_between_visits)) + gender + age + bmi + smoke + drink + diabetes + hypertension, family = poisson, data = d)
  dispersion_overall <- rbind(dispersion_overall, data.frame(indicator = i, model = "adjusted model 3", dispersion = calc_dispersion(fit3)))
  
  print(paste0(i, " overall dispersion done"))
}

# Format wide table
df_overall <- dispersion_overall %>%
  pivot_wider(names_from = model, values_from = dispersion)

# Save
write.xlsx(df_overall, "output/dispersion_overall.xlsx", rowNames = FALSE)

#-----------------------------------------------
# Gender-Stratified Dispersion Analysis
#-----------------------------------------------
dispersion_gender <- data.frame()

for (i in indicator$indicator) {
  data_gender_split <- split(data_transformed, data_transformed$gender)
  
  for (gender_label in names(data_gender_split)) {
    d <- data_gender_split[[gender_label]] %>%
      select(dah, value = all_of(i), status_next, days_between_visits, age, bmi, smoke, drink, hypertension, diabetes) %>%
      filter(!value %in% c("", "/", "-----")) %>%
      mutate(value = scale(as.numeric(value))) %>%
      filter(!is.na(value))
    
    # Unadjusted
    fit0 <- glm(status_next ~ value + offset(log(days_between_visits)), family = poisson, data = d)
    dispersion_gender <- rbind(dispersion_gender, data.frame(
      indicator = i,
      model = "unadjusted model",
      gender = gender_label,
      dispersion = calc_dispersion(fit0)
    ))
    
    # Adjusted model 1
    fit1 <- glm(status_next ~ value + offset(log(days_between_visits)) + age, family = poisson, data = d)
    dispersion_gender <- rbind(dispersion_gender, data.frame(
      indicator = i,
      model = "adjusted model 1",
      gender = gender_label,
      dispersion = calc_dispersion(fit1)
    ))
    
    # Adjusted model 2
    fit2 <- glm(status_next ~ value + offset(log(days_between_visits)) + age + drink + smoke, family = poisson, data = d)
    dispersion_gender <- rbind(dispersion_gender, data.frame(
      indicator = i,
      model = "adjusted model 2",
      gender = gender_label,
      dispersion = calc_dispersion(fit2)
    ))
    
    # Adjusted model 3
    fit3 <- glm(status_next ~ value + offset(log(days_between_visits)) + age + bmi + smoke + drink + diabetes + hypertension, family = poisson, data = d)
    dispersion_gender <- rbind(dispersion_gender, data.frame(
      indicator = i,
      model = "adjusted model 3",
      gender = gender_label,
      dispersion = calc_dispersion(fit3)
    ))
    
    print(paste0(i, " - ", gender_label, " group dispersion done"))
  }
}

# Split and pivot
df_male <- dispersion_gender %>%
  filter(gender == "male") %>%
  pivot_wider(names_from = model, values_from = dispersion)

df_female <- dispersion_gender %>%
  filter(gender == "female") %>%
  pivot_wider(names_from = model, values_from = dispersion)

# Save
write.xlsx(df_male, "output/dispersion_male.xlsx", rowNames = FALSE)
write.xlsx(df_female, "output/dispersion_female.xlsx", rowNames = FALSE)
