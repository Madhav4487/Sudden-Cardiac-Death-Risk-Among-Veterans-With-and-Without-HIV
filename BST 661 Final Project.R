### Data preparation part one
# Set the working directory
setwd("/Users/ritabajgain/Library/CloudStorage/OneDrive-UniversityofKentucky/BST 661")

# Read in the dataset named “vacs_scd_dat” and make the following changes
vacs <- read.csv("vacs_scd_dat.csv", header = TRUE)

# View the data
View(vacs)

str(vacs) ##list of variables with their types

# Install and load the klibraries
install.packages("lubridate")
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(survminer)
library(mice)
library(splines)
library(rms)
library(ggplot2)

### 1.(a) Calculate the follow up time in years
# Convert baseline and date_death to date format using lubridate
vacs$baseline<-mdy(vacs$baseline)
vacs$date_death<-mdy(vacs$date_death)

# Define December 31, 2014 as a end date
end_date <- mdy("12/31/2014")

# Calculate end of follow-up time 
vacs$end_followup <- ifelse(
  is.na(vacs$date_death),             
  end_date,                             
  pmin(vacs$date_death, end_date, na.rm = TRUE) 
)
vacs$end_followup <- as.Date(vacs$end_followup, origin = "1970-01-01")

# Calculate follow-up time in days
vacs$followup_days <- as.numeric(vacs$end_followup - vacs$baseline)

# Convert follow-up time to years
vacs$followup_years <- vacs$followup_days / 365.25

# View the first few rows of the dataset
head(vacs[, c("baseline", "date_death", "end_followup", "followup_days", "followup_years")])

### 1(b).  Create categories of CD4 cell count
# Create a new variable for CD4 categories
vacs$CD4_category <- ifelse(
  vacs$hiv == 0, 0,                           # 0: Veterans without HIV
  ifelse(
    vacs$hiv == 1 & vacs$cd4 >= 500, 1,        # 1: HIV with CD4 >= 500
    ifelse(
      vacs$hiv == 1 & vacs$cd4 >= 200 & vacs$cd4 < 500, 2,  # 2: HIV with 200 <= CD4 < 500
      ifelse(
        vacs$hiv == 1 & vacs$cd4 < 200 & !is.na(vacs$cd4), 3, # 3: HIV with CD4 < 200
        NA                        # NA for other cases, such as missing values
      )
    )
  )
)

### 1(c) Create categories of HIV viral load (VL)
# Create a new variable for VL categories
vacs$VL_category <- ifelse(
  vacs$hiv == 0, 0,                               # 0: Veterans without HIV
  ifelse(
    vacs$hiv == 1 & !is.na(vacs$vl) & vacs$vl < 500, 1,  # 1: HIV with VL < 500
    ifelse(
      vacs$hiv == 1 & vacs$vl >= 500, 2,          # 2: HIV with VL >= 500
      NA                                          # NA for other cases
    )
  )
)

### 1(d). Create dummy variables (binary indicators) for white race Black race, Other race, and Hispanic ethnicity
vacs$Black <- ifelse(vacs$racecomg == 2, 1, 0)      # 1 = Black, 0 = Not Black
vacs$Hispanic <- ifelse(vacs$racecomg == 3, 1, 0)   # 1 = Hispanic, 0 = Not Hispanic
vacs$Other <- ifelse(vacs$racecomg == 4, 1, 0)      # 1 = Other , 0= Not Other

# (just checking here) Calculating the missing values for each variable in the dataset
missing_values <- sapply(vacs, function(x) sum(is.na(x)))

# Convert the result into a data frame for better readability
missing_summary <- data.frame(
  Variable = names(missing_values),
  Missing_Count = missing_values,
  Missing_Percentage = round((missing_values / nrow(vacs)) * 100, 2)  
)

# View the summary of missing values
print(missing_summary)

# 2. Count the number of participants in each HIV status category
hiv_counts <- vacs %>%
  dplyr::group_by(hiv) %>%
  dplyr::summarize(Count = dplyr::n())
print(hiv_counts)

# 2.  Calculate the mean, standard deviation, median and quantilkes of the continuous variables
summary_all <- vacs %>%
  group_by(hiv) %>%
  summarize(
    mean_age = round(mean(age, na.rm = TRUE), 1),
    sd_age = round(sd(age, na.rm = TRUE), 1),
    median_age = round(median(age, na.rm = TRUE), 1),
    Q1_age = round(quantile(age, 0.25, na.rm = TRUE), 1),
    Q3_age = round(quantile(age, 0.75, na.rm = TRUE), 1),
    
    mean_egfr = round(mean(egfr, na.rm = TRUE), 1),
    sd_egfr = round(sd(egfr, na.rm = TRUE), 1),
    median_egfr = round(median(egfr, na.rm = TRUE), 1),
    Q1_egfr = round(quantile(egfr, 0.25, na.rm = TRUE), 1),
    Q3_egfr = round(quantile(egfr, 0.75, na.rm = TRUE), 1),
    
    mean_bmi = round(mean(bmi, na.rm = TRUE), 1),
    sd_bmi = round(sd(bmi, na.rm = TRUE), 1),
    median_bmi = round(median(bmi, na.rm = TRUE), 1),
    Q1_bmi = round(quantile(bmi, 0.25, na.rm = TRUE), 1),
    Q3_bmi = round(quantile(bmi, 0.75, na.rm = TRUE), 1),
    
    mean_hgb = round(mean(hgb, na.rm = TRUE), 1),
    sd_hgb = round(sd(hgb, na.rm = TRUE), 1),
    median_hgb = round(median(hgb, na.rm = TRUE), 1),
    Q1_hgb = round(quantile(hgb, 0.25, na.rm = TRUE), 1),
    Q3_hgb = round(quantile(hgb, 0.75, na.rm = TRUE), 1)
  ) %>%
  pivot_longer(
    cols = -hiv, 
    names_to = c(".value", "Variable"),
    names_pattern = "(.*)_(.*)"
  )

# Print the summary table
print(summary_all)

# Calculate the count and percentage of the categorical variables stratified by HIV
calculate_categorical_summary <- function(data, categorical_var) {
  data %>%
    group_by(hiv, !!sym(categorical_var)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    group_by(hiv) %>%
    mutate(percentage = round(100 * count / sum(count), 1)) %>%
    arrange(hiv, !!sym(categorical_var))
}

# List of categorical variables to analyze
categorical_vars <- c("female", "htn", "VL_category", "CD4_category", "Black", 
                      "Hispanic", "Other", "prev_cvd", "dmcom", "dyslipidemia", 
                      "smk_stat", "alcohol", "cocaine", "hcv3", "copd", 
                      "regimen", "inc_scd")

# Apply the function to each categorical variable
summary_list <- map(categorical_vars, ~ calculate_categorical_summary(vacs, .x))

# Combine summaries into a named list for easy access
names(summary_list) <- categorical_vars

# Print all summaries
for (var in categorical_vars) {
  cat("\nSummary for", var, ":\n")
  print(summary_list[[var]])
}

# 3.  Fit the Kaplan-Meier curves for time to SCD stratified by HIV status
# Create a survival object
Y <- Surv(time = vacs$followup_years, event = vacs$inc_scd == 1)

# Fit Kaplan-Meier curves
kmfit <- survfit(Y ~ hiv, data = vacs)

# View the summary of the model
summary(kmfit)

# Perform the log-rank test
logrank_result <- survdiff(Y ~ hiv, data = vacs)

# Print the results
print(logrank_result)

# 4 this portion is only on the text file.

# 5.a.  The number of people with and without hiv developed by SCD
# Summarize the number of people with and without HIV and how many developed SCD
table_hiv_scd <- vacs %>%
  group_by(hiv, inc_scd) %>%
  summarize(Count = n(), .groups = 'drop')

# Print the table
print(table_hiv_scd)

# Summarize the number of SCD cases stratified by CD4 category for veterans with HIV
table_cd4_scd <- vacs %>%
  filter(hiv == 1) %>%    # Include only veterans with HIV
  group_by(CD4_category, inc_scd) %>%
  summarize(Count = n(), .groups = 'drop')

# Print the table
print(table_cd4_scd)

# 5. b. Summarize the number of SCD cases stratified by VL category for veterans with HIV
table_vl_scd <- vacs %>%
  filter(hiv == 1) %>%    # Include only veterans with HIV
  group_by(VL_category, inc_scd) %>%
  summarize(Count = n(), .groups = 'drop')

# Print the table
print(table_vl_scd)

# 6. Poission reression for the incidence rate and the confidence interval
# Step 1: Poisson Regression Model for HIV Status
poisson_hiv <- glm(inc_scd ~ hiv + age + female + Black + Hispanic + Other + offset(log(followup_years)),
                   family = poisson(link = "log"), data = vacs)

# Predict incidence rates for HIV status per 100,000 person-years
hiv_summary <- vacs %>%
  group_by(hiv) %>%
  summarize(
    total_scd = sum(inc_scd),
    total_py = sum(followup_years)
  ) %>%
  mutate(
    rate = (total_scd / total_py) * 100000,
    ci_lower = (total_scd - 1.96 * sqrt(total_scd)) / total_py * 100000,
    ci_upper = (total_scd + 1.96 * sqrt(total_scd)) / total_py * 100000
  )
print("Incidence Rate by HIV Status:")
print(hiv_summary)

# 6. a. Stratify by CD4 categories for HIV-positive veterans
vacs_cd4 <- vacs %>% filter(hiv == 1)
poisson_cd4 <- glm(inc_scd ~ CD4_category + age + female + Black + Hispanic + Other + offset(log(followup_years)),
                   family = poisson(link = "log"), data = vacs_cd4)

# Predict incidence rates by CD4 categories
cd4_summary <- vacs_cd4 %>%
  group_by(CD4_category) %>%
  summarize(
    total_scd = sum(inc_scd),
    total_py = sum(followup_years)
  ) %>%
  mutate(
    rate = (total_scd / total_py) * 100000,
    ci_lower = (total_scd - 1.96 * sqrt(total_scd)) / total_py * 100000,
    ci_upper = (total_scd + 1.96 * sqrt(total_scd)) / total_py * 100000
  )
print("Incidence Rate Stratified by CD4 Category:")
print(cd4_summary)

# 6. b. Stratify by VL categories for HIV-positive veterans
poisson_vl <- glm(inc_scd ~ VL_category + age + female + Black + Hispanic + Other + offset(log(followup_years)),
                  family = poisson(link = "log"), data = vacs_cd4)

# Predict incidence rates by VL categories
vl_summary <- vacs_cd4 %>%
  group_by(VL_category) %>%
  summarize(
    total_scd = sum(inc_scd),
    total_py = sum(followup_years)
  ) %>%
  mutate(
    rate = (total_scd / total_py) * 100000,
    ci_lower = (total_scd - 1.96 * sqrt(total_scd)) / total_py * 100000,
    ci_upper = (total_scd + 1.96 * sqrt(total_scd)) / total_py * 100000
  )
print("Incidence Rate Stratified by VL Category:")
print(vl_summary)

# 7.Assessment of the Proportional Hazards Assumption (Plot the log-log survival curves)
ggsurvplot(
  kmfit, 
  data = vacs,
  fun = "cloglog",       
  legend.title = "HIV Status", 
  legend.labs = c("HIV-", "HIV+"),
  xlab = "Time (log scale)", 
  ylab = "Log-Log Survival", 
  ggtheme = theme_minimal()
)

# 8.a. Calculating the missing values for each variable in the dataset
missing_values <- sapply(vacs, function(x) sum(is.na(x)))

# Convert the result into a data frame for better readability
missing_summary <- data.frame(
  Variable = names(missing_values),
  Missing_Count = missing_values,
  Missing_Percentage = round((missing_values / nrow(vacs)) * 100, 2)  
)

# View the summary of missing values
print(missing_summary)

# 8. b. Handling missing values
# running the given codes here. first 
set.seed(1)
hiv_sample_ids<-sample(vacs[vacs$hiv==1, 'new_id'], 3333, replace=F)
hiv_sample<-subset(vacs, vacs$new_id %in% hiv_sample_ids)
set.seed(2)
no_hiv_sample_ids<-sample(vacs[vacs$hiv==0, 'new_id'], 6667, replace=F)
no_hiv_sample<-subset(vacs, vacs$new_id %in% no_hiv_sample_ids)
vacs_subsample<-rbind(hiv_sample, no_hiv_sample)

# 8. c. Imputation using mice for the subsample dataset (vacs_subsample) created in b
# Function to count missing values
miss_vals <- function(x) { table(is.na(x)) }

# Apply the function to all columns in the dataset
lapply(vacs_subsample, miss_vals)

# Subset relevant variables for imputation
vacs_subsample_2imp <- vacs_subsample[, colnames(vacs_subsample) %in% c(
  "htn", "dyslipimedia", "smk_stat", "egfr", "bmi", "hgb", "hcv3", "dyslipidemia",
  "CD4_category", "VL_category", "followup_years", "age", "female", "hiv", "Black", "Hispanic", 
  "Other", "prev_cvd", "dmcom", "alcohol", "cocaine", "copd","inc_scd", "regimen")]

# set an initial imputation method structure.
ini <- mice(vacs_subsample_2imp, maxit = 0, seed = 1)
meth<-ini$method
pred<-ini$predictorMatrix

# Run the imputation
vacs_subsample_imp <- mice(
  vacs_subsample_2imp, 
  method = meth, 
  predictorMatrix = pred, 
  m = 5,
  seed = 82508
)

# 9.  Hazard ratio on the total effect model as a whole
mod5<-with(vacs_subsample_imp, coxph(Surv(followup_years, inc_scd)~hiv+age+female+Black+Hispanic+Other))
mod_summary<-summary(pool(mod5))
HRs<-exp(mod_summary[,2])
lower95<-exp(mod_summary[,2]-1.96*(mod_summary[,3]))
upper95<-exp(mod_summary[,2]+1.96*(mod_summary[,3]))
results<-data.frame(variable=mod_summary[,1], HRs, lower95, upper95,
                    pvalue=mod_summary[,6])
results

# 9. a. Convert CD4_category to a factor and set reference category
vacs_subsample_imp$data$CD4_category <- factor(vacs_subsample_imp$data$CD4_category)
vacs_subsample_imp$data$CD4_category <- relevel(vacs_subsample_imp$data$CD4_category, ref = "0")

# Fit Cox model with CD4_category
mod_cd4 <- with(vacs_subsample_imp, 
                coxph(Surv(followup_years, inc_scd) ~ CD4_category + age + female + Black + Hispanic + Other))

# Pool the results across imputed datasets
mod_cd4_summary <- summary(pool(mod_cd4))

# Extract hazard ratios and 95% CI for CD4 categories
HRs_cd4 <- exp(mod_cd4_summary[, 2])  # HRs
lower95_cd4 <- exp(mod_cd4_summary[, 2] - 1.96 * mod_cd4_summary[, 3])  # Lower CI
upper95_cd4 <- exp(mod_cd4_summary[, 2] + 1.96 * mod_cd4_summary[, 3])  # Upper CI

# Combine results into a clean data frame
results_cd4 <- data.frame(
  variable = mod_cd4_summary[, 1],
  HRs = HRs_cd4,
  lower95 = lower95_cd4,
  upper95 = upper95_cd4,
  pvalue = mod_cd4_summary[, 6]
)

# Print results
print("Hazard Ratios and 95% CIs for CD4 Categories:")
print(results_cd4)

# Convert VL_category to a factor and set reference category
vacs_subsample_imp$data$VL_category <- factor(vacs_subsample_imp$data$VL_category)
vacs_subsample_imp$data$VL_category <- relevel(vacs_subsample_imp$data$VL_category, ref = "0")

# Fit Cox model with VL_category
mod_vl <- with(vacs_subsample_imp, 
               coxph(Surv(followup_years, inc_scd) ~ VL_category + age + female + Black + Hispanic + Other))

# Pool the results across imputed datasets
mod_vl_summary <- summary(pool(mod_vl))

# Extract hazard ratios and 95% CI for VL categories
HRs_vl <- exp(mod_vl_summary[, 2])  # HRs
lower95_vl <- exp(mod_vl_summary[, 2] - 1.96 * mod_vl_summary[, 3])  # Lower CI
upper95_vl <- exp(mod_vl_summary[, 2] + 1.96 * mod_vl_summary[, 3])  # Upper CI

# Combine results into a clean data frame
results_vl <- data.frame(
  variable = mod_vl_summary[, 1],
  HRs = HRs_vl,
  lower95 = lower95_vl,
  upper95 = upper95_vl,
  pvalue = mod_vl_summary[, 6]
)

# Print results
print("Hazard Ratios and 95% CIs for VL Categories:")
print(results_vl)

# 9. b. Hazard ratio on the direct effect model as a whole
mod6<-with(vacs_subsample_imp, coxph(Surv(followup_years, inc_scd)~hiv+age+female+Black+Hispanic+Other+
                                       alcohol+ cocaine +htn+dmcom+dyslipidemia+egfr+bmi+hgb+prev_cvd+hcv3))
mod_summary<-summary(pool(mod6))
HRs<-exp(mod_summary[,2])
lower95<-exp(mod_summary[,2]-1.96*(mod_summary[,3]))
upper95<-exp(mod_summary[,2]+1.96*(mod_summary[,3]))
results<-data.frame(variable=mod_summary[,1], HRs, lower95, upper95,
                    pvalue=mod_summary[,6])
results

# Convert CD4_category to a factor and set reference category
vacs_subsample_imp$data$CD4_category <- factor(vacs_subsample_imp$data$CD4_category)
vacs_subsample_imp$data$CD4_category <- relevel(vacs_subsample_imp$data$CD4_category, ref = "0")

# Fit Cox model with CD4_category
mod1_cd4 <- with(vacs_subsample_imp, 
                coxph(Surv(followup_years, inc_scd) ~ CD4_category + age + female + Black + Hispanic + Other +
                        alcohol+ cocaine +htn+dmcom+dyslipidemia+egfr+bmi+hgb+prev_cvd+hcv3 ))

# Pool the results across imputed datasets
mod1_cd4_summary <- summary(pool(mod1_cd4))

# Extract hazard ratios and 95% CI for CD4 categories
HRs_cd4 <- exp(mod1_cd4_summary[, 2])  # HRs
lower95_cd4 <- exp(mod1_cd4_summary[, 2] - 1.96 * mod1_cd4_summary[, 3])  # Lower CI
upper95_cd4 <- exp(mod1_cd4_summary[, 2] + 1.96 * mod1_cd4_summary[, 3])  # Upper CI

# Combine results into a clean data frame
results1_cd4 <- data.frame(
  variable = mod1_cd4_summary[, 1],
  HRs = HRs_cd4,
  lower95 = lower95_cd4,
  upper95 = upper95_cd4,
  pvalue = mod1_cd4_summary[, 6]
)

# Print results
print("Hazard Ratios and 95% CIs for CD4 Categories:")
print(results1_cd4)

# Convert VL_category to a factor and set reference category
vacs_subsample_imp$data$VL_category <- factor(vacs_subsample_imp$data$VL_category)
vacs_subsample_imp$data$VL_category <- relevel(vacs_subsample_imp$data$VL_category, ref = "0")

# Fit Cox model with VL_category
mod2_vl <- with(vacs_subsample_imp, 
               coxph(Surv(followup_years, inc_scd) ~ VL_category + age + female + Black + Hispanic + Other + 
                       alcohol+ cocaine +htn+dmcom+dyslipidemia+egfr+bmi+hgb+prev_cvd+hcv3 ))

# Pool the results across imputed datasets
mod2_vl_summary <- summary(pool(mod2_vl))

# Extract hazard ratios and 95% CI for VL categories
HRs_vl <- exp(mod2_vl_summary[, 2])  # HRs
lower95_vl <- exp(mod2_vl_summary[, 2] - 1.96 * mod2_vl_summary[, 3])  # Lower CI
upper95_vl <- exp(mod2_vl_summary[, 2] + 1.96 * mod2_vl_summary[, 3])  # Upper CI

# Combine results into a clean data frame
results2_vl <- data.frame(
  variable = mod2_vl_summary[, 1],
  HRs = HRs_vl,
  lower95 = lower95_vl,
  upper95 = upper95_vl,
  pvalue = mod2_vl_summary[, 6]
)

# Print results
print("Hazard Ratios and 95% CIs for VL Categories:")
print(results2_vl)


# QN 10, import the another dataset named “timevar_cd4_final”
timevar <- read.csv("timevar_cd4_final.csv", header =TRUE)
# View the dataset
View(timevar)

# 10. a. Create a new variable for CD4 categories
timevar$CD4_category <- ifelse(
  timevar$hiv == 0, 0,                           # 0: Veterans without HIV
  ifelse(
    timevar$hiv == 1 & timevar$cd4_count >= 500, 1,        # 1: HIV with CD4 >= 500
    ifelse(
      timevar$hiv == 1 & timevar$cd4_count >= 200 & timevar$cd4_count < 500, 2,  # 2: HIV with 200 <= CD4 < 500
      ifelse(
        timevar$hiv == 1 & timevar$cd4_count < 200 & !is.na(timevar$cd4_count), 3, # 3: HIV with CD4 < 200
        NA                        # NA for other cases, such as missing values
      )
    )
  )
)

# Create dummy variables (binary indicators) for white race, Black race, Other race, and Hispanic ethnicity
timevar$Black <- ifelse(timevar$racecomg == 2, 1, 0)      # 1 = Black, 0 = Not Black
timevar$Hispanic <- ifelse(timevar$racecomg == 3, 1, 0)   # 1 = Hispanic, 0 = Not Hispanic
timevar$Other <- ifelse(timevar$racecomg == 4, 1, 0)      # 1 = Other , 0= Not Other

# 10. b. Excluding the missing cd4_count (if HIV positive) as well as the condounders (age, female, Black, Other)
# Exclude rows with missing cd4_count for HIV-positive individuals
timevar <- timevar[!(timevar$hiv == 1 & is.na(timevar$CD4_category)), ] # cd4_count is renamed as CD4_ccategory here

# Specify the confounders and mediators
confounders <- c("age", "female", "Black", "Hispanic", "Other", "prev_cvd", "htn", 
                 "dmcom", "dyslipidemia", "egfr", "bmi", "hgb", "alcohol", "cocaine", 
                 "hcv3")
                   
# Exclude rows with missing values in confounders and mediators
timevar <- timevar[complete.cases(timevar[, confounders]), ]

# 11. Fit a time-varying Cox proportional hazards model to calculate the direct effect
# of HIV stratified by CD4_category on time to SCD.
# 11. a. # Convert CD4_category to factor and set a reference level
timevar$CD4_category <- factor(timevar$CD4_category)
timevar$CD4_category <- relevel(timevar$CD4_category, ref = "0")

# Fit Cox model stratified by CD4_category
cox_stratified_cd4 <- coxph(Surv(tstart, tstop, scd) ~ CD4_category+age+female+Black+Hispanic+Other+
                              alcohol+ cocaine +htn+dmcom+dyslipidemia+egfr+bmi+hgb+prev_cvd+hcv3, 
                            data = timevar)

# Summarize the model
summary_stratified_cd4 <- summary(cox_stratified_cd4)

# Print model summary
print(summary_stratified_cd4)


# 12. Assessment of Interactions
# 12.a. Fit a model that includes the main effects of HIV as well as all confounders and mediators, and 2-way interactions 
# between HIV and each of the confounders and mediators
# List of confounders
confounders <- c("age", "female", "Black", "Hispanic", "Other")

# List of mediators
mediators <- c("prev_cvd", "htn", "dmcom", "dyslipidemia", "alcohol", "cocaine", "bmi", "egfr", "hcv3", "hgb")

# Combine confounders and mediators
all_vars <- c(confounders, mediators)

# Generate interaction terms with HIV
interaction_terms <- paste("hiv *", all_vars, collapse = " + ")

# Create the full formula
formula_a <- as.formula(paste("Surv(followup_years, inc_scd) ~ hiv +", paste(all_vars, collapse = " + "), "+", interaction_terms))

# Fit the Cox model
model_a <- coxph(formula_a, data = vacs)

# View summary
summary(model_a)

# 12. b. Fit a Model with Main Effects of HIV, Confounders, and Mediators
# Create the formula for main effects only
formula_b <- as.formula(paste("Surv(followup_years, inc_scd) ~ hiv +", paste(all_vars, collapse = " + ")))

# Fit the Cox model
model_b <- coxph(formula_b, data = vacs)

# View summary
summary(model_b)

# 12. c. Perform a Chunk Test for the 2-Way Interactions
# Perform the chunk test
chunk_test <- anova(model_b, model_a, test = "LRT")

# View the results
print(chunk_test)

# Determine significance based on alpha = 0.15
if (chunk_test$`Pr(>|Chi|)`[2] < 0.15) {
  cat("The 2-way interactions are significantly associated with time to SCD (p <", chunk_test$`Pr(>|Chi|)`[2], ").\n")
} else {
  cat("The 2-way interactions are not significantly associated with time to SCD (p =", chunk_test$`Pr(>|Chi|)`[2], ").\n")
}

# d. ii. (12 i is skipped here)
# Extract the first imputed dataset
vacs_complete <- complete(vacs_subsample_imp, 1)

# Subset for veterans with HIV
with_hiv <- vacs_complete[vacs_complete$hiv == 1, ]

# Subset for veterans without HIV
without_hiv <- vacs_complete[vacs_complete$hiv == 0, ]

# Convcert these into factor on the dataset called with_HIV
# Convert hypertension  to a factor with none as the reference level
with_hiv$htn <- factor(with_hiv$htn, 
                                   levels = c(0,1, 2), 
                                   labels = c("None", "Controlled", "Uncontrolled"))

# Convert smoking status to a factor with Never as the reference level
with_hiv$smk_stat <- factor(with_hiv$smk_stat, 
                       levels = c(0,1, 2), 
                       labels = c("Never", "Current", "Former"))

# Convert Hepatitis to a factor with Negative as the reference level
with_hiv$hcv3 <- factor(with_hiv$hcv3, 
                            levels = c(0,1,2), 
                            labels = c("Negative", "Postivie", "Never Tested"))


# Model for veterans with HIV
cox_with_hiv<- coxph(Surv(followup_years, inc_scd) ~ 
                          age +female+Black+Hispanic+Other+ prev_cvd + htn + dmcom + dyslipidemia + smk_stat+
                          hcv3 + egfr + bmi + hgb + alcohol + cocaine + copd,
                        data = with_hiv)

summary(cox_with_hiv)

# Convcert these into factor on the dataset called without _HIV
# Convert hypertension  to a factor with none as the reference level
without_hiv$htn <- factor(without_hiv$htn, 
                       levels = c(0,1, 2), 
                       labels = c("None", "Controlled", "Uncontrolled"))

# Convert smoking status to a factor with Never as the reference level
without_hiv$smk_stat <- factor(without_hiv$smk_stat, 
                            levels = c(0,1, 2), 
                            labels = c("Never", "Current", "Former"))

# Convert Hepatitis to a factor with Negative as the reference level
without_hiv$hcv3 <- factor(without_hiv$hcv3, 
                        levels = c(0,1,2), 
                        labels = c("Negative", "Postivie", "Never Tested"))


# Model for veterans without  HIV
cox_without_hiv<- coxph(Surv(followup_years, inc_scd) ~ 
                        age +female+Black+Hispanic+Other+ prev_cvd + htn + dmcom + dyslipidemia + smk_stat+
                        hcv3 + egfr + bmi + hgb + alcohol + cocaine + copd,
                      data = without_hiv)

summary(cox_without_hiv)


# 13. Association of Number of CVD Risk Factors with time to SCD by HIV status
# a.Create binary indicators of the following conditions (1 if present, 0 if absent)
# Hypertension indicator
vacs$htn_present <- ifelse(vacs$htn %in% c(1, 2), 1, 0)

# Current smoking indicator
vacs$current_smoker <- ifelse(vacs$smk_stat == 1, 1, 0)

# Hepatitis C virus indicator
vacs$hcv_positive <- ifelse(vacs$hcv3 == 1, 1, 0)

# Low hemoglobin indicator
vacs$low_hgb <- ifelse(vacs$hgb < 10, 1, 0)

# 13.b. # Create cv_rfs variable
vacs$cv_rfs <- rowSums(
  vacs[, c("prev_cvd", "alcohol", "copd", "htn_present", "current_smoker", "hcv_positive", "low_hgb")], 
  na.rm = FALSE)

# 13.c.First, # Subset data by HIV status
vacs_hiv <- subset(vacs, hiv == 1)         # Veterans with HIV
vacs_no_hiv <- subset(vacs, hiv == 0)      # Veterans without HIV

# Calculate percentiles for cv_rfs for veterans with HIV
knots_hiv <- quantile(vacs_hiv$cv_rfs, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)

# Calculate percentiles for cv_rfs for veterans without HIV
knots_no_hiv <- quantile(vacs_no_hiv$cv_rfs, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)

# Fit the Cox model for veterans with HIV
cox_hiv <- coxph(Surv(followup_years, inc_scd) ~ ns(cv_rfs, knots = knots_hiv) +
                   age + female + Black + Hispanic + Other,
                 data = vacs_hiv)

# Summary of the model
summary(cox_hiv)

# Fit the Cox model for veterans without HIV
cox_no_hiv <- coxph(Surv(followup_years, inc_scd) ~ ns(cv_rfs, knots = knots_no_hiv) +
                      age + female + Black + Hispanic + Other,
                    data = vacs_no_hiv)

# Summary of the model
summary(cox_no_hiv)

# 13. d. The knots are defined as follows
# Subset data for veterans with HIV
vacs_hiv <- subset(vacs, hiv == 1)

# Set up data distribution
dd_hiv <- datadist(vacs_hiv)
options(datadist = 'dd_hiv')

# Determine knot locations for cv_rfs at 10th, 50th, and 90th percentiles
knots_hiv <- quantile(vacs_hiv$cv_rfs, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
print("Knots for cv_rfs (HIV+):")
print(knots_hiv)

# Subset data for veterans without HIV
vacs_no_hiv <- subset(vacs, hiv == 0)

# Set up data distribution
dd_no_hiv <- datadist(vacs_no_hiv)
options(datadist = 'dd_no_hiv')

# Determine knot locations for cv_rfs at 10th, 50th, and 90th percentiles
knots_no_hiv <- quantile(vacs_no_hiv$cv_rfs, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
print("Knots for cv_rfs (HIV-):")
print(knots_no_hiv)

# Fit the Cox PH model with restricted cubic spline for cv_rfs
mod_hiv <- cph(Surv(followup_years, inc_scd) ~ rcs(cv_rfs, knots = knots_hiv) + 
                 age + female + Black + Hispanic + Other,
               data = vacs_hiv, x = TRUE, y = TRUE)

# Print the model summary
print("Cox PH Model with RCS for Veterans with HIV:")
print(mod_hiv)

# Fit the Cox PH model with restricted cubic spline for cv_rfs
mod_no_hiv <- cph(Surv(followup_years, inc_scd) ~ rcs(cv_rfs, knots = knots_no_hiv) + 
                    age + female + Black + Hispanic + Other,
                  data = vacs_no_hiv, x = TRUE, y = TRUE)

# Print the model summary
print("Cox PH Model with RCS for Veterans without HIV:")
print(mod_no_hiv)

# Plot the spline effect of cv_rfs for veterans with HIV
plot(Predict(mod_hiv, cv_rfs), 
     main = "Effect of CV Risk Factors on SCD Risk (HIV+)",
     xlab = "Number of Cardiovascular Risk Factors (cv_rfs)",
     ylab = "Hazard Ratio (HR)",
     col = "purple")

# Plot the spline effect of cv_rfs for veterans without HIV
plot(Predict(mod_no_hiv, cv_rfs), 
     main = "Effect of CV Risk Factors on SCD Risk (HIV-)",
     xlab = "Number of Cardiovascular Risk Factors (cv_rfs)",
     ylab = "Hazard Ratio (HR)",
     col = "red")

# 14. Read the data named "recurrent_mi_vacs.csv"
recurrent <- read.csv("recurrent_mi_vacs.csv", header = TRUE)
# View the data
View(recurrent)

# Association of HIV with recurrent acute myocardial infarction (AMI)
# a. # Count total AMI events
total_ami_events <- sum(recurrent$ami, na.rm = TRUE)

# Print the result
cat("Total AMI events:", total_ami_events, "\n")

# b. # Aggregate the number of AMI events per individual
ami_counts <- aggregate(ami ~ new_id, data = recurrent, sum)

# Count individuals with more than 1 AMI event
individuals_with_multiple_ami <- sum(ami_counts$ami > 1)

# Print the result
cat("Number of people with more than 1 AMI:", individuals_with_multiple_ami, "\n")

# Check for missing values
sum(is.na(recurrent$tstart))
sum(is.na(recurrent$tstop))

# Identify problematic rows
problematic_rows <- recurrent[recurrent$tstop <= recurrent$tstart, ]

# View the problematic rows
print(problematic_rows)

# Exclude rows with tstop <= tstart
recurrent_clean <- recurrent[recurrent$tstop > recurrent$tstart, ]

# 14.c.  Fit the Cox proportional hazards model
# Convert racecomg to a factor with White as the reference level
recurrent_clean$racecomg <- factor(recurrent_clean$racecomg, 
                                   levels = c(1, 2, 3, 4), 
                                   labels = c("White", "Black", "Hispanic", "Other"))

recurrent_cox_model <- coxph(Surv(tstart, tstop, ami== 1) ~ hiv + age_bl + female + racecomg, 
                   data = recurrent_clean)
summary(recurrent_cox_model)

mod_summary = summary(recurrent_cox_model)

# Extract coefficients, standard errors, and p-values
coef_col <- "coef"
se_col <- "se(coef)"
pval_col <- "Pr(>|z|)"

HRs <- exp(mod_summary$coefficients[, coef_col])
lower95 <- exp(mod_summary$coefficients[, coef_col] - 1.96 * mod_summary$coefficients[, se_col])
upper95 <- exp(mod_summary$coefficients[, coef_col] + 1.96 * mod_summary$coefficients[, se_col])
pvalues <- mod_summary$coefficients[, pval_col]

# Create results data frame
results <- data.frame(
  variable = rownames(mod_summary$coefficients),
  HRs = HRs,
  lower95 = lower95,
  upper95 = upper95,
  pvalue = pvalues
)

# View results
print(results)