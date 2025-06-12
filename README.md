# 💊 Survival Analysis on VACS Dataset: HIV, SCD, and AMI Risk

This project analyzes the association between HIV status and adverse cardiovascular outcomes (sudden cardiac death [SCD] and acute myocardial infarction [AMI]) using data from the **Veterans Aging Cohort Study (VACS)**. It applies survival analysis methods to assess direct and indirect effects, including time-varying covariates and multiple imputation.

---

## 📚 Dataset Overview

- **Primary Dataset**: `vacs_scd_dat.csv`
- **Time-varying Dataset**: `timevar_cd4_final.csv`
- **AMI Dataset**: `recurrent_mi_vacs.csv`
- **Study Population**: Veterans with and without HIV
- **Primary Outcomes**: 
  - Sudden Cardiac Death (SCD)
  - Recurrent Acute Myocardial Infarction (AMI)

---

## 🧪 Analytical Methods

### ✅ Data Preprocessing
- Follow-up time calculated based on baseline and date of death (or December 31, 2014)
- Creation of CD4 and viral load (VL) categories
- Creation of race and ethnicity binary indicators
- Handling of missing data via multiple imputation using `mice`

### 📊 Descriptive Statistics
- Summary statistics stratified by HIV status
- Continuous variable summaries: mean, SD, median, IQR
- Categorical variable summaries: count, percentage

### ⚙️ Survival Analysis

#### Kaplan-Meier and Log-Rank
- Time to SCD stratified by HIV status
- Visualized using `survminer::ggsurvplot`

#### Poisson Regression
- Incidence rate per 100,000 person-years
- Models stratified by HIV, CD4, and VL categories

#### Cox Proportional Hazards
- Total effect and direct effect models
- Stratified models for HIV+ and HIV- groups
- Time-varying Cox model using start-stop format

### 🔄 Multiple Imputation
- Subsampled ~10,000 participants (HIV+ and HIV-)
- Imputed missing values for relevant predictors
- Pooled estimates from imputed datasets

## 📉 Key Outcomes

- **CD4 < 200 and VL ≥ 500** associated with increased SCD risk
- HIV+ status independently associated with elevated cardiovascular mortality
- Greater number of cardiovascular risk factors → higher SCD hazard
- Cox models indicate race, age, and clinical comorbidities as significant predictors
## 🔮 Future Work

- Extend to additional outcomes (e.g., stroke, all-cause mortality)
- Incorporate additional time-varying biomarkers (e.g., eGFR, BMI)
- Compare Cox models with machine learning survival methods (e.g., random survival forests)
- Create Shiny dashboard for exploring hazard ratios interactively



## 🧠 Key Learnings

- Applied survival analysis including Cox, KM, and Poisson models
- Managed missing data using multiple imputation (`mice`)
- Stratified model interpretations by HIV status, CD4, and VL
- Gained experience with spline modeling, time-varying covariates, and interaction terms



## 📬 Contact

**Madhav Dahal**  
📧 madhavdahal16@gmail.com  
🔗 [LinkedIn](https://www.linkedin.com/in/madhav-dahal-ms-9a1147b0)  
🔗 [GitHub](https://github.com/Madhav4487)
