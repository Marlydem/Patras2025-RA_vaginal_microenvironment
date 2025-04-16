setwd("/Users/marlydmejia/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Desktop\ Documents-\ Patras\ Lab/Proposals/Manuscripts/Mine/2024_Rheumatoid_Arthritis_Vaginal/VSP80_2024updated")

install.packages("gtsummary")
install.packages("ggpubr")
library(gtsummary)
install.packages(c("dplyr"))
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr) #for cytokine boxplots

install.packages("tidyverse")
install.packages("car")
install.packages("lmtest")
install.packages("nlme")
install.packages("ggridges")
library(tidyverse)
library(car)  # For variance inflation factors (VIFs) to check multicollinearity
library(lmtest)  # For model diagnostics
library(ggridges)
install.packages(c("dplyr", "ggplot2"))
install.packages("pROC")
install.packages("MASS")
install.packages("emmeans")
library(pROC)
library(nlme)
library(MASS)
library(emmeans)


# Install and load necessary packages
install.packages(c("boot"))
library(boot)
your_data_few <- read.csv("Metadata_RAstudy_master_categorized_noNANreally.csv") #limited adjust. variables for microbes
your_data_noNAN <- read.csv("Metadata_RAstudy_master_categorized_noNAN.csv") #alt noNAN
your_data_all <- read.csv("Metadata_RAstudy_master_categorized_v2excess.csv") #full
quartz()


##### Getting Table data as factors, ordering, and renaming #####
your_data_all$Menopausal_bin <- factor(your_data_all$Menopausal_bin,
                                  levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
#your_data_all$Diet_bin <- factor(your_data_all$Diet_bin,
#                                   levels = c("High", "Reduced", "Not disclosed"))
your_data_all$Diet_bin <- factor(your_data_all$Diet_bin,
                                  levels = c("1", "0", "2"))
your_data_all$Diet_bin <- revalue(your_data_all$Diet_bin, c("0" = "Reduced",
                                                                "1" = "High",
                                                                "2" = "Not disclosed"))
your_data_all$Alcohol_consumption <- factor(your_data_all$Alcohol_consumption,
                                  levels = c("Once daily", "Occasionally", "Rarely", "Not applicable", "Not disclosed"))
your_data_all$Cigarette_smoking <- factor(your_data_all$Cigarette_smoking,
                                  levels = c("Multiple times a day", "Occasionally", "Rarely", "Not applicable", "Not disclosed"))
#your_data_all$Education <- factor(your_data_all$Education,
#            levels = c("Advanced degree","College/University degree",  "Some college", "High School", "Elementary/Middle #School", "Not disclosed"))
your_data_all$Education <- factor(your_data_all$Education,
                                  levels = c("Advanced Degree","College/University Degree",  "Some College", "High School", "Elementary/Middle School", "Not disclosed"))
your_data_all$Education <- revalue(your_data_all$Education, 
                                   c("Advanced Degree" = "Advanced degree",
                                    "College/University Degree" = "College/University degree",
                                    "Some College" = "Some college",
                                    "High School" = "High school",
                                    "Elementary/Middle School" = "Elementary/Middle school",
                                    "Not disclosed" = "Not disclosed"))
your_data_all$Pregnancies <- factor(your_data_all$Pregnancies,
                                                     levels = c("2+", "1", "None", "Not disclosed"))
your_data_all$Pregnancies <- revalue(your_data_all$Pregnancies, 
                                                      c("None" = "0",
                                                        "1" = "1",
                                                        "2+" = "2+",
                                                        "Not disclosed" = "Not disclosed"))
your_data_all$Lifetime_Sexual_Partners_bin <- factor(your_data_all$Lifetime_Sexual_Partners_bin,
                                          levels = c("3+", "1 or 2", "0", "Unknown"))
your_data_all$Lifetime_Sexual_Partners_bin <- revalue(your_data_all$Lifetime_Sexual_Partners_bin, 
                                   c("0" = "0",
                                     "1 or 2" = "1 or 2",
                                     "3+" = "3+",
                                     "Unknown" = "Not disclosed"))
your_data_all$Sexually_active_6months <- factor(your_data_all$Sexually_active_6months,
                                                     levels = c("Yes", "No", "Not disclosed"))
your_data_all$Urogen_Sympt_Cat <- factor(your_data_all$Urogen_Sympt_Cat,
                       levels = c("2", "1", "0"))
your_data_all$Urogen_Sympt_Cat <- revalue(your_data_all$Urogen_Sympt_Cat, 
                                                      c("0" = "Not disclosed",
                                                        "1" = "No",
                                                        "2" = "Yes"))
your_data_all$Abx <- factor(your_data_all$Abx,
                                         levels = c("Yes", "No", "Unknown"))
your_data_all$Abx <- revalue(your_data_all$Abx, 
                                          c("No" = "No",
                                            "Yes" = "Yes",
                                            "Unknown" = "Not disclosed"))
your_data_all$Diabetes_status <- factor(your_data_all$Diabetes_status,
                            levels = c("T2D", "T1D", "Prediabetic", "None"))


##### Table 1 and packages/libraries to import #####
VSPall_triala <-
  your_data_all %>%
  dplyr::select(RA_status, Age, Clinic, Ethnicity, Menopausal_bin, Diet_bin, Alcohol_consumption, Cigarette_smoking, Education, Pregnancies, Lifetime_Sexual_Partners_bin, Sexually_active_6months, Diabetes_status, Urogen_Sympt_Cat, Abx, CST_assignment)
head("VSP cases and controls")
tbl_summary_1a <-  VSPall_triala %>%
  tbl_summary(
    by = RA_status,
    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} / {N} ({p}%)"),
    label = list(Alcohol_consumption ~ "Alcohol Consumption", Diet_bin ~ "Diet (Meat/Fat intake)", Cigarette_smoking ~ "Cigarette Smoking", Menopausal_bin ~ "Menopausal Status", Sexually_active_6months ~ "Sexually Active", Pregnancies ~ "Lifetime Pregnancies", Lifetime_Sexual_Partners_bin ~ "Lifeitime Sexual partners", Diabetes_status ~ "Diabetes Status", Urogen_Sympt_Cat ~ "Urogential Symptoms", Abx ~ "Antibiotic Use", CST_assignment ~ "CST"),
    digits = all_continuous() ~ 0) %>%
  add_overall () %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  modify_footnote(all_stat_cols() ~ "Median (range) or Frequency (%)") %>%
  modify_caption("**Table 1. Clinical demographics of study participants stratified by RA status**") %>%
  bold_labels()
tbl_summary_1a
#save html, then print to PDF, then export image



##### Dataframes of multiple conditions ####
your_data_cyt <- read.csv("Metadata_RAstudy_master_categorized_v2lessIMMUNE.csv") #partial USE FOR AUC Immune ONLY COMP
your_data_microbes <- read.csv("Metadata_RAstudy_master_categorized_v2less.csv") #partial USE FOR AUC MICROBES ONLY COMP
your_data_immune <- read.csv("Metadata_RAstudy_master_categorized_cytokines.csv") #partial USE FOR AUC Immune ONLY COMP
your_data_immune_abx <- read.csv("Metadata_RAstudy_master_categorized_cytokines_lessabx.csv") #makes no diff in ABX
df_RF <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagRF.csv") #ELISA
df_RFlsp <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagRFlsp.csv")
df_RFlspabx <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagRFlspabx.csv") #ELISA #makes no diff in ABX
df_ACCP <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagACCP.csv") #ELISA
df_ACCPlsp <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagACCPlsp.csv") #ELISA
df_ACCPlspabx <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagACCPlspabx.csv") #makes no diff in ABX
df_CRP <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRP.csv") #ELISA
df_CRPlsp <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRPlsp.csv") #ELISA
df_CRPlspabx <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRPlspabx.csv") #ELISA makes no diff in ABX
df_CRPserum <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRPserum.csv") #ELISA serum vs vag
df_CRPserum_ACPA <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRPserum_xACPA.csv") #ELISA serum vs vag

df_meds <- read.csv("Metadata_RAstudy_master_RAcategorized_meds_yn.csv") #for all except JSN
your_data_RA <- read.csv("Metadata_RAstudy_master_RAcategorized.csv")
your_data_RA_vagELISAS <- read.csv("Metadata_RAstudy_master_RAcategorized_vagELISAS.csv")
your_data_RA_vagELISASxACPA <- read.csv("Metadata_RAstudy_master_RAcategorized_vagELISASxACPA.csv")
your_data_RA_ACPA <- read.csv("Metadata_RAstudy_master_RAcategorized_xACPA.csv") #RA only _ full (CDAI)
your_data_RA_ACPAlspabx <- read.csv("Metadata_RAstudy_master_RAcategorized_xACPAlspabx.csv") #no diff without
your_data_4factors <- read.csv("Metadata_RAstudy_master_RAcategorized_4factors.csv")
your_data_medrank <- read.csv("Metadata_RAstudy_master_RAcategorized_meds.csv") #lm 
your_data_4factorslspabx <- read.csv("Metadata_RAstudy_master_RAcategorized_4factors_lspabx.csv") #JSN, RE, ACCP, RF

#### "Table 1" generation ####
VSPall_group <-
  your_data_microbes_fin %>%
  select(RA_status, Age, Clinic, Ethnicity, Menopausal_bin, Diet_bin, Alcohol_consumption, Cigarette_smoking, Education, Pregnancies, Lifetime_Sexual_Partners_bin, Sexually_active_6months, Urogen_Sympt_Cat, Abx, CST_assignment)
tbl_summary_1b <-  VSPall_group %>%
  tbl_summary(
    by = RA_status,
    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} / {N} ({p}%)"),
    label = list(Alcohol_consumption ~ "Alcohol Consumption", Diet_bin ~ "Diet (Meat/Fat intake)", Cigarette_smoking ~ "Cigarette Smoking", Menopausal_bin ~ "Menopausal Status", Sexually_active_6months ~ "Sexually Active", Pregnancies ~ "Lifetime Pregnancies", Lifetime_Sexual_Partners_bin ~ "Lifeitime Sexual partners", Urogen_Sympt_Cat ~ "Urogential Symptoms", Abx ~ "Antibiotic Use", CST_assignment ~ "CST"),
    digits = all_continuous() ~ 0) %>%
  add_overall () %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  modify_footnote(all_stat_cols() ~ "Median (range) or Frequency (%)") %>%
  modify_caption("**Table 1. Clinical Case-Control Demographics**") %>%
  bold_labels()
tbl_summary_1b

VSPall_RA <-
  your_data_RA_vagELISAS %>%
  select(AntiCCP_Ab_positive, Age, Clinic, Disease_Activity, Radiographic_Erosions, Joint_space_narrowing, RF_positive, ACPA_detectable, RF_detectable, C_reactive_detectable, Ethnicity, Menopausal_bin, Diet_bin, Alcohol_consumption, Cigarette_smoking, Education, Pregnancies, Lifetime_Sexual_Partners_bin, Sexually_active_6months, Urogen_Sympt_Cat, Abx, CST_assignment)
tbl_summary_1c <-  VSPall_RA %>%
  tbl_summary(
    by = AntiCCP_Ab_positive,
    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} / {N} ({p}%)"),
    label = list(Alcohol_consumption ~ "Alcohol Consumption", Diet_bin ~ "Diet (Meat/Fat intake)", Cigarette_smoking ~ "Cigarette Smoking", Menopausal_bin ~ "Menopausal Status", Sexually_active_6months ~ "Sexually Active", Pregnancies ~ "Lifetime Pregnancies", Lifetime_Sexual_Partners_bin ~ "Lifeitime Sexual partners", Urogen_Sympt_Cat ~ "Urogential Symptoms", Abx ~ "Antibiotic Use", CST_assignment ~ "CST"),
    digits = all_continuous() ~ 0) %>%
  add_overall () %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  modify_footnote(all_stat_cols() ~ "Median (range) or Frequency (%)") %>%
  modify_caption("**Table 1. Clinical Case-Control Demographics**") %>%
  bold_labels()
tbl_summary_1c

VSPall_RA_acpa <-
  your_data_RA_vagELISASxACPA %>%
  select(AntiCCP_Ab_positive, Age, Clinic, Disease_Activity, Radiographic_Erosions, Joint_space_narrowing, RF_positive, ACPA_detectable, RF_detectable, C_reactive_detectable, Ethnicity, Menopausal_bin, Diet_bin, Alcohol_consumption, Cigarette_smoking, Education, Pregnancies, Lifetime_Sexual_Partners_bin, Sexually_active_6months, Urogen_Sympt_Cat, Abx, CST_assignment)
tbl_summary_1d <-  VSPall_RA_acpa %>%
  tbl_summary(
    by = Radiographic_Erosions,
    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} / {N} ({p}%)"),
    label = list(Alcohol_consumption ~ "Alcohol Consumption", Diet_bin ~ "Diet (Meat/Fat intake)", Cigarette_smoking ~ "Cigarette Smoking", Menopausal_bin ~ "Menopausal Status", Sexually_active_6months ~ "Sexually Active", Pregnancies ~ "Lifetime Pregnancies", Lifetime_Sexual_Partners_bin ~ "Lifeitime Sexual partners", Urogen_Sympt_Cat ~ "Urogential Symptoms", Abx ~ "Antibiotic Use", CST_assignment ~ "CST"),
    digits = all_continuous() ~ 0) %>%
  add_overall () %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  modify_footnote(all_stat_cols() ~ "Median (range) or Frequency (%)") %>%
  modify_caption("**Table 1. Clinical Case-Control Demographics**") %>%
  bold_labels()
tbl_summary_1d


####FOR USE OFFICIALLY####
your_data_microbes_fin <- read.csv("Metadata_RAstudy_master_categorized_v2excess_less.csv") #partial set USE FOR AUC Immune ONLY COMP
your_data_immune_fin <- read.csv("Metadata_RAstudy_master_categorized_cytokines_less.csv") #partial set USE FOR AUC Immune ONLY COMP
df_ACCPlsp <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagACCPlsp.csv") #ELISA
df_RFlsp <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagRFlsp.csv") #ELISA #converted all 0 to 0.001
df_CRPlsp <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRPlsp.csv") #ELISA

your_data_4factors <- read.csv("Metadata_RAstudy_master_RAcategorized_4factors.csv") #JSN, RE, ACCP, RF
your_data_4factorsESR <- read.csv("Metadata_RAstudy_master_RAcategorized_4factors_ESR.csv")
your_data_RA_vagELISASxACPA <- read.csv("Metadata_RAstudy_master_RAcategorized_vagELISASxACPA.csv")
df_CRPserum_ACPA <- read.csv("Metadata_RAstudy_master_categorized_v2excess_vagCRPserum_xACPA.csv") #ELISA serum vs vag

cytokine_data <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few.csv") #for visualizations of cytokine data (log)
cytokine_data_MEDS <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few_RAonlyMEDS.csv") #for visualizations of cytokine data (log)
df_meds <- read.csv("Metadata_RAstudy_master_RAcategorized_meds_yn.csv") #yes and no (current, not - same outputs as using 0 and 1)

#include Menopausal_bin_ordered + Diet_bin in all models
# Find the optimal threshold (e.g., maximizing Youden Index)
optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
print(optimal_threshold)
zv


##### MICROBES #####

#Print p-val and coefficients of microbial impact on RA status
head(your_data_microbes_fin[212:251])
microbes <- colnames(your_data_microbes_fin)[212:251]
factors <- c("RA_status_cat")
p_value_micmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                         dimnames = list(microbes, factors))
coef_micmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                      dimnames = list(microbes, factors))
# Loop through microbes and factors
for (microbe in microbes) {
  for (factor in factors) {
    # Build the formula
    formula_mic <- as.formula(paste(microbe, "~ Diet_bin + Clinic_cat + Menopausal_bin_cat + Abx_Cat + Alohocl_cat +", factor))
    # Fit the model
    model_mic <- tryCatch(lm(formula_mic, data = your_data_microbes_fin), error = function(e) NULL)
    coef_summary <- summary(model_mic)$coefficients
    if (factor %in% rownames(coef_summary)) {
      p_value_micmatrix[microbe, factor] <- coef_summary[factor, "Pr(>|t|)"]
      coef_micmatrix[microbe, factor] <- coef_summary[factor, "Estimate"]
    }}}
# View the resulting p-value matrix
print(p_value_micmatrix)
print(coef_micmatrix)
write.csv(p_value_micmatrix, "Summary/p_value_matrix_RAefMicrobe.csv", row.names = TRUE)
write.csv(coef_micmatrix, "Summary/coef_matrixRAefMicrobe.csv", row.names = TRUE)

for (microbe in microbes) {
  for (factor in factors) {
    # Build the formula
    formula_mic <- as.formula(paste(microbe, "~ Diet_bin + Clinic_cat + Menopausal_bin_cat + Abx_Cat +", factor))
    # Fit the model
    model_mic <- tryCatch(lm(formula_mic, data = your_data_microbes_fin), error = function(e) NULL)
    coef_summary <- summary(model_mic)$coefficients
    if (factor %in% rownames(coef_summary)) {
      p_value_micmatrix[microbe, factor] <- coef_summary[factor, "Pr(>|t|)"]
      coef_micmatrix[microbe, factor] <- coef_summary[factor, "Estimate"]
  }}}
# View the resulting p-value matrix
print(p_value_micmatrix)
print(coef_micmatrix)
write.csv(p_value_micmatrix, "Summary/p_value_matrix_RAefMicrobe_noalcohol.csv", row.names = TRUE)
write.csv(coef_micmatrix, "Summary/coef_matrixRAefMicrobe_noalcohol.csv", row.names = TRUE)

#Print p-val and coefficients of microbial impact on RA status
p_value_micmatrix2 <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                             dimnames = list(factors, microbes))
coef_micmatrix2 <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                          dimnames = list(factors, microbes))
# Loop through factors and microbes
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(factor, "~ Diet_bin + Clinic_cat + Menopausal_bin_cat + Abx_Cat +", microbe))
    # Fit the model with error handling
    model_mic2 <- tryCatch(glm(formula_mic2, family = binomial(link = "logit"), data = your_data_microbes_fin), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary2 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        p_value_micmatrix2[factor, microbe] <- coef_summary2[microbe, "Pr(>|z|)"]
        coef_micmatrix2[factor, microbe] <- coef_summary2[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_micmatrix2)
print("Coefficient Matrix:")
print(coef_micmatrix2)
# Optionally save to CSV files
write.csv(p_value_micmatrix2, "Summary/p_value_micmatrix2_MICROBEefRA_glm.csv", row.names = TRUE)
write.csv(coef_micmatrix2, "Summary/coef_micmatrix2_MICROBEefRA_glm.csv", row.names = TRUE)


#attempting to generate a good model - resorted to stating the limitations in paper. sig but not linear
nl_formula <- function(beta0, beta1, beta2, beta3, RA_status_cat, Diet_bin_cat, Menopausal_bin_cat) {
  beta0 + beta1 * RA_status_cat + beta2 * Diet_bin_cat + beta3 * Menopausal_bin_cat + I(RA_status_cat^2)
}
# Model fitting
nlme_model <- nlme(
  model = Peptoniphilus_lacrimalis ~ nl_formula(beta0, beta1, beta2, beta3, RA_status_cat, Diet_bin_cat, Menopausal_bin_cat),
  data = your_data_microbes_fin,
  fixed = beta0 + beta1 + beta2 + beta3 ~ 1,
  random = beta0 ~ 1 | Abx_Cat/Alcohol_cat/Clinic_cat,
  start = c(beta0 = 1, beta1 = 0.00353, beta2 = -0.00244, beta3 = -0.000629),
  method = "ML"
)
# Summary of the model
summary(nlme_model)
vif(nlme_model)
par(mfrow = c(2, 2))
plot(nlme_model)
plot(nlme_model, which = 1])
residuals <- resid(nlme_model)
fitted_values <- fitted(nlme_model)


##### MICROBES model #####
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Total_Bifido + Total_Peptoniphilus + Lactobacillus_crispatus +
             Lactobacillus_gasseri +	Lactobacillus_jensenii +	Megasphaera_lornae +
               UBA629_sp005465875 +	Bifidobacterium_leopoldii	+ Bifidobacterium_breve	+ Bifidobacterium_	+
               Prevotella_timonensis	+	Prevotella_bivia	+ Sneathia_vaginalis	+ Peptoniphilus_lacrimalis	+
               Campylobacter_ureolyticus	+ Dialister_sp001553355	+ Dialister_propionicifaciens	+
               Dialister_micraerophilus +	Anaerococcus_vaginalis + observed_features + shannon_entropy +
               Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
             family = binomial(link = "logit"), data = your_data_microbes_fin)
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Total_Bifido +
              Megasphaera_lornae + Bifidobacterium_leopoldii	+
               Peptoniphilus_lacrimalis	+ Dialister_sp001553355 +
               Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
             family = binomial(link = "logit"), data = your_data_microbes_fin)
#AIC: 88.624 #AUC: 0.8925
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Total_Bifido +
              Megasphaera_lornae + Bifidobacterium_leopoldii	+
               Dialister_sp001553355 +
               Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
             family = binomial(link = "logit"), data = your_data_microbes_fin)
#AIC: 91.2 #AUC: 0.8785
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Total_Bifido + Total_Peptoniphilus +
              Megasphaera_lornae +
               Dialister_sp001553355	+
              Anaerococcus_vaginalis + 
               Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
             family = binomial(link = "logit"), data = your_data_microbes_fin)#AIC: 94.sbr59
#AIC: 86.447 #AUC: 0.9013
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Total_Bifido + Total_Peptoniphilus +
              Megasphaera_lornae +
              Bifidobacterium_leopoldii	+
              Dialister_sp001553355	+
              Anaerococcus_vaginalis + 
               Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
             family = binomial(link = "logit"), data = your_data_microbes_fin)#AIC: 94.sbr59
#AIC: 83.2416 #AUC: 0.9139 USE THIS
model_microbes <- glm(RA_status_cat ~ Total_Prevotella +
                        Megasphaera_lornae +
                        Dialister_sp001553355	+
                        Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
                      family = binomial(link = "logit"), data = your_data_microbes_fin)
#AIC: 94.045 #AUC: 0.8554
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Campylobacter_ureolyticus +
                        Megasphaera_lornae +
                        Dialister_sp001553355	+
                        Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
                      family = binomial(link = "logit"), data = your_data_microbes_fin)
#AIC: 91.968 #AUC: 0.8631


##### Microbes final model#####
model_microbes <- glm(RA_status_cat ~ Total_Prevotella + Campylobacter_ureolyticus +
                        Megasphaera_lornae +
                        Dialister_sp001553355	+
                        Menopausal_bin_ordered	+	Clinic_cat + Diet_bin + Abx_Cat + Alcohol_cat,
                      family = binomial(link = "logit"), data = your_data_microbes_fin)
#AIC: 91.968 #AUC: 0.8631


summary(model_microbes)
vif(model_microbes)
pred_probs <- predict(model_microbes, type = "response")
length(your_data_microbes_fin$RA_status)
# Convert probabilities to binary outcome (threshold = 0.5) %0.3 and 0.4878675 yielded same AUC, worse confusion matrices
pred_class <- ifelse(pred_probs > 0.5, 1, 0)
length(pred_class)
# Confusion matrix
table(Predicted = pred_class, Actual = your_data_microbes_fin$RA_status)
roc_curve <- roc(your_data_microbes_fin$RA_status_cat, pred_probs)
plot(roc_curve)
auc(roc_curve)




summary(model_microbes)
vif(model_microbes)
AIC(model_microbes)
pred_probsmic <- predict(model_microbes, type = "response")
pred_classmic <- ifelse(pred_probsmic > 0.5, 1, 0)
table(Predicted = pred_classmic, Actual = your_data_microbes_fin$RA_status_cat)
roc_curvemic <- roc(your_data_microbes_fin$RA_status_cat, pred_probsmic)
plot(roc_curvemic)
auc(roc_curvemic)








#### RA ONLY MICROBES ####
#Print p-val and coefficients of microbial impact on RA status
head(your_data_4factors[211:285])
microbes <- colnames(your_data_4factors)[211:286] #*include to 286 for lm)
factors1 <- c("CDAI_score", "Disease_Activity_cat", "AntiCCP_cat_Ab_positive",	"RF_cat_positive")
factors2 <- c("Radiographic_cat_Erosions", "Joint_cat_space_narrowing")
p_value_micRAmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors1),
                            dimnames = list(microbes, factors1))
coef_micRAmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors1),
                         dimnames = list(microbes, factors1))

#TESTING
best_model <- stepAIC(lm(Lactobacillus_gasseri ~ RF_cat_positive * (Menopausal_bin_cat + Diet_bin), data = your_data_4factors))
summary(best_model)
best_model <- stepAIC(lm(Peptoniphilus_harei ~ RF_cat_positive * (Menopausal_bin_cat + Diet_bin), data = your_data_4factors))
summary(best_model)
next_model <- lm(Peptoniphilus_harei ~ RF_cat_positive + Menopausal_bin_cat + Diet_bin, data = your_data_4factors)
summary(next_model)
best_model <- stepAIC(lm(Prevotella_colorans ~ RF_cat_positive * (Menopausal_bin_cat + Diet_bin), data = your_data_4factors))
summary(best_model)
best_model <- stepAIC(lm(CDAI_score ~ Bifidobacterium_vaginale * (Menopausal_bin_cat + Diet_bin), data = your_data_4factors))
summary(best_model)
best_model <- stepAIC(glm(CDAI_score ~ Bifidobacterium_vaginale * (Menopausal_bin_cat + Diet_bin), family=poisson(link="log"), data = your_data_4factors))
summary(best_model)
next_model <- glm(CDAI_score ~ Bifidobacterium_vaginale + Menopausal_bin_cat + Diet_bin, family=poisson(link="log"), data = your_data_4factors)
summary(next_model)
best_model <- stepAIC(glm(CDAI_score ~ Peptoniphilus_lacrimalis * (Menopausal_bin_cat + Diet_bin), family=poisson(link="log"), data = your_data_4factors))
summary(best_model)
next_model <- glm(CDAI_score ~ Peptoniphilus_lacrimalis + Menopausal_bin_cat + Diet_bin, family=poisson(link="log"), data = your_data_4factors)
summary(next_model)
best_model <- stepAIC(glm(CDAI_score ~ Fenollaria_massiliensis * (Menopausal_bin_cat + Diet_bin), family=poisson(link="log"), data = your_data_4factors))
summary(best_model)
next_model <- glm(CDAI_score ~ Fenollaria_massiliensis + Menopausal_bin_cat + Diet_bin, family=poisson(link="log"), data = your_data_4factors)
summary(next_model)
"
Lactobacillus_iners
Megasphaera_lornae
Sneathia_vaginalis
Anaerococcus_tetradius
Dialister_sp001553355
Prevotella_bivia
Peptoniphilus_lacrimalis
Bifidobacterium_vaginale
Fenollaria_massiliensis"

#END CONFIRMATION TESTING



# Loop through microbes and factors
for (microbe in microbes) {
  for (factor in factors1) {
    # Build the formula
    formula_micRAa <- as.formula(paste(microbe, "~ Diet_bin + Menopausal_bin_cat +", factor))
    # Fit the model
    model_micRA <- tryCatch(lm(formula_micRAa, data = your_data_4factors), error = function(e) NULL)
    coef_summary <- summary(model_micRA)$coefficients
    if (factor %in% rownames(coef_summary)) {
      p_value_micRAmatrix[microbe, factor] <- coef_summary[factor, "Pr(>|t|)"]
      coef_micRAmatrix[microbe, factor] <- coef_summary[factor, "Estimate"]
    }}}
# View the resulting p-value matrix
print(p_value_micRAmatrix)
print(coef_micRAmatrix)
write.csv(p_value_micRAmatrix, "Summary/p_value_matrix_CDAIACCPRF_efMicrobe.csv", row.names = TRUE)
write.csv(coef_micRAmatrix, "Summary/coef_matrix_CDAIACCPRF_efMicrobe.csv", row.names = TRUE)

p_value_micRAmatrix2 <- matrix(NA, nrow = length(microbes), ncol = length(factors2),
                              dimnames = list(microbes, factors2))
coef_micRAmatrix2 <- matrix(NA, nrow = length(microbes), ncol = length(factors2),
                           dimnames = list(microbes, factors2))
for (microbe in microbes) {
  for (factor in factors2) {
    # Build the formula
    formula_micRAb <- as.formula(paste(microbe, "~ Diet_bin + Menopausal_bin_cat + AntiCCP_cat_Ab_positive +", factor))
    # Fit the model
    model_micRA <- tryCatch(lm(formula_micRAb, data = your_data_4factors), error = function(e) NULL)
    coef_summary <- summary(model_micRA)$coefficients
    if (factor %in% rownames(coef_summary)) {
      p_value_micRAmatrix2[microbe, factor] <- coef_summary[factor, "Pr(>|t|)"]
      coef_micRAmatrix2[microbe, factor] <- coef_summary[factor, "Estimate"]
    }}}
# View the resulting p-value matrix
print(p_value_micRAmatrix2)
print(coef_micRAmatrix2)
write.csv(p_value_micRAmatrix2, "Summary/p_value_matrix_JSNRE_efMicrobe.csv", row.names = TRUE)
write.csv(coef_micRAmatrix2, "Summary/coef_matrix_JSNRE_efMicrobe.csv", row.names = TRUE)



#TESTING IT GOES THROUGH
if (!is.null(model_micRA2)) {
  coef_summary2 <- summary(model_micRA2)$coefficients
  
  # Ensure coef_summary2 is a matrix or data frame
  if (is.matrix(coef_summary2) || is.data.frame(coef_summary2)) {
    # Check if 'microbe' exists in the row names of coef_summary2
    if (any(grepl(microbe, rownames(coef_summary2)))) {
      # Extract the coefficient related to the microbe
      coef_summary2_filtered <- coef_summary2[grep(microbe, rownames(coef_summary2)), , drop = FALSE]  # drop=FALSE ensures it's still a data frame or matrix
      
      # Now access the p-value and estimate correctly
      p_value_micRAmatrix3[factor, microbe] <- coef_summary2_filtered[1, "Pr(>|z|)"]
      coef_micRAmatrix3[factor, microbe] <- coef_summary2_filtered[1, "Estimate"]
    } else {
      message(paste("Microbe", microbe, "not found in the model for factor", factor))
    }
  } else {
    message("Model summary coefficients are not in the expected format (matrix/data.frame).")
  }
} else {
  message(paste("Model failed for factor", factor, "and microbe", microbe))
}
install.packages("logistf")
library(logistf)
model_micRA2_firth <- tryCatch(
  logistf(formula_micRA2a, data = your_data_4factors),
  error = function(e) NULL
)
cor(your_data_4factors[, c("RF_cat_positive", "Diet_bin", "Menopausal_bin_cat", "shannon_entropy")], use = "complete.obs")

#Print p-val and coefficients of microbial impact on RA status
# Loop through factors and microbes
microbes_lim <- colnames(your_data_4factors)[211:285] #*include to 286 for lm)
model_microbes <- glm(CDAI_score ~ Diet_bin + Menopausal_bin_cat + Total_Lactobacillus,
                      family= poisson(link="log"), data = your_data_4factors)
summary(model_microbes)

microbes_lim <- colnames(your_data_4factors)[211:285] #*include to 286 for lm)
factors1ax <- c("CDAI_score") #Disease_Activity_cat
coef_micRAmatrix3a <- matrix(NA, nrow = length(factors1ax), ncol = length(microbes_lim),
                             dimnames = list(factors1ax, microbes_lim))
stderror_micRAmatrix3a <- matrix(NA, nrow = length(factors1ax), ncol = length(microbes_lim),
                                 dimnames = list(factors1ax, microbes_lim))
zvalue_micRAmatrix3a <- matrix(NA, nrow = length(factors1ax), ncol = length(microbes_lim),
                               dimnames = list(factors1ax, microbes_lim))
p_value_micRAmatrix3a <- matrix(NA, nrow = length(factors1ax), ncol = length(microbes_lim),
                                dimnames = list(factors1ax, microbes_lim))

for (factor in factors1ax) {
  for (microbe in microbes_lim) {
    # Build the formula
    formula_micRA2a <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", microbe))
    # Fit the model with error handling
    model_micRA2 <- tryCatch(glm(formula_micRA2a, family = poisson(link="log"), data = your_data_4factors), 
                             error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_micRA2)) {
      coef_summary2 <- summary(model_micRA2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        coef_micRAmatrix3a[factor, microbe] <- coef_summary2[microbe, "Estimate"]
        stderror_micRAmatrix3a[factor, microbe] <- coef_summary2[microbe, "Std. Error"]
        zvalue_micRAmatrix3a[factor, microbe] <- coef_summary2[microbe, "z value"]
        p_value_micRAmatrix3a[factor, microbe] <- coef_summary2[microbe, "Pr(>|z|)"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_micRAmatrix3a)
print("stderror Matrix:")
print(stderror_micRAmatrix3a)
print("z-value Matrix:")
print(zvalue_micRAmatrix3a)
print("Coefficient Matrix:")
print(coef_micRAmatrix3a)
# Optionally save to CSV files
write.csv(p_value_micRAmatrix3a, "Summary/p_value_micRAmatrix2_micRAROBEef_CDAI_redo.csv", row.names = TRUE)
write.csv(stderror_micRAmatrix3a, "Summary/stderror_micRAmatrix2_micRAROBEef_CDAI_redo.csv", row.names = TRUE)
write.csv(zvalue_micRAmatrix3a, "Summary/zvalue_micRAmatrix2_micRAROBEef_CDAI_redo.csv", row.names = TRUE)
write.csv(coef_micRAmatrix3a, "Summary/coef_micRAmatrix2_micRAROBEef_CDAI_redo.csv", row.names = TRUE)


p_value_micRAmatrix3alm <- matrix(NA, nrow = length(factors1ax), ncol = length(microbes_lim),
                                dimnames = list(factors1ax, microbes_lim))
coef_micRAmatrix3alm <- matrix(NA, nrow = length(factors1ax), ncol = length(microbes_lim),
                             dimnames = list(factors1ax, microbes_lim))
for (factor in factors1ax) {
  for (microbe in microbes_lim) {
    # Build the formula
    formula_micRA2a <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", microbe))
    # Fit the model with error handling
    model_micRA2 <- tryCatch(lm(formula_micRA2a, data = your_data_4factors), error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_micRA2)) {
      coef_summary2 <- summary(model_micRA2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        p_value_micRAmatrix3alm[factor, microbe] <- coef_summary2[microbe, "Pr(>|t|)"]
        coef_micRAmatrix3alm[factor, microbe] <- coef_summary2[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
write.csv(p_value_micRAmatrix3alm, "Summary/p_value_micRAmatrix2_micRAROBEef_CDAI_lm.csv", row.names = TRUE)
write.csv(coef_micRAmatrix3alm, "Summary/coef_micRAmatrix2_micRAROBEef_CDAI_lm.csv", row.names = TRUE)

for (factor in factors1ax) {
  for (microbe in microbes_lim) {
    # Build the formula
    formula_micRA2a <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", microbe))
    # Fit the model with error handling
    model_micRA2 <- tryCatch(glm(formula_micRA2a, family = poisson(link="log"), data = your_data_4factors), 
                             error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_micRA2)) {
      coef_summary2 <- summary(model_micRA2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        p_value_micRAmatrix3a[factor, microbe] <- coef_summary2[microbe, "Pr(>|z|)"]
        coef_micRAmatrix3a[factor, microbe] <- coef_summary2[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
print("P-Value Matrix:")
print(p_value_micRAmatrix3a)
print("Coefficient Matrix:")
print(coef_micRAmatrix3a)
write.csv(p_value_micRAmatrix3a, "Summary/p_value_micRAmatrix2_micRAROBEef_CDAI_redo.csv", row.names = TRUE)
write.csv(coef_micRAmatrix3a, "Summary/coef_micRAmatrix2_micRAROBEef_CDAI_redo.csv", row.names = TRUE)


factors1b <- c("AntiCCP_cat_Ab_positive",	"RF_cat_positive")
for (factor in factors1b) {
  for (microbe in microbes_lim) {
    # Build the formula
    formula_micRA2a <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", microbe))
    # Fit the model with error handling
    model_micRA2 <- tryCatch(glm(formula_micRA2a, family = binomial(link="logit"), data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_micRA2)) {
      coef_summary2 <- summary(model_micRA2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        p_value_micRAmatrix3b[factor, microbe] <- coef_summary2[microbe, "Pr(>|z|)"]
        coef_micRAmatrix3b[factor, microbe] <- coef_summary2[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
write.csv(p_value_micRAmatrix3b, "Summary/p_value_micRAmatrix2_micRAROBEef_ACCPRF.csv", row.names = TRUE)
write.csv(coef_micRAmatrix3b, "Summary/coef_micRAmatrix2_micRAROBEef_ACCPRF.csv", row.names = TRUE)


p_value_micRAmatrix4 <- matrix(NA, nrow = length(factors2), ncol = length(microbes),
                               dimnames = list(factors2, microbes))
coef_micRAmatrix4 <- matrix(NA, nrow = length(factors2), ncol = length(microbes),
                            dimnames = list(factors2, microbes))
for (factor in factors2) {
  for (microbe in microbes) {
    # Build the formula
    formula_micRA2b <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat + AntiCCP_cat_Ab_positive +", microbe))
    # Fit the model with error handling
    model_micRA2 <- tryCatch(glm(formula_micRA2b, family = binomial(link = "logit"), data = your_data_4factors), 
                             error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_micRA2)) {
      coef_summary2 <- summary(model_micRA2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        p_value_micRAmatrix4[factor, microbe] <- coef_summary2[microbe, "Pr(>|z|)"]
        coef_micRAmatrix4[factor, microbe] <- coef_summary2[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
write.csv(p_value_micRAmatrix4, "Summary/p_value_micRAmatrix2_micRAROBEef_JSNRF.csv", row.names = TRUE)
write.csv(coef_micRAmatrix4, "Summary/coef_micRAmatrix2_micRAROBEef_JSNRF.csv", row.names = TRUE)






#### RA ONLY MICROBES MODELS - Meds ####
par(mfrow = c(2, 2))
summary(model_microbesACCP)
vif(model_microbesACCP)
pred_probs <- predict(model_microbesACCP, type = "response")
length(your_data_microbes_fin$RA_status)
# Convert probabilities to binary outcome (threshold = 0.5) %0.3 and 0.4878675 yielded same AUC, worse confusion matrices
pred_class <- ifelse(pred_probs > 0.5, 1, 0)
length(pred_class)
# Confusion matrix
table(Predicted = pred_class, Actual = model_microbesACCP$RA_status)
roc_curve <- roc(model_microbesACCP$RA_status_cat, pred_probs)
plot(roc_curve)
auc(roc_curve)

head(your_data_4factors[210:286])
microbes <- colnames(your_data_4factors)[211:286] #*include to 286 for lm)
medications <- c("Hydroxychloroquine", "Sulfasalazine", "Methotrexate_dose_.mg.", "Prednisone_dose_.mg.", "Methotrexate", "Prednisone")
p_value_matrix <- matrix(NA, nrow = length(microbes), ncol = length(medications),
                         dimnames = list(microbes, medications))
coef_matrix <- matrix(NA, nrow = length(microbes), ncol = length(medications),
                      dimnames = list(microbes, medications))

# Loop through microbes and medications (not great because so few of bacteria present and no linear)
for (microbe in microbes) {
  for (med in medications) {
    # Build the formula
    formula <- as.formula(paste(microbe, "~ Diet_bin + Clinic_cat + Menopausal_bin_ordered +", med))
    
    # Fit the model
    model <- lm(formula, data = df_meds)
    
    # Extract the p-value for the medication variable
    p_value <- summary(model)$coefficients[med, "Pr(>|t|)"]
    estimate_coef <- summary(model)$coefficients[med, "Estimate"]
    
    # Store the p-value in the matrix
    p_value_matrix[microbes, med] <- p_value
    coef_matrix[microbes, med] <- estimate_coef
    
  }
}

# View the resulting p-value matrix
print(p_value_matrix)

# Optionally save to a CSV file
write.csv(p_value_matrix, "p_value_matrix_medsefmicrobes_yn.csv", row.names = TRUE)
write.csv(coef_matrix, "coef_matrix_medsefmicrobes_yn.csv", row.names = TRUE)


model_Sulf <- lm(shannon_entropy ~ Sulfasalazine, 
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ Sulfasalazine, 
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ Sulfasalazine + shannon_entropy,
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ shannon_entropy,
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ Sulfasalazine * shannon_entropy,
                 data = df_meds)
summary(model_Sulf)
vif(model_Sulf)
par(mfrow = c(2, 2))
plot(model_Sulf)

cor.test(df_meds$Methotrexate_dose_.mg., df_meds$Total_Lactobacillus)

model_methDose <- lm(Total_Lactobacillus ~ Methotrexate_dose_.mg.,
                     data = df_meds)
model_methDose <- lm(Total_Lactobacillus ~ Methotrexate_dose_.mg. +
                       Megasphaera_lornae +
                       Fannyhessea_vaginae +
                       Propionimicrobium_lymphophilum+
                       Bifidobacterium_vaginale +
                       Diet_cat + Clinic_cat + Menopausal_bin,
                     data = df_meds)
model_methDose <- lm(Methotrexate_dose_.mg. ~ Total_Lactobacillus +
                       Megasphaera_lornae +
                       Fannyhessea_vaginae +
                       Propionimicrobium_lymphophilum+
                       Bifidobacterium_vaginale +
                       Diet_cat + Clinic_cat + Menopausal_bin,
                     data = df_meds)
summary(model_methDose)
vif(model_methDose)
par(mfrow = c(2, 2))
plot(model_methDose)


model_predDose <- lm(Prednisone_dose_.mg. ~
                       Anaerococcus_tetradius +
                       Aerococcus_christensenii +
                       Diet_cat + Clinic_cat + Menopausal_bin,
                     data = df_meds)
model_predDose <- lm(Aerococcus_christensenii ~
                       Anaerococcus_tetradius *
                       Prednisone_dose_.mg. +
                       Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                     data = df_meds)
summary(model_predDose)
vif(model_predDose)
par(mfrow = c(2, 2))
plot(model_predDose)


model_Sulf <- lm(Lactobacillus_iners ~
                   Sulfasalazine +
                   Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                 data = df_meds)
#0.2721 (adjR 0.04719)
model_Sulf <- lm(Lactobacillus_iners ~
                   Sulfasalazine + shannon_entropy +
                   Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                 data = df_meds)
#0.01426 (adjR 0.2791) *USE THIS
model_Sulf <- lm(Lactobacillus_iners ~
                   Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                 data = df_meds)
#0.01022 (adjR 0.3192) *NO USE THIS!!! Sulf p=0.0566
model_Sulf <- lm(Lactobacillus_iners ~
                   shannon_entropy +
                   Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                 data = df_meds)
#0.0107 (adjR 0.2745)
model_Sulf <- lm(Lactobacillus_iners ~
                   Sulfasalazine * shannon_entropy +
                   Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                 data = df_meds)
#0.002515 (adjR 0.3959) *NO USE THIS!!! Sulf p=0.00852


model_Sulf <- lm(shannon_entropy ~ Sulfasalazine, 
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ Sulfasalazine, 
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ Sulfasalazine + shannon_entropy,
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ shannon_entropy,
                 data = df_meds)
model_Sulf <- lm(Lactobacillus_iners ~ Sulfasalazine * shannon_entropy,
                 data = df_meds)
summary(model_Sulf)
vif(model_Sulf)
par(mfrow = c(2, 2))
plot(model_Sulf)

summary(modelmed)
vif(modelmed)
AIC(modelmed)

# Predict probabilities
pred_probsmed <- predict(modelmed, type = "response")
length(df_meds$Sulfasalazine)
pred_classmed <- ifelse(pred_probsmed > 0.5, 1, 0)
length(pred_classmed)
table(Predicted = pred_classmed, Actual = df_meds$Sulfasalazine)
roc_curvemed <- roc(df_meds$Sulfasalazine, pred_probsmed)
plot(roc_curvemed)
auc(roc_curvemed)











##### IMMUNE #####
#Print p-val and coefficients of microbial impact on RA status
head(your_data_immune_fin[158:196])
cytokines <- colnames(your_data_immune_fin)[159:196]
factors <- c("RA_status_cat")
p_value_cytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                            dimnames = list(cytokines, factors))
coef_cytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                         dimnames = list(cytokines, factors))
# Loop through cytokines and factors
for (cytokine in cytokines) {
  for (factor in factors) {
    # Build the formula
    formula_cyt <- as.formula(paste(cytokine, "~ Diet_bin + Clinic_cat + Menopausal_bin_cat + Lifetime_Sexual_Partners_bin +", factor))
    # Fit the model
    model_cyt <- tryCatch(lm(formula_cyt, data = your_data_immune_fin), error = function(e) NULL)
    coef_summary <- summary(model_cyt)$coefficients
    if (factor %in% rownames(coef_summary)) {
      p_value_cytmatrix[cytokine, factor] <- coef_summary[factor, "Pr(>|t|)"]
      coef_cytmatrix[cytokine, factor] <- coef_summary[factor, "Estimate"]
    }}}
# View the resulting p-value matrix
print(p_value_cytmatrix)
print(coef_cytmatrix)
write.csv(p_value_cytmatrix, "Summary/p_value_cytmatrix_RAefcytokine.csv", row.names = TRUE)
write.csv(coef_cytmatrix, "Summary/coef_cytmatrix_RAefcytokine.csv", row.names = TRUE)

#Print p-val and coefficients of cytrobial impact on RA status
p_value_cytmatrix2 <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                             dimnames = list(factors, cytokines))
coef_cytmatrix2 <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                          dimnames = list(factors, cytokines))
p_value_cytmatrix2p <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                             dimnames = list(factors, cytokines))
coef_cytmatrix2p <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                          dimnames = list(factors, cytokines))
# Loop through factors and cytokines
#binomial
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Clinic_cat + Menopausal_bin_cat + Lifetime_Sexual_Partners_bin +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family = binomial(link ="logit"), data = your_data_immune_fin), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix2[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix2[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
#poisson
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2p <- as.formula(paste(factor, "~ Diet_bin + Clinic_cat + Menopausal_bin_cat + Lifetime_Sexual_Partners_bin +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2p, family = poisson(link ="log"), data = your_data_immune_fin), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix2p[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix2p[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_cytmatrix2)
print("Coefficient Matrix:")
print(coef_cytmatrix2)
# Optionally save to CSV files
write.csv(p_value_cytmatrix2p, "Summary/p_value_cytmatrix2_cytokineefRA_poisson.csv", row.names = TRUE)
write.csv(coef_cytmatrix2p, "Summary/coef_cytmatrix2_cytokineefRA_poisson.csv", row.names = TRUE)


#To see cytokines that correlate with menopausal status differently in RA and Control
factors <- c("Menopausal_bin_cat")
#poisson
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2p <- as.formula(paste(factor, "~ RA_status +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2p, family = poisson(link ="log"), data = your_data_immune_fin), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix2p[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix2p[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
# Print results
write.csv(p_value_cytmatrix2p, "Summary/p_value_cytmatrix2_menopauseefcytokineRA_poisson.csv", row.names = TRUE)
write.csv(coef_cytmatrix2p, "Summary/coef_cytmatrix2_menopauseefcytokineRA_poisson.csv", row.names = TRUE)

#just cytokine and menopausal state - no RA significance
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2p <- as.formula(paste(factor, "~", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2p, family = poisson(link ="log"), data = your_data_immune_fin), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix2p[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix2p[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
# Optionally save to CSV files
write.csv(p_value_cytmatrix2p, "Summary/p_value_cytmatrix2_menopauseefcytokine_poisson.csv", row.names = TRUE)
write.csv(coef_cytmatrix2p, "Summary/coef_cytmatrix2_menopauseefcytokine_poisson.csv", row.names = TRUE)

#TESTER TESTER
your_data_immune_fin$Menopausal_bin <- factor(your_data_immune_fin$Menopausal_bin, 
                                                     levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
your_data_immune_fin$RA_status <- factor(your_data_immune_fin$RA_status, 
                                              levels = c("Control", "Rheumatoid Arthritis"))
your_data_immune_fin$Diet_bin <- factor(your_data_immune_fin$Diet_bin, 
                                              levels = c("Reduced", "High"))
your_data_immune_fin$Clinic <- factor(your_data_immune_fin$Clinic, 
                                        levels = c("McNair", "Smith"))
your_data_immune_fin$Lifetime_Sexual_Partners_bin <- factor(your_data_immune_fin$Lifetime_Sexual_Partners_bin, 
                                      levels = c("0", "1 or 2", "3+"))

contrasts(your_data_immune_fin$RA_status)

#To see cytokines that correlate with factors differently in RA and Control
model_cyt1 <- lm(sCD40L ~
                  Lifetime_Sexual_Partners_bin + Menopausal_bin_cat + Clinic_cat + Diet_bin, data = your_data_immune_fin)
summary(model_cyt1)
model_cyt2 <- lm(sCD40L ~
                  RA_status + Diet_bin + Clinic_cat + Menopausal_bin_cat + Lifetime_Sexual_Partners_bin, data = your_data_immune_fin)
summary(model_cyt2)
model_cyt3 <- lm(sCD40L ~
                  RA_status:Diet_bin + Clinic_cat + Menopausal_bin_cat + Lifetime_Sexual_Partners_bin, data = your_data_immune_fin)
summary(model_cyt3)
model_cyt4 <- lm(sCD40L ~
                  RA_status:Menopausal_bin_cat + Clinic_cat + Diet_bin + Lifetime_Sexual_Partners_bin, data = your_data_immune_fin)
summary(model_cyt4)
anova(model_cyt1, model_cyt2, model_cyt3, model_cyt4)

install.packages("MASS")
library(MASS)
install.packages("emmeans")
library(emmeans)
your_data_immune_fin$Clinic <- relevel(your_data_immune_fin$Clinic, ref = "McNair")
your_data_immune_fin$Diet_bin <- relevel(your_data_immune_fin$Diet_bin, ref = "Reduced")
your_data_immune_fin$Menopausal_bin <- relevel(your_data_immune_fin$Menopausal_bin, ref = "Post-menopause")
your_data_immune_fin$Lifetime_Sexual_Partners_bin <- relevel(your_data_immune_fin$Lifetime_Sexual_Partners_bin, ref = "0")


best_model <- stepAIC(lm(TGFa ~ RA_status * (Menopausal_bin + Clinic + Diet_bin + Lifetime_Sexual_Partners_bin), data = your_data_immune_fin))
summary(best_model)
emmeans(best_model, pairwise ~ Menopausal_bin, adjust = "tukey")

model_bestX <- lm(formula = IL.12.p40. ~ RA_status + Menopausal_bin + Clinic + RA_status:Menopausal_bin + RA_status:Clinic + Lifetime_Sexual_Partners_bin + Diet_bin,
                  data = your_data_immune_fin)
summary(model_bestX)
#FINISH OF TESTING
emmeans(model_bestX, pairwise ~ Menopausal_bin, adjust = "tukey")


#Mann whitney to find notable cytokines in general
library(dplyr)
# Load your data (replace 'your_file.csv' with your actual filename)
cytokine_data <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few.csv")
# Ensure RA_status_cat is a factor
cytokine_data$RA_status_cat <- as.factor(cytokine_data$RA_status_cat)
# Filter columns prefixed with 'cyt_'
cytokine_columns <- grep("^cyt_", colnames(cytokine_data), value = TRUE)
# Initialize a results dataframe
mann_whitney_results <- data.frame(
  Cytokine = character(),
  W_Statistic = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)
# Loop through each cytokine column and perform Mann-Whitney U test
for (cytokine in cytokine_columns) {
  # Perform the test
  test_result <- wilcox.test(
    cytokine_data[[cytokine]] ~ cytokine_data$RA_status_cat,
    data = cytokine_data,
    exact = FALSE  # Use asymptotic approximation if data is large
  )
  # Save the results
  mann_whitney_results <- mann_whitney_results %>%
    add_row(
      Cytokine = cytokine,
      W_Statistic = test_result$statistic,
      P_Value = test_result$p.value
    )
}
# View the results
print(mann_whitney_results)
# Save the results to a file if needed
write.csv(mann_whitney_results, "Summary/RAefCyto_Mann_whitney_results.csv", row.names = FALSE)





#### Immune plots PCoA #### 
# Load necessary libraries
install.packages("vegan")
library(vegan)    # For distance calculation and PERMANOVA
#[vegan 2.6-8)]
library(ggplot2)  # For plotting
library(dplyr)    # For data manipulation
library(plyr)

cytokine_data_samples <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few.csv") #alt noNAN

# 1. Prepare the data
# Extract cytokine columns (e.g., columns prefixed with 'cyt_')
cytokine_data <- cytokine_data_samples %>% dplyr::select(starts_with("cyt_"))
metadata <- cytokine_data_samples$RA_status  # Replace with your metadata column
metadata_diet <- cytokine_data_samples$Diet_RAbin  # Replace with your metadata column
metadata_menopause <- cytokine_data_samples$RA_menopausebin  # Replace with your metadata column
metadata_LSP <- cytokine_data_samples$Lifetime_Sexual_Partners_binRA  # Replace with your metadata column

log_scale_cytokines = log(as.data.frame(cytokine_data))

# Ensure no missing values in the data
log_scale_cytokines <- na.omit(log_scale_cytokines)

# 2. Calculate distance matrix (e.g., euclidean, Bray-Curtis)
dist_matrix <- vegdist(log_scale_cytokines, method = "euclidean")

# 3. Perform PCoA
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)  # k = 2 for 2D plot
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PC1", "PC2")

# 4. Perform PERMANOVA
permanova_result <- adonis2(dist_matrix ~ metadata, data = log_scale_cytokines)

# Extract p-value and R-squared
p_value <- permanova_result$`Pr(>F)`[1]
r_squared <- permanova_result$R2[1]

# 5. Merge PCoA points with metadata
pcoa_points$Group <- metadata

# 6. Plot the PCoA
pcoa_plot <- ggplot(pcoa_points, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 1, linetype = "solid", alpha = 0.5) +  # Customize ellipses
  theme_minimal() +
  labs(
    title = "PCoA of Cytokine Profiles",
    x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")
  ) +
  theme(axis.title = element_text(family="Arial", size=14 ),
        axis.text = element_text(family = "Arial", size = 12),
  ) +
  theme(legend.position = "bottom",
        legend.text = element_text(family="Arial", size = 12), 
        legend.title.position = "bottom",
        legend.title = element_text(family="Arial", size = 14)) +
  scale_color_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +
  annotate(
    "text", x = max(pcoa_points$PC1) * 1.2, y = max(pcoa_points$PC2)*1.05,
    label = paste0("p-value: ", format.pval(p_value), "\nR²: ", round(r_squared, 3)),
    hjust = 1, vjust = 1, size = 4, color = "black"
  )
# Print the plot
print(pcoa_plot)

# Repeat for diet status
permanova_result <- adonis2(dist_matrix ~ metadata_diet, data = log_scale_cytokines)
p_value <- permanova_result$`Pr(>F)`[1]
r_squared <- permanova_result$R2[1]
metadata_diet <- revalue(metadata_diet, c("Control_High" = "Control: High",
                        "Control_Reduced" = "Control: Reduced",
                        "Rheumatoid Arthritis_High" = "RA: High",
                        "Rheumatoid Arthritis_Reduced" = "RA: Reduced"))
pcoa_points$Diet <- metadata_diet
pcoa_plot <- ggplot(pcoa_points, aes(x = PC1, y = PC2, color = Diet)) +
  geom_point(size = 4, alpha = 0.8) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 1, linetype = "solid", alpha = 0.5) +  # Customize ellipses
  theme_minimal() +
  labs(
    title = "PCoA of Cytokine Profiles",
    x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")
  ) +
  theme(axis.title = element_text(family="Arial", size=14 ),
        axis.text = element_text(family = "Arial", size = 12),
  ) +
  theme(legend.position = "bottom",
        legend.text = element_text(family="Arial", size = 12), 
        legend.title = element_text(family="Arial", size = 14),
        legend.title.position = "bottom") +
scale_color_manual(values = c("Control: High" = "#9ab7b1",
                                "Control: Reduced" = "#cbcbc1",
                                "RA: High" = "#8b0c43",
                                "RA: Reduced" = "#f1aecb")) +
  annotate(
    "text", x = max(pcoa_points$PC1) * 1.5, y = max(pcoa_points$PC2)*1.05,
    label = paste0("p-value: ", format.pval(p_value), "\nR²: ", round(r_squared, 3)),
    hjust = 1, vjust = 1, size = 4, color = "black"
  )
# Print the plot
print(pcoa_plot)


# Repeat for menopausal status
permanova_result <- adonis2(dist_matrix ~ metadata_menopause, data = log_scale_cytokines)
p_value <- permanova_result$`Pr(>F)`[1]
r_squared <- permanova_result$R2[1]
metadata_menopause <- revalue(metadata_menopause, c("Control_Pre-menopausal" = "Control: Pre-",
                                          "Control_Peri/Menopause" = "Control: Peri-/",
                                          "Control_Post-menopausal" = "Control: Post-",
                                          "Control_Post-menopause" = "Control: Post-",
                                          "Rheumatoid Arthritis_Pre-menopausal" = "RA: Pre-",
                                          "Rheumatoid Arthritis_Peri/Menopause" = "RA: Peri-/",
                                          "Rheumatoid Arthritis_Post-menopausal" = "RA: Post-"))

pcoa_points$Menopause <- metadata_menopause
pcoa_points$Menopause <- factor(pcoa_points$Menopause, 
                                            levels = c("Control: Pre-",
                                                       "Control: Peri-/",
                                                       "Control: Post-",
                                                       "RA: Pre-",
                                                       "RA: Peri-/",
                                                       "RA: Post-"))

pcoa_plot <- ggplot(pcoa_points, aes(x = PC1, y = PC2, color = Menopause)) +
  geom_point(size = 4, alpha = 0.8) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 1, linetype = "solid", alpha = 0.5) +  # Customize ellipses
  theme_minimal() +
  labs(
    title = "PCoA of Cytokine Profiles",
    x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")
  ) +
  theme(axis.title = element_text(family="Arial", size=14 ),
        axis.text = element_text(family = "Arial", size = 12),
  ) +
  theme(legend.position = "bottom",
        legend.text = element_text(family="Arial", size = 12), 
        legend.justification = "right",
        legend.title = element_text(family="Arial", size = 14),
        legend.title.position = "bottom") +
scale_color_manual(values = c("Control: Pre-" = "#d5e5e2",
                                "Control: Peri-/" = "#abcbc5",
                                "Control: Post-" = "#566663",
                                "RA: Pre-" = "#e19eb4",
                                "RA: Peri-/" = "#8b0c43",
                                "RA: Post-" = "#460622"
                               )) +
  annotate(
    "text", x = max(pcoa_points$PC1) * 1.3, y = max(pcoa_points$PC2)*1.435,
    label = paste0("p-value: ", format.pval(p_value), "\nR²: ", round(r_squared, 3)),
    hjust = 1, vjust = 1, size = 4, color = "black"
  )
# Print the plot
print(pcoa_plot)



# Repeat for lsp status
permanova_result <- adonis2(dist_matrix ~ metadata_LSP, data = log_scale_cytokines)
p_value <- permanova_result$`Pr(>F)`[1]
r_squared <- permanova_result$R2[1]
pcoa_points$LSP <- metadata_LSP
pcoa_plot <- ggplot(pcoa_points, aes(x = PC1, y = PC2, color = LSP)) +
  geom_point(size = 4, alpha = 0.8) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 1, linetype = "solid", alpha = 0.5) +  # Customize ellipses
  theme_minimal() +
  labs(
    title = "PCoA of Cytokine Profiles",
    x = paste0("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")
  ) +
  theme(axis.title = element_text(family="Arial", size=14 ),
        axis.text = element_text(family = "Arial", size = 12),
  ) +
  theme(legend.position = "bottom",
        legend.title.position = "bottom",
        legend.text = element_text(family="Arial", size = 12), 
        legend.title = element_text(family="Arial", size = 14)) +
  scale_color_manual(values = c("Control: 0" = "#d5e5e2",
                                "Control: 1 or 2" = "#abcbc5",
                                "Control: 3+" = "#566663",
                                "RA: 0" = "#e19eb4",
                                "RA: 1 or 2" = "#8b0c43",
                                "RA: 3+" = "#460622"
  ))+
  annotate(
    "text", x = max(pcoa_points$PC1) * 2.2, y = max(pcoa_points$PC2)*1.05,
    label = paste0("p-value: ", format.pval(p_value), "\nR²: ", round(r_squared, 3)),
    hjust = 1, vjust = 1, size = 4, color = "black"
  )
# Print the plot
print(pcoa_plot)



ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = RACST_data)) +
  geom_point(size = 4) +  # Plot the points
  stat_ellipse(type = "t", level = 0.95, linewidth = 1, linetype = "solid", alpha = 0.5) +  # Customize ellipses
  theme_minimal() +
  labs(
    title = paste("PCoA of Cytokine Data"),
    x = paste("PCoA1 (", round(variance_explained[1,1], 2), "%)"), y = paste("PCoA2 (", round(variance_explained[2,1], 2), "%)")) +
  theme(axis.title = element_text(family="Arial", size=16),
        axis.text = element_text(family = "Arial", size = 14),
        legend.text = element_text(family="Arial", size = 14), 
        legend.title = element_text(family="Arial", size = 14),
        legend.justification = "top") +
  scale_color_manual(values = c("Control_I" = "gold",
                                "Control_II" = "tan",
                                "Control_III" = "orange",
                                "Control_IV-A" = "lavender",
                                "Control_IV-C" = "green",
                                "Rheumatoid Arthritis_I" = "darkgoldenrod3",
                                "Rheumatoid Arthritis_II" = "brown",
                                "Rheumatoid Arthritis_III" = "darkorange2",
                                "Rheumatoid Arthritis_IV-A" = "plum",
                                "Rheumatoid Arthritis_IV-B" = "pink",
                                "Rheumatoid Arthritis_IV-C" = "darkgreen",
                                "Rheumatoid Arthritis_V" = "deepskyblue3"
  ))
)



#### Immune plots density #### 

cytokine_data <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few.csv")
# Reshape data to long format, where each cytokine is now a row
cytokine_data_long <- cytokine_data %>%
  pivot_longer(cols = starts_with("cyt"),  # Select cytokine columns (e.g., IL6, TNF, IL10)
               names_to = "Cytokine", values_to = "Value")
cytokine_data_long$Cytokine <- gsub("cyt_","",as.character(cytokine_data_long$Cytokine))
# Count occurrences of each cytokine value within RA status and Menopausal group
cytokine_counts <- cytokine_data_long %>%
  group_by(RA_status, Menopausal_bin, Cytokine, Value) %>%
  tally(name = "Count")  # Count occurrences of each cytokine value in each group
# View the resulting counts
head(cytokine_counts)
# Save the counts to a new CSV file
write.csv(cytokine_counts, "cytokine_counts_by_menopause.csv", row.names = FALSE)

# Select data for multiple cytokines, e.g., "IL6", "TNF", and "IL10"
#cytokine_data_long_selected <- cytokine_data_long %>%
#filter(Cytokine == "IL.a"))  # Use %in% to select multiple cytokines

# Select data for multiple cytokines, e.g., "IL6", "TNF", and "IL10"
cytokine_data_long_set1 <- cytokine_data_long %>%
  filter(Cytokine %in% c("IL.1RA", "IL.8", "IL.1a", "MIG.CXCL9", "M.CSF", "IL.18", "IL.1b"))  # Use %in% to select multiple cytokines
cytokine_data_long_set2 <- cytokine_data_long %>%
  filter(Cytokine %in% c("MCP.1", "G.CSF", "EGF", "RANTES", "VEGF.A", "IP.10", "PDGF.AB.BB", "FGF.2"))
cytokine_data_long_set3 <- cytokine_data_long %>%
  filter(Cytokine %in% c("IFNg", "sCD40L", "IL.27", "MIP.1a", "IL.6", "IL.9", "MCP.3", "Fractalkine", "Eotaxin", "IFNa2", "MDC", "TGFa", "IL.12.p40."))
cytokine_data_long_set4 <- cytokine_data_long %>%
  filter(Cytokine %in% c("GM.CSF", "TNFa", "PDGF.AA", "IL.13", "IL.10", "IL.17F", "FLT.3L", "IL.12.p70.", "IL.17A", "IL.15", "TNFb", "IL.3", "IL.4", "IL.2", "IL.5"))

cytokine_data_long_MENOmodels_sigs <- cytokine_data_long %>%
  filter(Cytokine %in% c("IFNa2",
                         "IFNg",
                         "IL.12.p40.",
                         "TGFa",
                          "IL.18",
                         "MCP.1",
                         "sCD40L"
                         ))
cytokine_data_long_MENOsigs <- cytokine_data_long %>%
  filter(Cytokine %in% c("Eotaxin",
                         "GM.CSF",
                         "IL.2",
                         "IL.9",
                         "IL.10",
                         "IL.22",
                         "IL.27"
  ))
cytokine_data_long_DIETmodels_sigs <- cytokine_data_long %>%
  filter(Cytokine %in% c("sCD40L",
                         "IL.9",
                         "IL.10",
                         "IL.27"
  ))
cytokine_data_long_DIETsigs <- cytokine_data_long %>%
  filter(Cytokine %in% c(
                         "Eotaxin",
                         "GM.CSF",
                         "IL.1a",
                         "IL.8",
                         "IL.22",
                         "TGFa"
  ))
cytokine_data_long_MENOsigs$Menopausal_bin <- factor(cytokine_data_long_MENOsigs$Menopausal_bin, 
                                                 levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
cytokine_data_long_MENOmodels_sigs$Menopausal_bin <- factor(cytokine_data_long_MENOmodels_sigs$Menopausal_bin, 
                                                           levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))

#TESTING
cytokine_data_RAvC_menodietTESTING <- cytokine_data_long %>%
  filter(Cytokine %in% c("sCD40L", "EGF",
                         "Eotaxin", 
                         "GM.CSF", 
                         "IFNa2", 
                         "IFNg", 
                         "IL.1a", 
                         "IL.1RA",
                         "IL.2",
                         "IL.8",
                         "IL.9",
                         "IL.10",
                         "IL.12p40",
                         "IL.18",
                         "IL.22",
                         "IL.27",
                         "MCP.1",
                         "M.CSF", 
                         "TGF.a"))
cytokine_data_RAvC_menodietTESTING$Menopausal_bin <- factor(cytokine_data_RAvC_menodietTESTING$Menopausal_bin, 
                                                            levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
ggplot(cytokine_data_RAvC_menodietTESTING, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=17, nrow = 2) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Diet bin and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options
ggplot(cytokine_data_RAvC_menodietTESTING, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Menopausal_bin + Cytokine, ncol=16, nrow = 3) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options

quartz()        
#Menopausal_bin_ordered
#Lifetime_Sexual_Partners_bin
# Create the density plots, faceted by Menopausal_status
ggplot(cytokine_data_long_set1, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=7, nrow = 2) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options
ggplot(cytokine_data_long_set2, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=7, nrow = 2) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options
ggplot(cytokine_data_long_set3, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=11, nrow = 2) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options
ggplot(cytokine_data_long_set4, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=11, nrow = 2) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options


ggplot(cytokine_data_long_MENOsigs, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Menopausal_bin + Cytokine, ncol=7, nrow = 3) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options
ggplot(cytokine_data_long_MENOmodels_sigs, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Menopausal_bin + Cytokine, ncol=7, nrow = 3) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options

ggplot(cytokine_data_long_DIETmodels_sigs, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=4, nrow = 3) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options
ggplot(cytokine_data_long_DIETsigs, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=6, nrow = 3) +  # Facet by Menopausal_status, 3 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal Status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options




cytokine_data_long_bestmodelx_DIET <- cytokine_data_long %>%
  filter(Cytokine %in% c("sCD40L",  
                         "Eotaxin",
                         "IL.8",
                         "IL.10",
                         "IL.27"))
cytokine_data_long_bestmodelx_DIET$Cytokine <- factor(cytokine_data_long_bestmodelx_DIET$Cytokine, 
          levels = c("sCD40L","Eotaxin","IL.8","IL.10","IL.27"))
cytokine_data_long_bestmodelx_MENO <- cytokine_data_long %>%
  filter(Cytokine %in% c("GM.CSF",  
                         "IL.2",
                         "IL.12.p40.",
                         "TGFa"))
cytokine_data_long_bestmodelx_MENO$Cytokine <- factor(cytokine_data_long_bestmodelx_MENO$Cytokine, 
           levels = c("GM.CSF","IL.2","IL.12.p40.","TGFa"))
cytokine_data_long_bestmodelx_MENO$Menopausal_bin <- factor(cytokine_data_long_bestmodelx_MENO$Menopausal_bin, 
           levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
cytokine_data_long_bestmodelx_DIETMENO <- cytokine_data_long %>%
  filter(Cytokine %in% c("IFNg",  
                         "IL.1a",
                         "IL.9",
                         "IL.22"))
cytokine_data_long_bestmodelx_DIETMENO$Cytokine <- factor(cytokine_data_long_bestmodelx_DIETMENO$Cytokine, 
           levels = c("IFNg","IL.1a","IL.9","IL.22"))
cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin <- revalue(cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin, c("Control_High" = "Control: High",
  "Control_Reduced" = "Control: Reduced",
  "Rheumatoid Arthritis_High" = "RA: High",
  "Rheumatoid Arthritis_Reduced" = "RA: Reduced"))
cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin <- factor(cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin,
           levels = c("Control: Reduced", "Control: High", "RA: Reduced", "RA: High"))
cytokine_data_long_bestmodelx_DIETMENO$Menopausal_bin <- factor(cytokine_data_long_bestmodelx_DIETMENO$Menopausal_bin, 
            levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))

ggplot(cytokine_data_long_bestmodelx_DIET, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Diet_bin + Cytokine, ncol=5, nrow = 2) +  # Facet by Diet_bin, 5 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Diet and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options

ggplot(cytokine_data_long_bestmodelx_MENO, aes(x = Value, fill = RA_status)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Menopausal_bin + Cytokine, ncol=4, nrow = 3) +  # Facet by Diet_bin, 5 columns
  scale_fill_manual(values = c("Control" = "#abcbc5", "Rheumatoid Arthritis" = "#8b0c43")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Menopausal status and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options

ggplot(cytokine_data_long_bestmodelx_DIETMENO, aes(x = Value, fill = Diet_RAbin)) +
  geom_density(alpha = 0.5) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~ Menopausal_bin + Cytokine, ncol=4, nrow = 3) +  # Facet by Diet_bin, 5 columns
  scale_fill_manual(values = c("Control: High" = "#9ab7b1",
                                 "Control: Reduced" = "#cbcbc1",
                                 "RA: High" = "#8b0c43",
                                 "RA: Reduced" = "#f1aecb")) +
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by Diet and RA Status") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options

##### RA Immune box plots diet, meno, diet + meno ####

cytokine_data <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few.csv")
# Reshape data to long format, where each cytokine is now a row
cytokine_data_long <- cytokine_data %>%
  pivot_longer(cols = starts_with("cyt"),  # Select cytokine columns (e.g., IL6, TNF, IL10)
               names_to = "Cytokine", values_to = "Value")
cytokine_data_long$Cytokine <- gsub("cyt_","",as.character(cytokine_data_long$Cytokine))

cytokine_data_long_bestmodelx_DIET <- cytokine_data_long %>%
  filter(Cytokine %in% c("sCD40L",  
                         "Eotaxin",
                         "IL.8",
                         "IL.10",
                         "IL.27"))
cytokine_data_long_bestmodelx_DIET$Cytokine <- factor(cytokine_data_long_bestmodelx_DIET$Cytokine, 
                                                      levels = c("sCD40L","Eotaxin","IL.8","IL.10","IL.27"))
cytokine_data_long_bestmodelx_MENO <- cytokine_data_long %>%
  filter(Cytokine %in% c("GM.CSF",  
                         "IL.2",
                         "IL.12.p40.",
                         "TGFa"))
cytokine_data_long_bestmodelx_MENO$Cytokine <- factor(cytokine_data_long_bestmodelx_MENO$Cytokine, 
                                                      levels = c("GM.CSF","IL.2","IL.12.p40.","TGFa"))
cytokine_data_long_bestmodelx_MENO$Menopausal_bin <- factor(cytokine_data_long_bestmodelx_MENO$Menopausal_bin, 
                                                            levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
cytokine_data_long_bestmodelx_DIETMENO <- cytokine_data_long %>%
  filter(Cytokine %in% c("IFNg",  
                         "IL.1a",
                         "IL.9",
                         "IL.22"))
cytokine_data_long_bestmodelx_DIETMENO$Cytokine <- factor(cytokine_data_long_bestmodelx_DIETMENO$Cytokine, 
                                                          levels = c("IFNg","IL.1a","IL.9","IL.22"))
cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin <- revalue(cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin, c("Control_High" = "Control: High",
                                                                                                                  "Control_Reduced" = "Control: Reduced",
                                                                                                                  "Rheumatoid Arthritis_High" = "RA: High",
                                                                                                                  "Rheumatoid Arthritis_Reduced" = "RA: Reduced"))
cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin <- factor(cytokine_data_long_bestmodelx_DIETMENO$Diet_RAbin,
                                                            levels = c("Control: Reduced", "RA: Reduced", "Control: High", "RA: High"))
cytokine_data_long_bestmodelx_DIETMENO$Menopausal_bin <- factor(cytokine_data_long_bestmodelx_DIETMENO$Menopausal_bin, 
                                                                levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
cytokine_data_long_bestmodelx_DIET$RA_status <- factor(cytokine_data_long_bestmodelx_DIET$RA_status, 
                                                      levels = c("Control","Rheumatoid Arthritis"))
cytokine_data_long_bestmodelx_DIET$RA_status <- revalue(cytokine_data_long_bestmodelx_DIET$RA_status, c("Control" = "Control", "Rheumatoid Arthritis" = "RA"))
cytokine_data_long_bestmodelx_DIET$RA_status <- factor(cytokine_data_long_bestmodelx_DIET$RA_status, 
                                                       levels = c("Control","RA"))
cytokine_data_long_bestmodelx_MENO$RA_status <- factor(cytokine_data_long_bestmodelx_MENO$RA_status, 
                                                       levels = c("Control","Rheumatoid Arthritis"))
cytokine_data_long_bestmodelx_MENO$RA_status <- revalue(cytokine_data_long_bestmodelx_MENO$RA_status, c("Control" = "Control", "Rheumatoid Arthritis" = "RA"))
cytokine_data_long_bestmodelx_DIETMENO$RA_status <- factor(cytokine_data_long_bestmodelx_DIETMENO$RA_status, 
                                                       levels = c("Control","Rheumatoid Arthritis"))
cytokine_data_long_bestmodelx_DIETMENO$RA_status <- revalue(cytokine_data_long_bestmodelx_DIETMENO$RA_status, c("Control" = "Control", "Rheumatoid Arthritis" = "RA"))

# Load required libraries
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(dplyr)


# Ensure RA_status is a factor
cytokine_data_long_bestmodelx_DIET$RA_status <- factor(
  cytokine_data_long_bestmodelx_DIET$RA_status, levels = c("Control", "RA")
)

# Define colors
color_mapRA <- c("Control" = "#abcbc5", "RA" = "#8b0c43")  # Teal & Burgundy
color_mapRA_diet <- c("Control: High" = "#9ab7b1",
                      "Control: Reduced" = "#cbcbc1",
                      "RA: High" = "#8b0c43",
                      "RA: Reduced" = "#f1aecb")
ggplot(cytokine_data_long_bestmodelx_DIET, aes(x = RA_status, y = Value, fill = RA_status)) +
  scale_fill_manual(values = c("Control" = "#abcbc5", "RA" = "#8b0c43")) +  
  scale_color_manual(values = c("Control" = "#abcbc5", "RA" = "#8b0c43")) + 
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options


# Perform Mann-Whitney test (Wilcoxon rank-sum test)
stat_tests <- cytokine_data_long_bestmodelx_DIET %>%
  group_split(Cytokine, Diet_bin) %>%
  map_df(~ {
    test_result <- tryCatch(
      wilcox.test(Value ~ RA_status, data = .x)$p.value,
      error = function(e) NA  # If test fails, return NA
    )
    
    tibble(
      Cytokine = unique(.x$Cytokine),
      Diet_bin = unique(.x$Diet_bin),
      p_value = test_result
    )
  }) %>%
  mutate(p_label = paste0("p = ", format.pval(p_value, digits = 3, eps = 0.001)))

print(stat_tests)

plot <- ggplot(cytokine_data_long_bestmodelx_DIET, aes(x = RA_status, y = Value, fill = RA_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.75, color = "black", fill = "white", linewidth = 0.5) +  # White boxes, thin borders
  geom_jitter(aes(color = RA_status, fill = RA_status), shape = 21,  
              size = 2.75, alpha = 0.7, width = 0.25, stroke = 0.5, color = "black") +  # Thin black outline for dots
  scale_fill_manual(values = color_mapRA) +  
  scale_color_manual(values = color_mapRA) +  
  scale_y_log10() +  
  facet_wrap(Diet_bin ~ Cytokine, nrow = 2) +  
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.x = element_text(size = 14, family = "Arial", color = "black"),
    axis.text = element_text(size = 14, family = "Arial", color = "black"),
    panel.spacing = unit(0.5, "lines")  
  ) +
  labs(
    title = "Cytokine Levels by RA and Diet",
    x = "Group",
    y = "Cytokine Concentration (log scale)"
  )
# Only add statistical annotations if there are significant results
if (nrow(stat_tests) > 0) {
  plot <- plot +
    geom_text(
      data = stat_tests,
      aes(x = 1.5, y = max(Value, na.rm = TRUE) * 1.2, label = p_label),
      inherit.aes = FALSE,
      size = 4,
      vjust = 1
    ) +
    geom_segment(
      data = stat_tests,
      aes(x = 1, xend = 2, y = max(Value, na.rm = TRUE) * 1.15, yend = max(Value, na.rm = TRUE) * 1.15),
      inherit.aes = FALSE,
      color = "black"
    )
}
# Print the final plot
print(plot)


# Perform Mann-Whitney test (Wilcoxon rank-sum test)
stat_tests_meno <- cytokine_data_long_bestmodelx_MENO %>%
  group_split(Cytokine, Menopausal_bin) %>%
  map_df(~ {
    test_result <- tryCatch(
      wilcox.test(Value ~ RA_status, data = .x)$p.value,
      error = function(e) NA  # If test fails, return NA
    )
    
    tibble(
      Cytokine = unique(.x$Cytokine),
      Menopausal_bin = unique(.x$Menopausal_bin),
      p_value = test_result
    )
  }) %>%
  mutate(p_label = paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001)))

print(stat_tests_meno)

#  filter(!is.na(p_label))  # Remove non-significant values
# Base plot
plot <- ggplot(cytokine_data_long_bestmodelx_MENO, aes(x = RA_status, y = Value, fill = RA_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.75, color = "black", fill = "white", linewidth = 0.5) +  # White boxes, thin borders
  geom_jitter(aes(color = RA_status, fill = RA_status), shape = 21,  
              size = 2.75, alpha = 0.7, width = 0.25, stroke = 0.5, color = "black") +  # Thin black outline for dots
  scale_fill_manual(values = color_mapRA) +  
  scale_color_manual(values = color_mapRA) +  
  scale_y_log10() +  
  facet_wrap(Menopausal_bin ~ Cytokine, nrow = 3) +  
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.x = element_text(size = 14, family = "Arial", color = "black"),
    axis.text = element_text(size = 14, family = "Arial", color = "black"),
    panel.spacing = unit(0.5, "lines")  
  ) +
  labs(
    title = "Cytokine Levels by RA and Menopausal status",
    x = "Group",
    y = "Cytokine Concentration (log scale)"
  )
# Only add statistical annotations if there are significant results
if (nrow(stat_tests_meno) > 0) {
  plot <- plot +
    geom_text(
      data = stat_tests_meno,
      aes(x = 1.5, y = max(Value, na.rm = TRUE) * 1.2, label = p_label),
      inherit.aes = FALSE,
      size = 4,
      vjust = 1
    ) +
    geom_segment(
      data = stat_tests_meno,
      aes(x = 1, xend = 2, y = max(Value, na.rm = TRUE) * 1.15, yend = max(Value, na.rm = TRUE) * 1.15),
      inherit.aes = FALSE,
      color = "black"
    )
}
# Print the final plot
print(plot)



# Perform Mann-Whitney test (Wilcoxon rank-sum test)
stat_tests_DM <- cytokine_data_long_bestmodelx_DIETMENO %>%
  group_split(Cytokine, Diet_bin, Menopausal_bin) %>%
  map_df(~ {
    test_result <- tryCatch(
      wilcox.test(Value ~ RA_status, data = .x)$p.value,
      error = function(e) NA  # Handle errors gracefully
    )
    
    tibble(
      Cytokine = unique(.x$Cytokine),
      Diet_bin = unique(.x$Diet_bin),
      Menopausal_bin = unique(.x$Menopausal_bin),
      p_value = test_result
    )
  }) %>%
  filter(p_value < 0.05) %>%  # Keep only significant comparisons
  mutate(p_label = paste0("p = ", format.pval(p_value, digits = 3, eps = 0.001)))

# Base plot
plot <- ggplot(cytokine_data_long_bestmodelx_DIETMENO, aes(x = Diet_RAbin, y = Value, fill = Diet_RAbin)) +
  geom_boxplot(outlier.shape = NA, width = 0.75, color = "black", fill = "white", linewidth = 0.5) +  # White boxes, thin borders
  geom_jitter(aes(color = Diet_RAbin, fill = Diet_RAbin), shape = 21,  
              size = 2.75, alpha = 0.7, width = 0.25, stroke = 0.5, color = "black") +  # Thin black outline for dots
  scale_fill_manual(values = color_mapRA_diet) +  
  scale_color_manual(values = color_mapRA_diet) +  
  scale_y_log10() +  
  facet_wrap(Menopausal_bin ~ Cytokine, nrow =3) +  
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.x = element_text(size = 14, family = "Arial", color = "black"),
    axis.text = element_text(size = 14, family = "Arial", color = "black"),
    panel.spacing = unit(0.5, "lines")  
  ) +
  labs(
    title = "Cytokine Levels by RA and Diet",
    x = "Group",
    y = "Cytokine Concentration (log scale)"
  )

# Print the final plot
print(plot)






#############################
#### RA ONLY IMMUNE  ####
#Diet, Men ---- ACCP
head(your_data_4factors[158:205])
cytokines <- colnames(your_data_4factors)[158:205]
factors <- c("Joint_cat_space_narrowing", "Radiographic_cat_Erosions")
factors_2 <- c("AntiCCP_cat_Ab_positive", "RF_cat_positive")
p_value_cytmatrix_RA2 <- matrix(NA, nrow = length(factors_2), ncol = length(cytokines),
                            dimnames = list(factors_2, cytokines))
coef_cytmatrix_RA2 <- matrix(NA, nrow = length(factors_2), ncol = length(cytokines),
                         dimnames = list(factors_2, cytokines))

#repeat replacing factors_2 for factors
for (factor in factors_2) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family = binomial(link ="logit"), data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA2[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix_RA2[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA2, "Summary/p_value_cytRAmatrix2_cytefACCPRF.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA2, "Summary/coef_cytRAmatrix2_cytefACCPRF.csv", row.names = TRUE)

#revise factors and factors_2, rerun matrix, then analysis before saving as below
p_value_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                               dimnames = list(factors, cytokines))
coef_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                            dimnames = list(factors, cytokines))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat + AntiCCP_cat_Ab_positive +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family = binomial(link ="logit"), data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefJSNRE.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefJSNRE.csv", row.names = TRUE)

#poisson AntiCCP_cat_Ab_positive - not good use of poisson, exclude
for (factor in factors_2) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat + AntiCCP_cat_Ab_positive +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family = poisson(link ="log"), data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefJSNRE_poiss.csv", row.names = TRUE) #exclude
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefJSNRE_poiss.csv", row.names = TRUE) #exclude
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefACCPRF_poiss.csv", row.names = TRUE) #exclude
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefACCPRF_poiss.csv", row.names = TRUE) #exclude

#SWITCH TO FACTOR FX CYTOKINE
p_value_cytmatrix_RAcyt <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                                  dimnames = list(cytokines, factors))
coef_cytmatrix_RAcyt <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                               dimnames = list(cytokines, factors))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Menopausal_bin_cat + AntiCCP_cat_Ab_positive +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RAcyt[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_cytmatrix_RAcyt[cytokine, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RAcyt, "Summary/p_value_cytRAmatrix2_JSNREefcyt.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RAcyt, "Summary/coef_cytRAmatrix2_JSNREefcyt.csv", row.names = TRUE)

p_value_cytmatrix_RAcyt <- matrix(NA, nrow = length(cytokines), ncol = length(factors_2),
                                dimnames = list(cytokines, factors_2))
coef_cytmatrix_RAcyt <- matrix(NA, nrow = length(cytokines), ncol = length(factors_2),
                             dimnames = list(cytokines, factors_2))
for (factor in factors_2) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Menopausal_bin_cat +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RAcyt[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_cytmatrix_RAcyt[cytokine, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RAcyt, "Summary/p_value_cytRAmatrix2_ACCPRFefcyt.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RAcyt, "Summary/coef_cytRAmatrix2_ACCPRFefcyt.csv", row.names = TRUE)



#CDAI
cytokines <- colnames(your_data_4factors)[158:205]
factors <- c("CDAI_score")
p_value_cytmatrix_RAx <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                                dimnames = list(cytokines, factors))
coef_cytmatrix_RAx <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                             dimnames = list(cytokines, factors))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Menopausal_bin_cat +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RAx[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_cytmatrix_RAx[cytokine, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RAx, "Summary/p_value_cytRAmatrix2_CDAIefcyt.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RAx, "Summary/coef_cytRAmatrix2_CDAIefcyt.csv", row.names = TRUE)

p_value_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                               dimnames = list(factors, cytokines))
coef_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                            dimnames = list(factors, cytokines))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|t|)"]
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefCDAI.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefCDAI.csv", row.names = TRUE)


#Troubleshoot out of bounds error
rownames(coef_summary)
if (factor %in% rownames(coef_summary)) {
  microbe[factor] <- coef_summary[factor, ]
} else {
  print("Factor not found in coef_summary")
}
any(is.na(factor))  # Check if the factor variable has any NA values
any(is.na(coef_summary)) 


coef_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                            dimnames = list(factors, cytokines))
stderror_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                                dimnames = list(factors, cytokines))
zvalue_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                              dimnames = list(factors, cytokines))
p_value_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                               dimnames = list(factors, cytokines))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family=poisson(link=log), data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
        stderror_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Std. Error"]
        zvalue_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "z value"]
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefCDAI_poiss_redo.csv", row.names = TRUE)
write.csv(stderror_cytmatrix_RA, "Summary/stderror_cytRAmatrix2_cytefCDAI_poiss_redo.csv", row.names = TRUE)
write.csv(zvalue_cytmatrix_RA, "Summary/zvalue_cytRAmatrix2_cytefCDAI_poiss_redo.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefCDAI_poiss_redo.csv", row.names = TRUE)

your_data_4factors_rounded <- read.csv("Metadata_RAstudy_master_RAcategorized_4factors_rounded.csv") #CDAI 40.5 rounded to 41
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family=poisson(link=log), data = your_data_4factors_rounded), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
        stderror_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Std. Error"]
        zvalue_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "z value"]
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefCDAI_poiss_redorounded.csv", row.names = TRUE)
write.csv(stderror_cytmatrix_RA, "Summary/stderror_cytRAmatrix2_cytefCDAI_poiss_redorounded.csv", row.names = TRUE)
write.csv(zvalue_cytmatrix_RA, "Summary/zvalue_cytRAmatrix2_cytefCDAI_poiss_redorounded.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefCDAI_poiss_redorounded.csv", row.names = TRUE)


factors <- c("Disease_Activity_cat")
p_value_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                               dimnames = list(factors, cytokines))
coef_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                            dimnames = list(factors, cytokines))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(glm(formula_cyt2, family=poisson(link=log), data = your_data_4factors), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|z|)"]
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefDisAct_poiss.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefDisAct_poiss.csv", row.names = TRUE)



cytokine_data_RAonly <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few_RAonly.csv")
# Reshape data to long format, where each cytokine is now a row
cytokine_data_long_RAonly <- cytokine_data_RAonly %>%
  pivot_longer(cols = starts_with("cyt"),  # Select cytokine columns (e.g., IL6, TNF, IL10)
               names_to = "Cytokine", values_to = "Value")
cytokine_data_long_RAonly$Cytokine <- gsub("cyt_","",as.character(cytokine_data_long_RAonly$Cytokine))

# Select data for multiple cytokines, e.g., "IL6", "TNF", and "IL10"
cytokine_data_CDAI_sigs <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("M.CSF", "TNFb",	"IL.27",	"IL.5",	"MIG.CXCL9",	"Eotaxin",	"IL.1b",	"MDC",	"TNFa",	"FLT.3L",	"MCP.3",	"IL.1a",	"MIP.1b",	"IL.12.p40.",	"sCD40L",	"IL.10",	"IFNa2",	"IP.10",	"IL.2",	"IL.12.p70.",	"IL.18",	"IFNg",	"IL.15",	"MCP.1",	"GM.CSF",	"IL.4",	"VEGF.A",	"RANTES",	"IL.6",	"TGFa",	"G.CSF",	"FGF.2"
  ))
#### RA ONLY IMMUNE plots cdai log ####
cytokines <- colnames(cytokine_data_CDAI_sigs)[362]
cytokine_data_CDAI_sigs$Predictions <- NA
for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste("CDAI_score ~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model
    model_cyt2 <- glm(formula_cyt2, family= poisson(log), data = cytokine_data_CDAI_sigs)
    # Generate predictions (for the given dataset)
    cytokine_data_CDAI_sigs$Predictions <- predict(model_cyt2, newdata = cytokine_data_CDAI_sigs, type = "response")
    }

gg <- ggplot(cytokine_data_CDAI_sigs, aes(x = Value, y = CDAI_score)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.5)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), colour = "black", linewidth = 0.5, alpha=0.4) +
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4) +  
  scale_x_log10() +
  scale_color_manual(values = c("High Activity" = "#460622", "Moderate Activity" = "#8b0c43", "Low Activity" = "#ae557b", "Remission" = "#dcb6c7")) +
  ylim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI", 
    x = "log(pg/mL)", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: RA and Controls"
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )
# Print the plot
print(gg)

#### RA ONLY IMmune plots cdai axes AND MEDS ####
cytokine_data_CDAI_sigslow <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("FLT.3L",	"IL.12.p70.",	"IL.2",	"IL.4",	"IL.5",	"TNFb"
  ))
cytokine_data_CDAI_sigsmed <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("Eotaxin", "IL.6", "IL.10",	"GM.CSF",	"IFNa2",	"IL.12.p40.",	"MCP.3",	"MDC",	"sCD40L",	"TGFa", "TNFa"
  ))
cytokine_data_CDAI_sigshigh <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("FGF.2",	"G.CSF",	"IFNg", "IL.18",	"IL.27", "IP.10",	"MCP.1", "VEGF.A",	"RANTES"
  ))
cytokine_data_CDAI_sigshigher <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("IL.1b", "M.CSF",	"MIG.CXCL9"
  ))
cytokine_data_CDAI_sigshighest <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("IL.1a"
  ))
cytokine_data_CDAI_sigs <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("M.CSF", "TNFb",	"IL.27",	"IL.5",	"MIG.CXCL9",	"Eotaxin",	"IL.1b",	"MDC",	"TNFa",	"FLT.3L",	"MCP.3",	"IL.1a",	"MIP.1b",	"IL.12.p40.",	"sCD40L",	"IL.10",	"IFNa2",	"IP.10",	"IL.2",	"IL.12.p70.",	"IL.18",	"IFNg",	"IL.15",	"MCP.1",	"GM.CSF",	"IL.4",	"VEGF.A",	"RANTES",	"IL.6",	"TGFa",	"G.CSF",	"FGF.2"))

gg <- ggplot(cytokine_data_CDAI_groupedPaper, aes(x = Value, y = CDAI_score)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.5)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), colour = "black", linewidth = 0.5, alpha=0.4) +
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4) +  
  scale_color_manual(values = c("High Activity" = "#460622", "Moderate Activity" = "#8b0c43", "Low Activity" = "#ae557b", "Remission" = "#dcb6c7")) +
  ylim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI", 
    x = "pg/mL", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title = element_text(size=10),
    legend.text = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )
# Print the plot
print(gg)

#general
gg <- ggplot(cytokine_data_CDAI_groupedPaper, aes(x = Value, y = CDAI_score)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.5)) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), colour = "black", linewidth = 0.5, alpha=0.4) +
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4) +  
  scale_color_manual(values = c("High Activity" = "#460622", "Moderate Activity" = "#8b0c43", "Low Activity" = "#ae557b", "Remission" = "#dcb6c7")) +
  ylim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI", 
    x = "log(pg/mL)", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: RA and Controls"
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )
# Print the plot
print(gg)





cytokine_data_RAonly <- read.csv("Metadata_RAstudy_categorized_cytokinesLABELED_few_RAonly.csv")
# Reshape data to long format, where each cytokine is now a row
cytokine_data_long_RAonly <- cytokine_data_RAonly %>%
  pivot_longer(cols = starts_with("cyt"),  # Select cytokine columns (e.g., IL6, TNF, IL10)
               names_to = "Cytokine", values_to = "Value")
cytokine_data_long_RAonly$Cytokine <- gsub("cyt_","",as.character(cytokine_data_long_RAonly$Cytokine))

#Grouped by appearance on paper BEAUTY
cytokine_data_CDAI_groupedPaper <- cytokine_data_long_RAonly %>%
  filter(Cytokine %in% c("M.CSF", "TNFb",	"IL.27",	"IL.5",	"MIG.CXCL9",	"Eotaxin",	"IL.1b",	"MDC",	"TNFa",	"FLT.3L",	"MCP.3",	"IL.1a",	"MIP.1b",	"IL.12.p40.",	"sCD40L",	"IL.10",	"IFNa2",	"IP.10",	"IL.2",	"IL.12.p70.",	"IL.18",	"IFNg",	"IL.15",	"MCP.1",	"GM.CSF",	"IL.4",	"VEGF.A",	"RANTES",	"IL.6",	"TGFa",	"G.CSF",	"FGF.2"))
cytokine_data_CDAI_groupedPaper$Cytokine <- factor(cytokine_data_CDAI_groupedPaper$Cytokine,
              levels = c("IL.18",	"TNFa",	"IL.12.p70.", "sCD40L", "Eotaxin", "IL.10", "IL.27",	"GM.CSF",	"IL.2", "IL.12.p40.", "TGFa",	"IFNg", "IL.1a", "MCP.3", "MDC", "TNFb",
"FGF.2", "FLT.3L", "G.CSF",	"IFNa2",	"IL.1b",	"IL.4",	"IL.5",	"IL.6",	"IL.15",	"IP.10", "MCP.1", "M.CSF", "MIG.CXCL9",	"MIP.1b",	"RANTES",	"VEGF.A"))
cytokine_data_CDAI_groupedPaper$Disease_Activity <- factor(cytokine_data_CDAI_groupedPaper$Disease_Activity,
              levels = c("Remission", "Low Activity", "Moderate Activity", "High Activity"))
cytokine_data_CDAI_groupedPaper$Cytokine <- revalue(cytokine_data_CDAI_groupedPaper$Cytokine,
                                                    c("IL.18"="IL-18",
                                                      "TNFa"="TNF",
                                                      "IL.12.p70."="IL-12p(70)",
                                                      "sCD40L"="sCD40L",
                                                      "Eotaxin"="Eotaxin",
                                                      "IL.10"="IL-10",
                                                      "IL.27"="IL-27",
                                                      "GM.CSF"="GM-CSF",
                                                      "IL.2"="IL-2",
                                                      "IL.12.p40."="IL-12p(40)",
                                                      "TGFa"="TGF-a",
                                                      "IFNg"="IFN-g",
                                                      "IL.1a"="IL-1a",
                                                      "MDC"="MDC",
                                                      "TNFb"="LT-a",
                                                      "MCP.3"="MCP-3",
                                                      "FGF.2"="FGF-2",
                                                      "FLT.3L"="FLT-3L",
                                                      "G.CSF"="G-CSF",
                                                      "IFNa2"="IFN-a2",
                                                      "IL.1b"="IL-1B",
                                                      "IL.4"="IL-4",
                                                      "IL.5"="IL-5",
                                                      "IL.6"="IL-6",
                                                      "IL.15"="IL-15",
                                                      "IP.10"="IP-10",
                                                      "MCP.1"="MCP-1",
                                                      "M.CSF"="M-CSF",
                                                      "MIG.CXCL9"="MIG (CXCL9)",
                                                      "MIP.1b"="MIP-1B",
                                                      "RANTES"="RANTES",
                                                      "VEGF.A"="VEGF-a"))
gg <- ggplot(cytokine_data_CDAI_groupedPaper, aes(x = CDAI_score, y = Value)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.7)) + 
  geom_smooth(
    method = "glm", 
    method.args = list(family = "poisson"), se= FALSE,
    colour = "black", 
    linewidth = 0.5, 
    alpha = 0.4
  ) +
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4, scales = "free_y") +  
  scale_color_manual(values = c(
    "High Activity" = "#380019", 
    "Moderate Activity" = "#af2460", 
    "Low Activity"= "#e89e9f",
    "Remission" = "#dcb6c7")) +
  xlim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "pg/mL", 
    x = "CDAI score", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: Influence of cytokines on predicting CDAI"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title= element_text(size=10),
    legend.text = element_text(size=8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 8)  # Size of facet labels
  )
print(gg)





#medications
cytokines <- colnames(df_meds)[158:197]
medications <- c("Hydroxychloroquine", "Sulfasalazine", "Methotrexate_dose_.mg.", "Prednisone_dose_.mg.", "Methotrexate", "Prednisone")
p_value_matrix <- matrix(NA, nrow = length(cytokines), ncol = length(medications),
                         dimnames = list(cytokines, medications))
coef_matrix <- matrix(NA, nrow = length(cytokines), ncol = length(medications),
                      dimnames = list(cytokines, medications))
stderror_matrix <- matrix(NA, nrow = length(cytokines), ncol = length(medications),
                      dimnames = list(cytokines, medications))
t_value_matrix <- matrix(NA, nrow = length(cytokines), ncol = length(medications),
                      dimnames = list(cytokines, medications))

# Loop through cytokines and medications
for (cytokine in cytokines) {
  for (med in medications) {
    # Build the formula
    formula <- as.formula(paste(cytokine, "~ Diet_bin + Clinic_cat + Menopausal_bin_ordered +", med))
    
    # Fit the model
    model <- lm(formula, data = df_meds)
    
    # Extract the p-value for the medication variable
    p_value <- summary(model)$coefficients[med, "Pr(>|t|)"]
    estimate_coef <- summary(model)$coefficients[med, "Estimate"]
    stderror <- summary(model)$coefficients[med, "Std. Error"]
    t_value <- summary(model)$coefficients[med, "t value"]
    
    # Store the p-value in the matrix
    p_value_matrix[cytokine, med] <- p_value
    coef_matrix[cytokine, med] <- estimate_coef
    stderror_matrix[cytokine, med] <- stderror
    t_value_matrix[cytokine, med] <- t_value
  }
}

# View the resulting p-value matrix
print(p_value_matrix)

# Optionally save to a CSV file
write.csv(p_value_matrix, "Summary/p_value_matrix_meds_yn.csv", row.names = TRUE)
write.csv(coef_matrix, "Summary/coef_matrix_meds_yn.csv", row.names = TRUE)
write.csv(stderror_matrix, "Summary/stderror_matrix_meds_yn.csv", row.names = TRUE)
write.csv(t_value_matrix, "Summary/tvalue_matrix_meds_yn.csv", row.names = TRUE)





model_med <- lm(IFNa2 ~
                  Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                data = df_meds)
#0.04064 (adjR 0.229)
model_Hydroxy <- lm(G.CSF ~
                      Hydroxychloroquine +
                      Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                    data = df_meds)
#0.04064 (adjR 0.229)model_Hydroxy <- lm(sCD40L ~
Hydroxychloroquine +
  Diet_bin + Clinic_cat + Menopausal_bin_ordered,
data = your_data_meds)
#0.04064 (adjR 0.229)
model_Sulfa <- lm(IFNa2 ~
                    Sulfasalazine +
                    Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                  data = df_meds)
#0.04064 (adjR 0.229)
model_Metdose <- lm(IFNa2 ~
                      Methotrexate_dose_.mg. +
                      Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                    data = df_meds)
#0.04064 (adjR 0.229)
model_Preddose <- lm(IFNa2 ~
                       Prednisone_dose_.mg. +
                       Diet_bin + Clinic_cat + Menopausal_bin_ordered,
                     data = df_meds)
#0.04064 (adjR 0.229)
summary(model_med)
summary(model_Hydroxy)
summary(model_Sulfa)
summary(model_Metdose)
summary(model_Preddose)


vif(model_Hydroxy)
plot(model_Hydroxy)



#### Increasing value Immune plots ####
# Step 1: Pre-calculate concentration ranges and slope directions for cytokines
cytokine_stats <- cytokine_data_CDAI_sigs %>%
  group_by(Cytokine) %>%
  summarise(
    concentration_range = max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE),
  ) %>%
  arrange(desc(concentration_range)) %>%
  mutate(Cytokine = factor(Cytokine, levels = Cytokine))

# Step 2: Update cytokine_data_CDAI_sigs with sorted Cytokine levels
cytokine_data_CDAI_sigs <- cytokine_data_CDAI_sigs %>%
  mutate(Cytokine = factor(Cytokine, levels = cytokine_stats$Cytokine))

# Step 3: Plot with facets having free scales and ordered by calculated stats
gg <- ggplot(cytokine_data_CDAI_sigs, aes(x = Value, y = CDAI_score)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.7)) + 
  geom_smooth(
    method = "glm", 
    method.args = list(family = "poisson"), 
    colour = "black", 
    linewidth = 0.5, 
    alpha = 0.4
  ) +
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4, scales = "free_x") +  
  scale_color_manual(values = c(
    "High Activity" = "#380019", 
    "Moderate Activity" = "#af2460", 
    "Low Activity"= "#e89e9f",
    "Remission" = "#dcb6c7")) +
  ylim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI score", 
    x = "pg/mL", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: RA and Controls"
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title= element_text(size=10),
    legend.text = element_text(size=8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 8)  # Size of facet labels
  )
# Print the plot
print(gg)
#2a0414
#601234
#5e002a
#EXTRA FACETSS AND FIT
gg <- ggplot(cytokine_data_CDAI_sigs, aes(x = CDAI_score, y = Value)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.7)) + 
  geom_smooth(
    method = "loess", 
    span=3, 
    colour = "black", 
    linewidth = 0.5, 
    alpha = 0.4
  ) +
  facet_wrap(~ Disease_Activity + Cytokine, ncol = 32, nrow = 4,) +  
  scale_color_manual(values = c(
    "High Activity" = "#380019", 
    "Moderate Activity" = "#af2460", 
    "Low Activity"= "#e89e9f",
    "Remission" = "#dcb6c7")) +
  xlim(c(0, 60)) +
  scale_y_log10() +
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI score", 
    x = "pg/mL", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: RA and Controls"
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title= element_text(size=10),
    legend.text = element_text(size=8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 8)  # Size of facet labels
  )
# Print the plot
print(gg)


ggplot(cytokine_data_CDAI_sigs, aes(y= Disease_Activity, x = Value, fill = Disease_Activity)) +
  geom_density_ridges(scale = 1.5, alpha = 0.5, ) +  # Add density plot, adjust alpha for transparency (also like 0.4)
  facet_wrap(~Cytokine, ncol=32) +
  scale_fill_manual(values = c("Remission" = "#FEFFFF", "Low Activity" = "#d19db4", "Moderate Activity" = "#8b0c43", "High Activity" = "#38041b")) +  # Customize colors for RA_status
  scale_x_log10() +
  labs(x = "Cytokine Value", y = "Density", title = "Density Plots of Cytokine Levels by CDAI") +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title= element_text(size=10),
    legend.text = element_text(size=8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )  # Additional theming options


gg <- ggplot(cytokine_data_CDAI_sigs, aes(x = CDAI_score, y = Value)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.7)) + 
  geom_smooth(
    method = "glm", 
    method.args = list(family = "poisson"), 
    colour = "black", 
    linewidth = 0.5, 
    alpha = 0.4
  ) +
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4, scales = "free_y") +  
  scale_color_manual(values = c(
    "High Activity" = "#380019", 
    "Moderate Activity" = "#af2460", 
    "Low Activity"= "#e89e9f",
    "Remission" = "#dcb6c7")) +
  xlim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI score", 
    x = "pg/mL", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: RA and Controls"
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title= element_text(size=10),
    legend.text = element_text(size=8),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 8)  # Size of facet labels
  )
# Print the plot
print(gg)
#2a0414
#601234
#5e002a












quartz()

####ORIGINAL: ATTEMPTED GRAPHS####
# install.packages("ggplot2")
# load package and data
options(scipen=999)  # turn-off scientific notation like 1e+48
library(ggplot2)
install.packages("ggExtra")
library(ggExtra)
theme_set(theme_minimal())  # pre-set the bw theme.
data("cytokine_data_long", package = "ggplot2")
# midwest <- read.csv("http://goo.gl/G1K41K")  # bkup data source

# Scatterplot
gg <- ggplot(cytokine_data_CDAI_sigs, aes(x=Value, y=CDAI_score)) + 
  geom_point(aes(col=Disease_Activity)) + 
  geom_smooth(method = "glm",
              method.args = list(family = poisson(link = "log")),
              formula = y ~ x,  # Poisson regression requires y ~ x by default
              se = TRUE,  # Display confidence interval
              color = "grey"))  +
  facet_wrap(~  Cytokine, ncol=8, nrow = 4) +  
  scale_x_log10() +
  ylim(c(0, 60)) + 
  labs(subtitle="Cytokines across CDAI scores", 
       y="CDAI", 
       x="pg/mL", 
       title="Scatterplot", 
       caption = "Source: RA and Controls") +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
    )
plot(gg)

cytokine_data_CDAI_sigs$Predictions <- NA

for (cytokine in cytokines) {
  # Build the formula
  formula_cyt2 <- as.formula(paste("CDAI_score ~ Diet_bin + Menopausal_bin_cat +", cytokine))
  
  # Fit the model
  model_cyt2 <- glm(formula_cyt2, family = poisson(link = "log"), data = cytokine_data_CDAI_sigs)
  
  # Generate predictions (for the given dataset)
  cytokine_data_CDAI_sigs$Predictions <- predict(model_cyt2, newdata = cytokine_data_CDAI_sigs, type = "response")}

#### Plot the results
gg <- ggplot(cytokine_data_CDAI_sigs, aes(x = Value, y = CDAI_score)) + 
  geom_point(aes(col = Disease_Activity, alpha = 0.5)) + 
  geom_line(aes(y = Predictions), color = "black", alpha = 0.4) +  # Plot predictions as lines
  facet_wrap(~ Cytokine, ncol = 8, nrow = 4) +  
  scale_x_log10() +
  scale_color_manual(values = c("High Activity" = "#460622", "Moderate Activity" = "#8b0c43", "Low Activity" = "#ae557b", "Remission" = "#dcb6c7")) +
  ylim(c(0, 60)) + 
  labs(
    subtitle = "Cytokines across CDAI scores", 
    y = "CDAI", 
    x = "log(pg/mL)", 
    title = "Cytokine Concentrations as Potential Modifiers of CDAI", 
    caption = "Source: RA and Controls"
  ) +
  theme(
    legend.position = "bottom",  # Position the legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for better readability
    strip.text.x = element_text(size = 10)  # Size of facet labels
  )

# Print the plot
print(gg)






















#### VAGINAL VARIABLES Mic AND RA [ACPA] ####
#Print p-val and coefficients of microbial impact on RA status
head(df_ACCPlsp[203:263])
microbes <- colnames(df_ACCPlsp)[203:263]
factors <- c("ACPA_ng_mL")
p_value_micACCPmatrix <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                            dimnames = list(factors, microbes))
coef_micACCPmatrix <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                         dimnames = list(factors, microbes))
p_value_ACCPmicmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                                dimnames = list(microbes, factors))
coef_ACCPmicmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                             dimnames = list(microbes, factors))
# Loop through factors and microbes
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(factor, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", microbe))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_mic2, data = df_ACCPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary1 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary1)) {
        p_value_micACCPmatrix[factor, microbe] <- coef_summary1[microbe, "Pr(>|t|)"]
        coef_micACCPmatrix[factor, microbe] <- coef_summary1[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(microbe, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", factor))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_mic2, data = df_ACCPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary2 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_ACCPmicmatrix[microbe, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_ACCPmicmatrix[microbe, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_micACCPmatrix)
print("Coefficient Matrix:")
print(coef_micACCPmatrix)
# Optionally save to CSV files
write.csv(p_value_micACCPmatrix, "Summary/p_value_micmatrix2_MICROBEefACCP.csv", row.names = TRUE)
write.csv(coef_micACCPmatrix, "Summary/coef_micmatrix2_MICROBEefACCP.csv", row.names = TRUE)
write.csv(p_value_ACCPmicmatrix, "Summary/p_value_micmatrix2_ACCPefMICROBE.csv", row.names = TRUE)
write.csv(coef_ACCPmicmatrix, "Summary/coef_micmatrix2_ACCPefMICROBE.csv", row.names = TRUE)





#RA interaction version
# Define matrices to store p-values and coefficients for interactions
interaction_names <- c("RA_statusControl", "RA_statusRheumatoid Arthritis")
p_value_micACCPmatrixRA <- matrix(NA, nrow = length(microbes), ncol = length(interaction_names),
                         dimnames = list(microbes, interaction_names))
coef_micACCPmatrixRA <- matrix(NA, nrow = length(microbes), ncol = length(interaction_names),
                      dimnames = list(microbes, interaction_names))

# Loop through RA_status:microbes
for (microbe in microbes) {
  # Build the formula for the model
  formula_mic <- as.formula(paste("ACPA_ng_mL ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +",
                                  paste0("RA_status:", microbe)))
  # Fit the model with error handling
  model_mic <- tryCatch(lm(formula_mic, data = df_ACCPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_mic)) {
    coef_summary <- summary(model_mic)$coefficients
    
    # Loop through interactions to extract p-values and coefficients
    for (interaction in interaction_names) {
      interaction_term <- paste0(interaction, ":", microbe)
      
      if (interaction_term %in% rownames(coef_summary)) {
        p_value_micACCPmatrixRA[microbe, interaction] <- coef_summary[interaction_term, "Pr(>|t|)"]
        coef_micACCPmatrixRA[microbe, interaction] <- coef_summary[interaction_term, "Estimate"]
      } else {
        message(paste("Interaction", interaction_term, "not found for microbe", microbe))
      }
    }
  } else {
    # Model did not fit
    message(paste("Model failed for microbe", microbe))
  }
}

# View results
print("P-Value Matrix:")
print(p_value_micACCPmatrixRA)
print("Coefficient Matrix:")
print(coef_micACCPmatrixRA)
# Optionally save to CSV files
write.csv(p_value_micACCPmatrixRA, "Summary/p_value_micmatrix2_MICROBEefACCP_RA.csv", row.names = TRUE)
write.csv(coef_micACCPmatrixRA, "Summary/coef_micmatrix2_MICROBEefACCP_RA.csv", row.names = TRUE)




# Loop through microbes for RA_status:ACPA
interaction_names <- c("RA_statusControl:ACPA_ng_mL", "RA_statusRheumatoid Arthritis:ACPA_ng_mL")
p_value_ACCPmicmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(microbes),
                                  dimnames = list(interaction_names, microbes))
coef_ACCPmicmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(microbes),
                               dimnames = list(interaction_names, microbes))

# Loop through microbes
for (microbe in microbes) {
  # Build the formula for the model
  formula_mic4 <- as.formula(paste(microbe, " ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat + RA_status:ACPA_ng_mL"))
  # Fit the model with error handling
  model_mic4 <- tryCatch(lm(formula_mic4, data = df_ACCPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_mic4)) {
    coef_summary4 <- summary(model_mic4)$coefficients
    
    # Check if the microbe exists in the coefficient names
    for (interaction in interaction_names) {
      if (interaction %in% rownames(coef_summary4)) {      
      p_value_ACCPmicmatrixRA[interaction_names, microbe] <- coef_summary4[interaction_names, "Pr(>|t|)"]
      coef_ACCPmicmatrixRA[interaction_names, microbe] <- coef_summary4[interaction_names, "Estimate"]
    } else {
      # Microbe not found in the model coefficients
      message(paste("Microbe", interaction_names, "not found in the model for factor", factor))
    }
  }} else {
    # Model did not fit
    message(paste("Model failed for microbe", microbe, "and interaction", interaction_names))
  }}


write.csv(p_value_ACCPmicmatrixRA, "Summary/p_value_micmatrix2_ACCPefMICROBE_RA.csv", row.names = TRUE)
write.csv(coef_ACCPmicmatrixRA, "Summary/coef_micmatrix2_ACCPefMICROBE_RA.csv", row.names = TRUE)


#### VAGINAL VARIABLES Mic AND RA [RF] ####
#Print p-val and coefficients of microbial impact on RA status
head(df_RFlsp[203:268])
microbes <- colnames(df_RFlsp)[203:268]
factors <- c("RF_U_mL")
p_value_micRFmatrix <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                                dimnames = list(factors, microbes))
coef_micRFmatrix <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                             dimnames = list(factors, microbes))
p_value_RFmicmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                                dimnames = list(microbes, factors))
coef_RFmicmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                             dimnames = list(microbes, factors))
# Loop through factors and microbes
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(factor, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", microbe))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_mic2, data = df_RFlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary1 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary1)) {
        p_value_micRFmatrix[factor, microbe] <- coef_summary1[microbe, "Pr(>|t|)"]
        coef_micRFmatrix[factor, microbe] <- coef_summary1[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(microbe, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", factor))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_mic2, data = df_RFlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary2 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_RFmicmatrix[microbe, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_RFmicmatrix[microbe, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_micRFmatrix)
print("Coefficient Matrix:")
print(coef_micRFmatrix)
# Optionally save to CSV files
write.csv(p_value_micRFmatrix, "Summary/p_value_micmatrix2_MICROBEefRF.csv", row.names = TRUE)
write.csv(coef_micRFmatrix, "Summary/coef_micmatrix2_MICROBEefRF.csv", row.names = TRUE)
write.csv(p_value_RFmicmatrix, "Summary/p_value_micmatrix2_RFefMICROBE.csv", row.names = TRUE)
write.csv(coef_RFmicmatrix, "Summary/coef_micmatrix2_RFefMICROBE.csv", row.names = TRUE)





#RA interaction version
# Define matrices to store p-values and coefficients for interactions
interaction_names <- c("RA_statusControl", "RA_statusRheumatoid Arthritis")
p_value_micRFmatrixRA <- matrix(NA, nrow = length(microbes), ncol = length(interaction_names),
                                  dimnames = list(microbes, interaction_names))
coef_micRFmatrixRA <- matrix(NA, nrow = length(microbes), ncol = length(interaction_names),
                               dimnames = list(microbes, interaction_names))

# Loop through RA_status:microbes
for (microbe in microbes) {
  # Build the formula for the model
  formula_mic <- as.formula(paste("RF_U_mL ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +",
                                  paste0("RA_status:", microbe)))
  # Fit the model with error handling
  model_mic <- tryCatch(lm(formula_mic, data = df_RFlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_mic)) {
    coef_summary <- summary(model_mic)$coefficients
    
    # Loop through interactions to extract p-values and coefficients
    for (interaction in interaction_names) {
      interaction_term <- paste0(interaction, ":", microbe)
      
      if (interaction_term %in% rownames(coef_summary)) {
        p_value_micRFmatrixRA[microbe, interaction] <- coef_summary[interaction_term, "Pr(>|t|)"]
        coef_micRFmatrixRA[microbe, interaction] <- coef_summary[interaction_term, "Estimate"]
      } else {
        message(paste("Interaction", interaction_term, "not found for microbe", microbe))
      }
    }
  } else {
    # Model did not fit
    message(paste("Model failed for microbe", microbe))
  }
}

# View results
print("P-Value Matrix:")
print(p_value_micRFmatrixRA)
print("Coefficient Matrix:")
print(coef_micRFmatrixRA)
# Optionally save to CSV files
write.csv(p_value_micRFmatrixRA, "Summary/p_value_micmatrix2_MICROBEefRF_RA.csv", row.names = TRUE)
write.csv(coef_micRFmatrixRA, "Summary/coef_micmatrix2_MICROBEefRF_RA.csv", row.names = TRUE)




# Loop through microbes for RA_status:ACPA
interaction_names <- c("RA_statusControl:RF_U_mL", "RA_statusRheumatoid Arthritis:RF_U_mL")
p_value_RFmicmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(microbes),
                                  dimnames = list(interaction_names, microbes))
coef_RFmicmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(microbes),
                               dimnames = list(interaction_names, microbes))

# Loop through microbes
for (microbe in microbes) {
  # Build the formula for the model
  formula_mic4 <- as.formula(paste(microbe, " ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat + RA_status:RF_U_mL"))
  # Fit the model with error handling
  model_mic4 <- tryCatch(lm(formula_mic4, data = df_RFlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_mic4)) {
    coef_summary4 <- summary(model_mic4)$coefficients
    
    # Check if the microbe exists in the coefficient names
    for (interaction in interaction_names) {
      if (interaction %in% rownames(coef_summary4)) {      
        p_value_RFmicmatrixRA[interaction_names, microbe] <- coef_summary4[interaction_names, "Pr(>|t|)"]
        coef_RFmicmatrixRA[interaction_names, microbe] <- coef_summary4[interaction_names, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", interaction_names, "not found in the model for factor", factor))
      }
    }} else {
      # Model did not fit
      message(paste("Model failed for microbe", microbe, "and interaction", interaction_names))
    }}


write.csv(p_value_RFmicmatrixRA, "Summary/p_value_micmatrix2_RFefMICROBE_RA.csv", row.names = TRUE)
write.csv(coef_RFmicmatrixRA, "Summary/coef_micmatrix2_RFefMICROBE_RA.csv", row.names = TRUE)

#### VAGINAL VARIABLES Mic AND RA [CRP] ####
#Print p-val and coefficients of microbial impact on RA status
head(df_CRPlsp[203:263])
microbes <- colnames(df_CRPlsp)[203:263]
factors <- c("C_reactive_pg_mL")
p_value_micCRPmatrix <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                                dimnames = list(factors, microbes))
coef_micCRPmatrix <- matrix(NA, nrow = length(factors), ncol = length(microbes),
                             dimnames = list(factors, microbes))
p_value_CRPmicmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                                dimnames = list(microbes, factors))
coef_CRPmicmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors),
                             dimnames = list(microbes, factors))
# Loop through factors and microbes
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(factor, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", microbe))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_mic2, data = df_CRPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary1 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (microbe %in% rownames(coef_summary1)) {
        p_value_micCRPmatrix[factor, microbe] <- coef_summary1[microbe, "Pr(>|t|)"]
        coef_micCRPmatrix[factor, microbe] <- coef_summary1[microbe, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(microbe, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", factor))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_mic2, data = df_CRPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_mic2)) {
      coef_summary2 <- summary(model_mic2)$coefficients
      # Check if the microbe exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_CRPmicmatrix[microbe, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_CRPmicmatrix[microbe, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and microbe", microbe))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_micCRPmatrix)
print("Coefficient Matrix:")
print(coef_micCRPmatrix)
# Optionally save to CSV files
write.csv(p_value_micCRPmatrix, "Summary/p_value_micmatrix2_MICROBEefCRP.csv", row.names = TRUE)
write.csv(coef_micCRPmatrix, "Summary/coef_micmatrix2_MICROBEefCRP.csv", row.names = TRUE)
write.csv(p_value_CRPmicmatrix, "Summary/p_value_micmatrix2_CRPefMICROBE.csv", row.names = TRUE)
write.csv(coef_CRPmicmatrix, "Summary/coef_micmatrix2_CRPefMICROBE.csv", row.names = TRUE)





#RA interaction version
# Define matrices to store p-values and coefficients for interactions
interaction_names <- c("RA_statusControl", "RA_statusRheumatoid Arthritis")
p_value_micCRPmatrixRA <- matrix(NA, nrow = length(microbes), ncol = length(interaction_names),
                                  dimnames = list(microbes, interaction_names))
coef_micCRPmatrixRA <- matrix(NA, nrow = length(microbes), ncol = length(interaction_names),
                               dimnames = list(microbes, interaction_names))

# Loop through RA_status:microbes
for (microbe in microbes) {
  # Build the formula for the model
  formula_mic <- as.formula(paste("C_reactive_pg_mL ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +",
                                  paste0("RA_status:", microbe)))
  # Fit the model with error handling
  model_mic <- tryCatch(lm(formula_mic, data = df_CRPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_mic)) {
    coef_summary <- summary(model_mic)$coefficients
    
    # Loop through interactions to extract p-values and coefficients
    for (interaction in interaction_names) {
      interaction_term <- paste0(interaction, ":", microbe)
      
      if (interaction_term %in% rownames(coef_summary)) {
        p_value_micCRPmatrixRA[microbe, interaction] <- coef_summary[interaction_term, "Pr(>|t|)"]
        coef_micCRPmatrixRA[microbe, interaction] <- coef_summary[interaction_term, "Estimate"]
      } else {
        message(paste("Interaction", interaction_term, "not found for microbe", microbe))
      }
    }
  } else {
    # Model did not fit
    message(paste("Model failed for microbe", microbe))
  }
}

# View results
print("P-Value Matrix:")
print(p_value_micCRPmatrixRA)
print("Coefficient Matrix:")
print(coef_micCRPmatrixRA)
# Optionally save to CSV files
write.csv(p_value_micCRPmatrixRA, "Summary/p_value_micmatrix2_MICROBEefCRP_RA.csv", row.names = TRUE)
write.csv(coef_micCRPmatrixRA, "Summary/coef_micmatrix2_MICROBEefCRP_RA.csv", row.names = TRUE)




# Loop through microbes for RA_status:ACPA
interaction_names <- c("RA_statusControl:C_reactive_pg_mL", "RA_statusRheumatoid Arthritis:C_reactive_pg_mL")
p_value_CRPmicmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(microbes),
                                  dimnames = list(interaction_names, microbes))
coef_CRPmicmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(microbes),
                               dimnames = list(interaction_names, microbes))

# Loop through microbes
for (microbe in microbes) {
  # Build the formula for the model
  formula_mic4 <- as.formula(paste(microbe, " ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat + RA_status:C_reactive_pg_mL"))
  # Fit the model with error handling
  model_mic4 <- tryCatch(lm(formula_mic4, data = df_CRPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_mic4)) {
    coef_summary4 <- summary(model_mic4)$coefficients
    
    # Check if the microbe exists in the coefficient names
    for (interaction in interaction_names) {
      if (interaction %in% rownames(coef_summary4)) {      
        p_value_CRPmicmatrixRA[interaction_names, microbe] <- coef_summary4[interaction_names, "Pr(>|t|)"]
        coef_CRPmicmatrixRA[interaction_names, microbe] <- coef_summary4[interaction_names, "Estimate"]
      } else {
        # Microbe not found in the model coefficients
        message(paste("Microbe", interaction_names, "not found in the model for factor", factor))
      }
    }} else {
      # Model did not fit
      message(paste("Model failed for microbe", microbe, "and interaction", interaction_names))
    }}


write.csv(p_value_CRPmicmatrixRA, "Summary/p_value_micmatrix2_CRPefMICROBE_RA.csv", row.names = TRUE)
write.csv(coef_CRPmicmatrixRA, "Summary/coef_micmatrix2_CRPefMICROBE_RA.csv", row.names = TRUE)




#### VAGINAL VARIABLES AND RA [ACPA] ####
#Print p-val and coefficients of cytrobial impact on RA status
head(df_ACCPlsp[159:197])
cytokines <- colnames(df_ACCPlsp)[159:197]
factors <- c("ACPA_ng_mL")
p_value_cytACCPmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                                dimnames = list(factors, cytokines))
coef_cytACCPmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                             dimnames = list(factors, cytokines))
stderror_cytACCPmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                                dimnames = list(factors, cytokines))
tvalue_cytACCPmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                             dimnames = list(factors, cytokines))
p_value_ACCPcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                                dimnames = list(cytokines, factors))
stderror_ACCPcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                             dimnames = list(cytokines, factors))
tvalue_ACCPcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                                dimnames = list(cytokines, factors))
coef_ACCPcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                             dimnames = list(cytokines, factors))
# Loop through factors and cytokines
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = df_ACCPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary1 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary1)) {
        p_value_cytACCPmatrix[factor, cytokine] <- coef_summary1[cytokine, "Pr(>|t|)"]
        coef_cytACCPmatrix[factor, cytokine] <- coef_summary1[cytokine, "Estimate"]
        tvalue_cytACCPmatrix[factor, cytokine] <- coef_summary1[cytokine, "t value"]
        stderror_cytACCPmatrix[factor, cytokine] <- coef_summary1[cytokine, "Std. Error"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = df_ACCPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_ACCPcytmatrix[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_ACCPcytmatrix[cytokine, factor] <- coef_summary2[factor, "Estimate"]
        tvalue_ACCPcytmatrix[cytokine, factor] <- coef_summary2[factor, "t value"]
        stderror_ACCPcytmatrix[cytokine, factor] <- coef_summary2[factor, "Std. Error"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_cytACCPmatrix)
print("Coefficient Matrix:")
print(coef_cytACCPmatrix)
# Optionally save to CSV files
write.csv(p_value_cytACCPmatrix, "Summary/p_value_cytmatrix2_cytokineefACCP.csv", row.names = TRUE)
write.csv(coef_cytACCPmatrix, "Summary/coef_cytmatrix2_cytokineefACCP.csv", row.names = TRUE)
write.csv(p_value_ACCPcytmatrix, "Summary/p_value_cytmatrix2_ACCPefcytokine.csv", row.names = TRUE)
write.csv(coef_ACCPcytmatrix, "Summary/coef_cytmatrix2_ACCPefcytokine.csv", row.names = TRUE)
write.csv(stderror_cytACCPmatrix, "Summary/stderror_cytmatrix2_cytokineefACCP.csv", row.names = TRUE)
write.csv(tvalue_cytACCPmatrix, "Summary/tvalue_cytmatrix2_cytokineefACCP.csv", row.names = TRUE)
write.csv(stderror_ACCPcytmatrix, "Summary/stderror_cytmatrix2_ACCPefcytokine.csv", row.names = TRUE)
write.csv(tvalue_ACCPcytmatrix, "Summary/tvalue_cytmatrix2_ACCPefcytokine.csv", row.names = TRUE)




#RA interaction version
# Define matrices to store p-values and coefficients for interactions
interaction_names <- c("RA_statusControl", "RA_statusRheumatoid Arthritis")
p_value_cytACCPmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                                  dimnames = list(cytokines, interaction_names))
coef_cytACCPmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                               dimnames = list(cytokines, interaction_names))
tvalue_cytACCPmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                                  dimnames = list(cytokines, interaction_names))
stderror_cytACCPmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                               dimnames = list(cytokines, interaction_names))

# Loop through RA_status:cytokines
for (cytokine in cytokines) {
  # Build the formula for the model
  formula_cyt <- as.formula(paste("ACPA_ng_mL ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +",
                                  paste0("RA_status:", cytokine)))
  # Fit the model with error handling
  model_cyt <- tryCatch(lm(formula_cyt, data = df_ACCPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_cyt)) {
    coef_summary <- summary(model_cyt)$coefficients
    
    # Loop through interactions to extract p-values and coefficients
    for (interaction in interaction_names) {
      interaction_term <- paste0(interaction, ":", cytokine)
      
      if (interaction_term %in% rownames(coef_summary)) {
        p_value_cytACCPmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Pr(>|t|)"]
        coef_cytACCPmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Estimate"]
        tvalue_cytACCPmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "t value"]
        stderror_cytACCPmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Std. Error"]
      } else {
        message(paste("Interaction", interaction_term, "not found for cytokine", cytokine))
      }
    }
  } else {
    # Model did not fit
    message(paste("Model failed for cytokine", cytokine))
  }
}

# View results
print("P-Value Matrix:")
print(p_value_cytACCPmatrixRA)
print("Coefficient Matrix:")
print(coef_cytACCPmatrixRA)
# Optionally save to CSV files
write.csv(p_value_cytACCPmatrixRA, "Summary/p_value_cytmatrix2_cytokineefACCP_RA.csv", row.names = TRUE)
write.csv(coef_cytACCPmatrixRA, "Summary/coef_cytmatrix2_cytokineefACCP_RA.csv", row.names = TRUE)
write.csv(tvalue_cytACCPmatrixRA, "Summary/tvalue_cytmatrix2_cytokineefACCP_RA.csv", row.names = TRUE)
write.csv(stderror_cytACCPmatrixRA, "Summary/stderror_cytmatrix2_cytokineefACCP_RA.csv", row.names = TRUE)



# Loop through cytokines for RA_status:ACPA
interaction_names <- c("RA_statusControl:ACPA_ng_mL", "RA_statusRheumatoid Arthritis:ACPA_ng_mL")
p_value_ACCPcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                                  dimnames = list(interaction_names, cytokines))
coef_ACCPcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                               dimnames = list(interaction_names, cytokines))
tvalue_ACCPcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                                  dimnames = list(interaction_names, cytokines))
stderror_ACCPcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                               dimnames = list(interaction_names, cytokines))

# Loop through cytokines
for (cytokine in cytokines) {
  # Build the formula for the model
  formula_cyt4 <- as.formula(paste(cytokine, " ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat + RA_status:ACPA_ng_mL"))
  # Fit the model with error handling
  model_cyt4 <- tryCatch(lm(formula_cyt4, data = df_ACCPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_cyt4)) {
    coef_summary4 <- summary(model_cyt4)$coefficients
    
    # Check if the cytokine exists in the coefficient names
    for (interaction in interaction_names) {
      if (interaction %in% rownames(coef_summary4)) {      
        p_value_ACCPcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Pr(>|t|)"]
        coef_ACCPcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Estimate"]
        tvalue_ACCPcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "t value"]
        stderror_ACCPcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Std. Error"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", interaction_names, "not found in the model for factor", factor))
      }
    }} else {
      # Model did not fit
      message(paste("Model failed for cytokine", cytokine, "and interaction", interaction_names))
    }}

write.csv(p_value_ACCPcytmatrixRA, "Summary/p_value_cytmatrix2_ACCPefcytokine_RA.csv", row.names = TRUE)
write.csv(coef_ACCPcytmatrixRA, "Summary/coef_cytmatrix2_ACCPefcytokine_RA.csv", row.names = TRUE)
write.csv(tvalue_ACCPcytmatrixRA, "Summary/tvalue_cytmatrix2_ACCPefcytokine_RA.csv", row.names = TRUE)
write.csv(stderror_ACCPcytmatrixRA, "Summary/stderror_cytmatrix2_ACCPefcytokine_RA.csv", row.names = TRUE)



#### VAGINAL VARIABLES AND RA [RF] ####
#Print p-val and coefficients of cytrobial impact on RA status
head(df_RFlsp[159:197])
cytokines <- colnames(df_RFlsp)[159:197]
factors <- c("RF_U_mL")
p_value_cytRFmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                              dimnames = list(factors, cytokines))
coef_cytRFmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                           dimnames = list(factors, cytokines))
p_value_RFcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                              dimnames = list(cytokines, factors))
coef_RFcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                           dimnames = list(cytokines, factors))
# Loop through factors and cytokines
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = df_RFlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary1 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary1)) {
        p_value_cytRFmatrix[factor, cytokine] <- coef_summary1[cytokine, "Pr(>|t|)"]
        coef_cytRFmatrix[factor, cytokine] <- coef_summary1[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = df_RFlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_RFcytmatrix[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_RFcytmatrix[cytokine, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_cytRFmatrix)
print("Coefficient Matrix:")
print(coef_cytRFmatrix)
# Optionally save to CSV files
write.csv(p_value_cytRFmatrix, "Summary/p_value_cytmatrix2_cytokineefRF.csv", row.names = TRUE)
write.csv(coef_cytRFmatrix, "Summary/coef_cytmatrix2_cytokineefRF.csv", row.names = TRUE)
write.csv(p_value_RFcytmatrix, "Summary/p_value_cytmatrix2_RFefcytokine.csv", row.names = TRUE)
write.csv(coef_RFcytmatrix, "Summary/coef_cytmatrix2_RFefcytokine.csv", row.names = TRUE)





#RA interaction version
# Define matrices to store p-values and coefficients for interactions
interaction_names <- c("RA_statusControl", "RA_statusRheumatoid Arthritis")
p_value_cytRFmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                                dimnames = list(cytokines, interaction_names))
coef_cytRFmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                             dimnames = list(cytokines, interaction_names))

# Loop through RA_status:cytokines
for (cytokine in cytokines) {
  # Build the formula for the model
  formula_cyt <- as.formula(paste("RF_U_mL ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +",
                                  paste0("RA_status:", cytokine)))
  # Fit the model with error handling
  model_cyt <- tryCatch(lm(formula_cyt, data = df_RFlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_cyt)) {
    coef_summary <- summary(model_cyt)$coefficients
    
    # Loop through interactions to extract p-values and coefficients
    for (interaction in interaction_names) {
      interaction_term <- paste0(interaction, ":", cytokine)
      
      if (interaction_term %in% rownames(coef_summary)) {
        p_value_cytRFmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Pr(>|t|)"]
        coef_cytRFmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Estimate"]
      } else {
        message(paste("Interaction", interaction_term, "not found for cytokine", cytokine))
      }
    }
  } else {
    # Model did not fit
    message(paste("Model failed for cytokine", cytokine))
  }
}

# View results
print("P-Value Matrix:")
print(p_value_cytRFmatrixRA)
print("Coefficient Matrix:")
print(coef_cytRFmatrixRA)
# Optionally save to CSV files
write.csv(p_value_cytRFmatrixRA, "Summary/p_value_cytmatrix2_cytokineefRF_RA.csv", row.names = TRUE)
write.csv(coef_cytRFmatrixRA, "Summary/coef_cytmatrix2_cytokineefRF_RA.csv", row.names = TRUE)




# Loop through cytokines for RA_status:ACPA
interaction_names <- c("RA_statusControl:RF_U_mL", "RA_statusRheumatoid Arthritis:RF_U_mL")
p_value_RFcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                                dimnames = list(interaction_names, cytokines))
coef_RFcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                             dimnames = list(interaction_names, cytokines))

# Loop through cytokines
for (cytokine in cytokines) {
  # Build the formula for the model
  formula_cyt4 <- as.formula(paste(cytokine, " ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat + RA_status:RF_U_mL"))
  # Fit the model with error handling
  model_cyt4 <- tryCatch(lm(formula_cyt4, data = df_RFlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_cyt4)) {
    coef_summary4 <- summary(model_cyt4)$coefficients
    
    # Check if the cytokine exists in the coefficient names
    for (interaction in interaction_names) {
      if (interaction %in% rownames(coef_summary4)) {      
        p_value_RFcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Pr(>|t|)"]
        coef_RFcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", interaction_names, "not found in the model for factor", factor))
      }
    }} else {
      # Model did not fit
      message(paste("Model failed for cytokine", cytokine, "and interaction", interaction_names))
    }}


write.csv(p_value_RFcytmatrixRA, "Summary/p_value_cytmatrix2_RFefcytokine_RA.csv", row.names = TRUE)
write.csv(coef_RFcytmatrixRA, "Summary/coef_cytmatrix2_RFefcytokine_RA.csv", row.names = TRUE)

#### VAGINAL VARIABLES AND RA [CRP] ####
#Print p-val and coefficients of cytrobial impact on RA status
head(df_CRPlsp[159:197])
cytokines <- colnames(df_CRPlsp)[159:197]
factors <- c("C_reactive_pg_mL")
p_value_cytCRPmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                               dimnames = list(factors, cytokines))
coef_cytCRPmatrix <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                            dimnames = list(factors, cytokines))
p_value_CRPcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                               dimnames = list(cytokines, factors))
coef_CRPcytmatrix <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                            dimnames = list(cytokines, factors))
# Loop through factors and cytokines
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = df_CRPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary1 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary1)) {
        p_value_cytCRPmatrix[factor, cytokine] <- coef_summary1[cytokine, "Pr(>|t|)"]
        coef_cytCRPmatrix[factor, cytokine] <- coef_summary1[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = df_CRPlsp), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_CRPcytmatrix[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_CRPcytmatrix[cytokine, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
# View results
print("P-Value Matrix:")
print(p_value_cytCRPmatrix)
print("Coefficient Matrix:")
print(coef_cytCRPmatrix)
# Optionally save to CSV files
write.csv(p_value_cytCRPmatrix, "Summary/p_value_cytmatrix2_cytokineefCRP.csv", row.names = TRUE)
write.csv(coef_cytCRPmatrix, "Summary/coef_cytmatrix2_cytokineefCRP.csv", row.names = TRUE)
write.csv(p_value_CRPcytmatrix, "Summary/p_value_cytmatrix2_CRPefcytokine.csv", row.names = TRUE)
write.csv(coef_CRPcytmatrix, "Summary/coef_cytmatrix2_CRPefcytokine.csv", row.names = TRUE)





#RA interaction version
# Define matrices to store p-values and coefficients for interactions
interaction_names <- c("RA_statusControl", "RA_statusRheumatoid Arthritis")
p_value_cytCRPmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                                 dimnames = list(cytokines, interaction_names))
coef_cytCRPmatrixRA <- matrix(NA, nrow = length(cytokines), ncol = length(interaction_names),
                              dimnames = list(cytokines, interaction_names))

# Loop through RA_status:cytokines
for (cytokine in cytokines) {
  # Build the formula for the model
  formula_cyt <- as.formula(paste("C_reactive_pg_mL ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat +",
                                  paste0("RA_status:", cytokine)))
  # Fit the model with error handling
  model_cyt <- tryCatch(lm(formula_cyt, data = df_CRPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_cyt)) {
    coef_summary <- summary(model_cyt)$coefficients
    
    # Loop through interactions to extract p-values and coefficients
    for (interaction in interaction_names) {
      interaction_term <- paste0(interaction, ":", cytokine)
      
      if (interaction_term %in% rownames(coef_summary)) {
        p_value_cytCRPmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Pr(>|t|)"]
        coef_cytCRPmatrixRA[cytokine, interaction] <- coef_summary[interaction_term, "Estimate"]
      } else {
        message(paste("Interaction", interaction_term, "not found for cytokine", cytokine))
      }
    }
  } else {
    # Model did not fit
    message(paste("Model failed for cytokine", cytokine))
  }
}

# View results
print("P-Value Matrix:")
print(p_value_cytCRPmatrixRA)
print("Coefficient Matrix:")
print(coef_cytCRPmatrixRA)
# Optionally save to CSV files
write.csv(p_value_cytCRPmatrixRA, "Summary/p_value_cytmatrix2_cytokineefCRP_RA.csv", row.names = TRUE)
write.csv(coef_cytCRPmatrixRA, "Summary/coef_cytmatrix2_cytokineefCRP_RA.csv", row.names = TRUE)




# Loop through cytokines for RA_status:ACPA
interaction_names <- c("RA_statusControl:C_reactive_pg_mL", "RA_statusRheumatoid Arthritis:C_reactive_pg_mL")
p_value_CRPcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                                 dimnames = list(interaction_names, cytokines))
coef_CRPcytmatrixRA <- matrix(NA, nrow = length(interaction_names), ncol = length(cytokines),
                              dimnames = list(interaction_names, cytokines))

# Loop through cytokines
for (cytokine in cytokines) {
  # Build the formula for the model
  formula_cyt4 <- as.formula(paste(cytokine, " ~ Diet_bin + Age + Lifetime_Sexual_Partners_bin_cat + RA_status:C_reactive_pg_mL"))
  # Fit the model with error handling
  model_cyt4 <- tryCatch(lm(formula_cyt4, data = df_CRPlsp), error = function(e) NULL)
  
  # Check if the model fit successfully
  if (!is.null(model_cyt4)) {
    coef_summary4 <- summary(model_cyt4)$coefficients
    
    # Check if the cytokine exists in the coefficient names
    for (interaction in interaction_names) {
      if (interaction %in% rownames(coef_summary4)) {      
        p_value_CRPcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Pr(>|t|)"]
        coef_CRPcytmatrixRA[interaction_names, cytokine] <- coef_summary4[interaction_names, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", interaction_names, "not found in the model for factor", factor))
      }
    }} else {
      # Model did not fit
      message(paste("Model failed for cytokine", cytokine, "and interaction", interaction_names))
    }}


write.csv(p_value_CRPcytmatrixRA, "Summary/p_value_cytmatrix2_CRPefcytokine_RA.csv", row.names = TRUE)
write.csv(coef_CRPcytmatrixRA, "Summary/coef_cytmatrix2_CRPefcytokine_RA.csv", row.names = TRUE)












cytokines <- colnames(your_data_4factorsESR)[158:205]
factors <- c("Erythrocyte_sedementation_rate_.mm.hr.")
p_value_cytmatrix_RAx <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                                dimnames = list(cytokines, factors))
coef_cytmatrix_RAx <- matrix(NA, nrow = length(cytokines), ncol = length(factors),
                             dimnames = list(cytokines, factors))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(cytokine, "~ Diet_bin + Menopausal_bin_cat +", factor))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = your_data_4factorsESR), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (factor %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RAx[cytokine, factor] <- coef_summary2[factor, "Pr(>|t|)"]
        coef_cytmatrix_RAx[cytokine, factor] <- coef_summary2[factor, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RAx, "Summary/p_value_cytRAmatrix2_ESRefcyt.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RAx, "Summary/coef_cytRAmatrix2_ESRefcyt.csv", row.names = TRUE)

p_value_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                               dimnames = list(factors, cytokines))
coef_cytmatrix_RA <- matrix(NA, nrow = length(factors), ncol = length(cytokines),
                            dimnames = list(factors, cytokines))
for (factor in factors) {
  for (cytokine in cytokines) {
    # Build the formula
    formula_cyt2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", cytokine))
    # Fit the model with error handling
    model_cyt2 <- tryCatch(lm(formula_cyt2, data = your_data_4factorsESR), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (cytokine %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Pr(>|t|)"]
        coef_cytmatrix_RA[factor, cytokine] <- coef_summary2[cytokine, "Estimate"]
      } else {
        # cytokine not found in the model coefficients
        message(paste("cytokine", cytokine, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", cytokine))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_cytefESR.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_cytefESR.csv", row.names = TRUE)




head(your_data_4factorsESR[211:215])
microbes <- colnames(your_data_4factors)[211:215] #*include to 286 for lm)
factors1 <- c("Erythrocyte_sedementation_rate_.mm.hr.")
p_value_micRAmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors1),
                              dimnames = list(microbes, factors1))
coef_micRAmatrix <- matrix(NA, nrow = length(microbes), ncol = length(factors1),
                           dimnames = list(microbes, factors1))
for (factor in factors) {
  for (microbe in microbes) {
    # Build the formula
    formula_mic2 <- as.formula(paste(factor, "~ Diet_bin + Menopausal_bin_cat +", microbe))
    # Fit the model with error handling
    model_mic2 <- tryCatch(lm(formula_cyt2, data = your_data_4factorsESR), 
                           error = function(e) NULL)
    # Check if the model fit successfully
    if (!is.null(model_cyt2)) {
      coef_summary2 <- summary(model_cyt2)$coefficients
      # Check if the cytokine exists in the coefficient names
      if (microbe %in% rownames(coef_summary2)) {
        p_value_cytmatrix_RA[factor, microbe] <- coef_summary2[microbe, "Pr(>|t|)"]
        coef_cytmatrix_RA[factor, microbe] <- coef_summary2[microbe, "Estimate"]
      } else {
        # microbe not found in the model coefficients
        message(paste("microbe", microbe, "not found in the model for factor", factor))
      }
    } else {
      # Model did not fit
      message(paste("Model failed for factor", factor, "and cytokine", microbe))
    }}}
write.csv(p_value_cytmatrix_RA, "Summary/p_value_cytRAmatrix2_microbeefESR.csv", row.names = TRUE)
write.csv(coef_cytmatrix_RA, "Summary/coef_cytRAmatrix2_microbeefESR.csv", row.names = TRUE)


your_data_model <- read.csv("Metadata_RAstudy_master_categorized_v2excess_FINALMODEL.csv")

modelRA_microbes <- glm(RA_status_cat ~ Total_Prevotella +
                        ACPA_detectable + Diet_bin +
                        Menopausal_bin_ordered	+	Clinic_cat + Abx_Cat,
                      family = binomial(link = "logit"), data = your_data_model)

summary(modelRA_microbes)
plot(modelRA_microbes)
vif(modelRA_microbes)
# Predict probabilities
pred_probsmicrobes <- predict(modelRA_microbes, type = "response")
length(your_data_model$RA_status_cat)
pred_classmicrobes <- ifelse(pred_probsmicrobes > 0.5, 1, 0)
length(pred_classmicrobes)
table(Predicted = pred_classmicrobes, Actual = your_data_model$RA_status_cat)
roc_curvemicrobes <- roc(your_data_model$RA_status_cat, pred_probsmicrobes)
plot(roc_curvemicrobes)
auc(roc_curvemicrobes)




your_data_immune_fin$Menopausal_bin <- factor(your_data_immune_fin$Menopausal_bin, 
                                              levels = c("Pre-menopause", "Peri/menopause", "Post-menopause"))
your_data_immune_fin$RA_status <- factor(your_data_immune_fin$RA_status, 
                                         levels = c("Control", "Rheumatoid Arthritis"))
your_data_immune_fin$Diet_bin <- factor(your_data_immune_fin$Diet_bin, 
                                        levels = c("Reduced", "High"))
your_data_immune_fin$Clinic <- factor(your_data_immune_fin$Clinic, 
                                      levels = c("McNair", "Smith"))
your_data_immune_fin$Lifetime_Sexual_Partners_bin <- factor(your_data_immune_fin$Lifetime_Sexual_Partners_bin, 
                                                            levels = c("0", "1 or 2", "3+"))

#various models
modelcyt <- lm(TGFa ~ Lifetime_Sexual_Partners_bin + Clinic + Menopausal_bin + Diet_bin,
               data = your_data_immune_fin)
summary(modelcyt)
modelcyt <- lm(IL.1a ~ RA_status + Lifetime_Sexual_Partners_bin + Clinic + Menopausal_bin + Diet_bin,
               data = your_data_immune_fin)
summary(modelcyt)

#clinic, menopausal, lsp
your_data_immune_fin$Diet_bin <- relevel(your_data_immune_fin$Diet_bin, ref = "Reduced")
your_data_immune_fin$Menopausal_bin <- relevel(your_data_immune_fin$Menopausal_bin, ref = "Pre-menopause")
modelcyt <- lm(IL.1a ~ RA_status + Menopausal_bin + Clinic + Diet_bin + RA_status:Menopausal_bin + RA_status:Clinic + RA_status:Diet_bin + Lifetime_Sexual_Partners_bin,
               data = your_data_immune_fin)
summary(modelcyt)
your_data_immune_fin$Menopausal_bin <- relevel(your_data_immune_fin$Menopausal_bin, ref = "Peri/menopause")
modelcyt <- lm(TGFa ~ RA_status + Menopausal_bin + Clinic + Diet_bin + RA_status:Menopausal_bin + RA_status:Diet_bin + Lifetime_Sexual_Partners_bin,
               data = your_data_immune_fin)
summary(modelcyt)
your_data_immune_fin$Menopausal_bin <- relevel(your_data_immune_fin$Menopausal_bin, ref = "Post-menopause")
modelcyt <- lm(TGFa ~ RA_status + Menopausal_bin + Clinic + Diet_bin + RA_status:Menopausal_bin + RA_status:Diet_bin + Lifetime_Sexual_Partners_bin,
               data = your_data_immune_fin)
summary(modelcyt)

modelcyt <- lm(IL.12.p40. ~ RA_status + Menopausal_bin_cat + Clinic_cat + RA_status:Menopausal_bin_cat + RA_status:Clinic_cat + Lifetime_Sexual_Partners_bin_cat + Diet_bin,
               data = your_data_immune_fin)
summary(modelcyt)

###### Define a function to calculate AUC #START HERE####
your_data_model <- read.csv("Metadata_RAstudy_master_categorized_v2excess_FINALMODEL.csv")
# Function to compute AUC for a given dataset
calculate_auc <- function(data, indices) {
  # Resample the data using the indices
  bootstrap_sample <- data[indices, ]
  # Fit the logistic regression model on the bootstrap sample
  model_bootstrap <- glm(RA_status_cat ~ Total_Prevotella + Total_Peptoniphilus +
                           EGF + sCD40L * Diet_bin +
                           Menopausal_bin_ordered	+	Clinic_cat + Abx_Cat,
                         family = binomial(link = "logit"), data = your_data_model)
  # Compute predicted probabilities
  pred_probs <- predict(model_bootstrap, type = "response", newdata = bootstrap_sample)
  # Calculate AUC
  auc_value <- roc(bootstrap_sample$RA_status, pred_probs)$auc
  return(auc_value)
}

####### Perform Bootsrap resampling ########
# Number of bootstrap iterations
n_iterations <- 1000
# Set up the bootstrap resampling
set.seed(123)  # For reproducibility
# Perform the bootstrap resampling
boot_results1 <- boot(data = your_data_model, statistic = calculate_auc, R = n_iterations)
# View the bootstrap results
boot_results1

model_bootstrap <- glm(RA_status_cat ~ Total_Prevotella + Total_Peptoniphilus +
                          EGF + sCD40L * Diet_bin +
                          Menopausal_bin_ordered	+	Clinic_cat + Abx_Cat,
                        family = binomial(link = "logit"), data = your_data_model)

summary(model_bootstrap)
plot(model_bootstrap)
vif(model_bootstrap)
# Predict probabilities
pred_probsmicrobes <- predict(model_bootstrap, type = "response")
length(your_data_model$RA_status_cat)
pred_classmicrobes <- ifelse(pred_probsmicrobes > 0.5, 1, 0)
length(pred_classmicrobes)
table(Predicted = pred_classmicrobes, Actual = your_data_model$RA_status_cat)
roc_curvemicrobes <- roc(your_data_model$RA_status_cat, pred_probsmicrobes)
plot(roc_curvemicrobes)
auc(roc_curvemicrobes)

####### Summarize ########

# Mean AUC from bootstrap resampling
mean_auc <- mean(boot_results1$t)
# Standard deviation of the AUC values
sd_auc <- sd(boot_results1$t)
# Confidence intervals for the AUC (e.g., 95% CI)
conf_interval1 <- boot.ci(boot_results1, type = "perc")
# Output results
cat("Mean AUC: ", mean_auc, "\n")
cat("Standard Deviation of AUC: ", sd_auc, "\n")
cat("95% Confidence Interval for AUC: ", conf_interval1$percent[4:5], "\n")
# Plot AUC distribution
hist(boot_results1$t, main = "Bootstrap AUC Distribution", 
     xlab = "AUC", col = "#8b0c43", border = "black")



#No ACPA
your_data_model_subsetcontrols2 <- read.csv("Metadata_RAstudy_master_categorized_subset.csv")

modelRA_microbes <- glm(RA_status_cat ~ Total_Prevotella + Total_Peptoniphilus +
                          EGF + sCD40L * Diet_bin +
                          Menopausal_bin_ordered	+	Clinic_cat + Abx_Cat,
                        family = binomial(link = "logit"), data = your_data_model_subsetcontrols2)
summary(modelRA_microbes)
plot(modelRA_microbes)
vif(modelRA_microbes)
# Predict probabilities
pred_probsmicrobes <- predict(modelRA_microbes, type = "response")
length(your_data_model_subsetcontrols2$RA_status_cat)
pred_classmicrobes <- ifelse(pred_probsmicrobes > 0.5, 1, 0)
length(pred_classmicrobes)
table(Predicted = pred_classmicrobes, Actual = your_data_model_subsetcontrols2$RA_status_cat)
roc_curvemicrobes <- roc(your_data_model_subsetcontrols2$RA_status_cat, pred_probsmicrobes)
plot(roc_curvemicrobes)
auc(roc_curvemicrobes)

###### Define a function to calculate AUC #####
# Function to compute AUC for a given dataset
calculate_auc <- function(data, indices) {
  # Resample the data using the indices
  bootstrap_sample <- data[indices, ]
  # Fit the logistic regression model on the bootstrap sample
  model_bootstrap <- glm(RA_status_cat ~ Total_Prevotella + Total_Peptoniphilus +
                           EGF + sCD40L * Diet_bin +
                           Menopausal_bin_ordered	+	Clinic_cat + Abx_Cat,
                         family = binomial(link = "logit"), data = your_data_model_subsetcontrols2)
  # Compute predicted probabilities
  pred_probs <- predict(model_bootstrap, type = "response", newdata = bootstrap_sample)
  # Calculate AUC
  auc_value <- roc(bootstrap_sample$RA_status, pred_probs)$auc
  return(auc_value)
}

####### Perform Bootsrap resampling ########
# Number of bootstrap iterations
n_iterations <- 1000
# Set up the bootstrap resampling
set.seed(123)  # For reproducibility
# Perform the bootstrap resampling
boot_results1 <- boot(data = your_data_model_subsetcontrols2, statistic = calculate_auc, R = n_iterations)
# View the bootstrap results
boot_results1

####### Summarize ########

# Mean AUC from bootstrap resampling
mean_auc <- mean(boot_results1$t)
# Standard deviation of the AUC values
sd_auc <- sd(boot_results1$t)
# Confidence intervals for the AUC (e.g., 95% CI)
conf_interval1 <- boot.ci(boot_results1, type = "perc")
# Output results
cat("Mean AUC: ", mean_auc, "\n")
cat("Standard Deviation of AUC: ", sd_auc, "\n")
cat("95% Confidence Interval for AUC: ", conf_interval1$percent[4:5], "\n")
# Plot AUC distribution
hist(boot_results1$t, main = "Bootstrap AUC Distribution", 
     xlab = "AUC", col = "#8b0c43", border = "black")




modelcyt <- lm(IL.15 ~ Diet_bin + Clinic_cat + Menopausal_bin_ordered + Methotrexate_dose_.mg.,
               data = df_meds)
summary(modelcyt)

