library(survival)       
library(ggsurvfit)
library(ggbreak)
library(survminer)

KTR_cohort_characteristics # phenotype data

# Classifying OTA exposure group
KTR_cohort_characteristics$OTA_level <- ifelse(KTR_cohort_characteristics$OTA < 0.14,"Low", "High")
KTR_cohort_characteristics$OTA_level <- ifelse(KTR_cohort_characteristics$OTA < 0.03,"Low", ifelse(KTR_cohort_characteristics$OTA < 0.14, "Median","High"))
KTR_cohort_characteristics$OTA_level <- factor(KTR_cohort_characteristics$OTA_level, levels = c("High", "Low", "Median"))
table(KTR_cohort_characteristics$OTA_level)

# Fit the Cox proportional hazards model to associate OTA exposure to survival rate
cox_model_OTA <- coxph(Surv(Duration_death, Death) ~ Sex+ Age + OTA + KREA_BLO+ ALB_URV + BMI + Hemoglobin + `Calcineurin inhibitors` + `PPI use` + `Antibiotics` , data = KTR_cohort_characteristics )
summary(cox_model_OTA)

cox_model_OTA_level <- coxph(Surv(Duration_death, Death) ~ Sex+ Age + OTA_level +KREA_BLO+ ALB_URV + BMI + Hemoglobin + `Calcineurin inhibitors` + `PPI use` + `Antibiotics` , data = KTR_cohort_characteristics )
summary(cox_model_OTA_level)

km_fit_OTA <- survfit(Surv(Duration_death, Death) ~ OTA_level, data = KTR_cohort_characteristics)

## Plot Kaplan-Meier curves of OTA exposure group
KM_OTA_all_groups <- 
  survfit2(Surv(Duration_death, Death) ~ OTA_level, data = KTR_cohort_characteristics) %>% 
  ggsurvfit(linewidth=2, alpha = 0.8) +
  add_censor_mark(color="grey30", size=0.9, shape=3) +
  add_risktable(risktable_stats = "n.risk") +
  labs(y="Survival probability (%)", x="Survival time (days)") +
  scale_y_continuous(expand = c(0.025, 0), limits = c(0, 1), label = scales::label_percent()) +
  scale_x_continuous(expand = c(0.035, 0)) +
  scale_color_manual(values=c("#c2a5cf","#5aae61","#e08214")) + 
  theme(
    legend.position = "right",
    legend.text=element_text(size=12, color="black"),
    strip.background=element_blank(),
    strip.text=element_text(size=12, color="black"),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(fill="white"),
    panel.border=element_rect(colour="black", fill=NA, size=1),
    plot.title=element_text(size=15, color="black"),
    axis.title=element_text(size=15, color="black"),
    axis.text.y=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12,color="black")) 
KM_OTA_all_groups

# Fit the Cox proportional hazards model to associate balance value to survival rate

balance_df_OTA # This is the dataset of phenotype data + computed balance value
balance_df_OTA$death_event <- as.numeric(balance_df_OTA$censored_2022_RECOVERL1)
balance_df_OTA$death_event <- ifelse(balance_df_OTA$death_event == 1,0,1)

KM_balance <- 
  survfit2(Surv(Duration_death, death_event) ~ balance_cat, data = balance_df_OTA) %>% 
  ggsurvfit(linewidth=2) +
  add_censor_mark(color="black", size=0.9, shape=3) +
  add_risktable(risktable_stats = "n.risk") +
  labs(y="Survival probability (%)", x="Survival time (days)") +
  scale_y_continuous(expand = c(0.025, 0), limits = c(0, 1), label = scales::label_percent()) +
  scale_x_continuous(expand = c(0.035, 0)) +
  scale_color_manual(values=c("#fdb863","#b2abd2")) + 
  theme(
    legend.position = "none",
    legend.text=element_text(size=12, color="black"),
    strip.background=element_blank(),
    strip.text=element_text(size=12, color="black"),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(fill="white"),
    panel.border=element_rect(colour="black", fill=NA, size=1),
    plot.title=element_text(size=15, color="black"),
    axis.title=element_text(size=15, color="black"),
    axis.text.y=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12,color="black")) 
KM_balance

## Cox regression for balance value
fit_OTA <- coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_value + sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh + ALB_URV + KREA_BLO, data = balance_df_OTA, id = CODE_BIOBANK)
summary(fit_OTA)

## Cox regression for balance group (High - Low)
fit_cat_OTA <- coxph(Surv(Duration_death, censored_2022_RECOVERL1) ~ balance_cat + sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh + ALB_URV + KREA_BLO, data = balance_df_OTA, id = CODE_BIOBANK)
summary(fit_cat_OTA)



