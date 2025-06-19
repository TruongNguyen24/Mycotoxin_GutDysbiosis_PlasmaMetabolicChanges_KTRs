compute_balance_baseline <- function(Y=Y, Z=Z, G1=G1, G2=G2, metadata=metadata) {
  
  A_ftbl <- driver::miniclo(Y)
  A_ftbl <- as.matrix(zCompositions::cmultRepl(A_ftbl, label=0, method="CZM", z.delete = F))
  A_ftbl 
  B_ftbl <- driver::miniclo(Z)
  B_ftbl <- as.matrix(zCompositions::cmultRepl(B_ftbl, label=0, method="CZM", z.delete = F))
  B_ftbl
  
  # G1
  A_ftbl_G1 <- A_ftbl[,colnames(A_ftbl) %in% G1]
  B_ftbl_G1 <- B_ftbl[,colnames(B_ftbl) %in% G1]

  f_G1 <- cbind(A_ftbl_G1, B_ftbl_G1)
  #f_G1 <- B_ftbl_G1
  #f_G1 <- A_ftbl_G1
  
  # G2
  A_ftbl_G2 <- A_ftbl[,colnames(A_ftbl) %in% G2]
  B_ftbl_G2 <- B_ftbl[,colnames(B_ftbl) %in% G2]

  f_G2 <- cbind(A_ftbl_G2, B_ftbl_G2)
  #f_G2 <- B_ftbl_G2
  #f_G2 <- A_ftbl_G2

  # Balance 
  term_G1 <- apply(f_G1, 1, function(x) mean(log2(x)))
  term_G2 <- apply(f_G2, 1, function(x) mean(log2(x)))

  #numerator_size <- ncol(f_G1)
  #denominator_size <- ncol(f_G2)

  #normalization_factor <- sqrt((numerator_size * denominator_size) / (numerator_size + denominator_size))

  B_value <- term_G1 - term_G2
  
  #B_value <- normalization_factor * (term_G1 - term_G2)
  
  out <- data.frame(balance_value=B_value)
  out$CODE_BIOBANK <- rownames(out)
  out$balance_cat_median <- ifelse(out$balance_value>quantile(out$balance_value)[3], "high", "low")
  #out$balance_cat2 <- ifelse(out$balance_value>quantile(out$balance_value)[4], "high", "low")
  #out_m <- out %>% left_join(metadata_surv, by = join_by(CODE_BIOBANK)) %>% drop_na()
  out_m  <- left_join(out, metadata_surv, by = "CODE_BIOBANK")
  rownames(out_m) <- out_m$CODE_BIOBANK
  cutpoint <- survminer::surv_cutpoint(out_m, time="Duration_death", event="censored_2022_RECOVERL1", variables="balance_value", minprop = 0.30)

  out <- survminer::surv_categorize(cutpoint, labels=c("low","high")) 
    
  outout <- bind_cols( (out_m %>% dplyr::select(-`Duration_death`, -`censored_2022_RECOVERL1`)), (data.frame(out) %>% rename(balance_cat=balance_value)) ) 
    
  outout$cutpoint <- cutpoint$cutpoint[,1]
  outout$median <- median(outout$balance_value)
  
  return(outout)
}

metab # metabolomics data
microbe # metagenomics data
metadata_surv # phenotype data

V(g_largest_components) # This is the modules created by Leiden Algorithm from Network analyis of taxa-metabolites covariance matrix
g <- g_largest_components
n_modules <- length(table(V(g)$modules))

pairs <- as_tibble(t(combn(min(as.numeric(names(table(V(g)$modules)))):max(as.numeric(names(table(V(g)$modules)))),2)))
pairs <- pairs %>% mutate(median=NA,
                          cutpoint=NA,
                          lm_OTA_estimate=NA,
                          lm_OTA_pval=NA,
                          lm_EnnB_estimate=NA,
                          lm_EnnB_pval=NA,
                          lm_TeA_estimate=NA,
                          lm_TeA_pval=NA,
                          pval_Cox_cutpoint=NA,
                          pval_Cox_median=NA,
                          fdr_Cox_cutpoint=NA,
                          fdr_Cox_median=NA,
                          pval_LR_cutpoint=NA,
                          pval_LR_median=NA,
                          fdr_LR_cutpoint=NA,
                          fdr_LR_median=NA,
                          HR_cutpoint_adj=NA,
                          HR_cutpoint_nonadj=NA,
                          HR_median_adj=NA,
                          HR_median_nonadj=NA,
                          plot_cutpoint=list(NA),
                          plot_median=list(NA)) %>% 
  rename(numerator=V1, denominator=V2)

for(row in 1:nrow(pairs)) {
  
  print(row) 
  
  Numerator <- V(g)[modules==as.character(pairs[row,]$numerator)]$name
  Denominator <- V(g)[modules==as.character(pairs[row,]$denominator)]$name
  
  mdat_f_mbs_balance <- compute_balance_baseline(Y = microbe, Z = metab, G1 = Numerator, G2 = Denominator, metadata = metadata_surv)
 
  pairs[row,]$median <- mdat_f_mbs_balance$median[1]
  pairs[row,]$cutpoint <- mdat_f_mbs_balance$cutpoint[1]
  
  # 2) Fit the linear model for balance_value ~ toxins + covariates
  lm_fit <- lm(
    balance_value ~ OTA_encoded + EnnB_encoded + TeA_encoded +
      sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh,
    data = mdat_f_mbs_balance
  )
  
  lm_summary <- summary(lm_fit)
  
  # 3) Extract coefficients/p-values for each toxin of interest
  #    Make sure these variable names match exactly the columns in your data
  pairs[row, ]$lm_OTA_estimate <- lm_summary$coefficients["OTA_encoded1", "Estimate"]
  pairs[row, ]$lm_OTA_pval     <- lm_summary$coefficients["OTA_encoded1", "Pr(>|t|)"]
  pairs[row, ]$lm_EnnB_estimate <- lm_summary$coefficients["EnnB_encoded1", "Estimate"]
  pairs[row, ]$lm_EnnB_pval     <- lm_summary$coefficients["EnnB_encoded1", "Pr(>|t|)"]
  pairs[row, ]$lm_TeA_estimate  <- lm_summary$coefficients["TeA_encoded1", "Estimate"]
  pairs[row, ]$lm_TeA_pval      <- lm_summary$coefficients["TeA_encoded1", "Pr(>|t|)"]

  mdat_f_mbs_balance$censored_2022_RECOVERL1 <- ifelse(mdat_f_mbs_balance$censored_2022_RECOVERL1 == 1, 1, 0)

  pairs[row,]$plot_cutpoint[[1]] <-
    survfit2(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_cat, data = mdat_f_mbs_balance) 
  
  pairs[row,]$plot_median[[1]] <- 
    survfit2(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_cat_median, data = mdat_f_mbs_balance)  
  
  pairs[row,]$pval_Cox_cutpoint <- summary(coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_value + sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh, data = mdat_f_mbs_balance, id = CODE_BIOBANK ))$coefficients[1,"Pr(>|z|)"]
  
  pairs[row,]$pval_Cox_median <- summary(coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_cat_median + sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh, data = mdat_f_mbs_balance, id = CODE_BIOBANK))$coefficients[1,"Pr(>|z|)"]
  
  pairs[row,]$pval_LR_cutpoint <- survdiff(Surv(Duration_death, censored_2022_RECOVERL1) ~ balance_cat, data = mdat_f_mbs_balance)$pval
  
  pairs[row,]$pval_LR_median <- survdiff(Surv(mdat_f_mbs_balance$Duration_death, mdat_f_mbs_balance$censored_2022_RECOVERL1) ~ mdat_f_mbs_balance$balance_cat_median)$pval
  
  pairs[row,]$HR_cutpoint_adj <- summary(coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_cat_median + sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh , data = mdat_f_mbs_balance, id = CODE_BIOBANK))$coefficients[1,"exp(coef)"]
  
  pairs[row,]$HR_cutpoint_nonadj <- summary(coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_value, data = mdat_f_mbs_balance, id = CODE_BIOBANK))$coefficients[1,"exp(coef)"]

  pairs[row,]$HR_median_adj <- summary(coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_cat_median + sex + age + BMI + dum_ppi + dum_antibiotics + dum_calcineurin.inh , data = mdat_f_mbs_balance, id = CODE_BIOBANK))$coefficients[1,"exp(coef)"]
  
  pairs[row,]$HR_median_nonadj <- summary(coxph(Surv(Duration_death,  censored_2022_RECOVERL1) ~ balance_cat_median, data = mdat_f_mbs_balance, id = CODE_BIOBANK))$coefficients[1,"exp(coef)"]

}
pairs

# FDR correction
                   
pairs$fdr_Cox_cutpoint <- p.adjust(pairs$pval_Cox_cutpoint, method = "BH")
pairs$fdr_Cox_median <- p.adjust(pairs$pval_Cox_median, method = "BH")
pairs$fdr_LR_cutpoint <- p.adjust(pairs$pval_LR_cutpoint, method = "BH")
pairs$fdr_LR_median <- p.adjust(pairs$pval_LR_median, method = "BH")
pairs$lm_OTA_fdr <- p.adjust(pairs$lm_OTA_pval, method = "BH")
pairs$lm_EnnB_fdr <- p.adjust(pairs$lm_EnnB_pval, method = "BH")
pairs$lm_TeA_fdr <- p.adjust(pairs$lm_TeA_pval, method = "BH")

view(pairs)

                   
