library(tidyverse)
library(phyloseq)
library(codacore)
library(driver)


#############
set.seed(456)

microbe # metagenonomics data
metab # metabolomics data
metadata_surv # phenotype data

# Preparing data for training and testing
Y_pretx <- otu_table(microbe, taxa_are_rows=F)
Y_pretx <- as.matrix(zCompositions::cmultRepl(Y_pretx, label=0, method="CZM", z.delete = F)) # multiplicative 0 replacement 

Z_pretx <- otu_table(metab, taxa_are_rows=F) + 1

Y_pretx <- miniclo(Y_pretx)

Z_pretx <- miniclo(Z_pretx)

Y_Z_pretx <- cbind(Y_pretx,Z_pretx)

train_num_r_pretx <- floor(table(sample_data(metadata_surv)$OTA_encoded)*0.8)["1"]
train_num_nr_pretx <- floor(table(sample_data(metadata_surv)$OTA_encoded)*0.8)["0"]

train_r_index_pretx <- sample(1:length(sample_data(metadata_surv)$OTA_encoded[sample_data(metadata_surv)$OTA_encoded=="1"]), train_num_r_pretx)

train_nr_index_pretx <- sample(1:length(sample_data(metadata_surv)$OTA_encoded[sample_data(metadata_surv)$OTA_encoded=="0"]), train_num_nr_pretx)

## Preparing outcomes
yTrain_pretx <- rbind(sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="1"][train_r_index_pretx,],
                      sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="0"][train_nr_index_pretx,]) %>% pull(OTA_encoded)

yTrain_pretx <- factor(yTrain_pretx, levels = c(0,1))

yTest_pretx <- rbind(sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="1"][-train_r_index_pretx,],
                     sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="0"][-train_nr_index_pretx,]) %>% pull(OTA_encoded)

yTest_pretx <- factor(yTest_pretx, levels = c(0,1))

## Preparing predictors
xTrain_pretx <- Y_Z_pretx[(rbind(sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="1"][train_r_index_pretx,],
                                 sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="0"][train_nr_index_pretx,]) %>% 
                             pull(CODE_BIOBANK)),]

xTest_pretx <- Y_Z_pretx[(rbind(sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="1"][-train_r_index_pretx,],
                                sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="0"][-train_nr_index_pretx,]) %>% 
                            pull(CODE_BIOBANK)),]

# Adjust for covariates
dfTrain_pretx <- rbind(sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="1"][train_r_index_pretx,],
                       sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="0"][train_nr_index_pretx,])
partial_pretx <- glm(OTA_encoded ~ BMI + age + sex + KREA_BLO, data=dfTrain_pretx, family='binomial')

dfTest_pretx <- rbind(sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="1"][-train_r_index_pretx,],
                      sample_data(metadata_surv)[sample_data(metadata_surv)$OTA_encoded=="0"][-train_nr_index_pretx,])

# Training the model with 100-times repeated 5-folds CV
codacore_replicates_pretx <- purrr::rerun(100, codacore(xTrain_pretx, yTrain_pretx, logRatioType='balances', offset=predict(partial_pretx), lambda=1, maxBaseLearners=1))

codacore_replicates_test_auc_pretx <- lapply(seq_len(length(codacore_replicates_pretx)), function(x) pROC::auc(pROC::roc(yTest_pretx, predict(partial_pretx, newdata=dfTest_pretx)+predict(codacore_replicates_pretx[[x]], xTest_pretx, isLogits=TRUE), quiet=TRUE))[1])

codacore_to_plot_pretx <- vector("list", 100)

for(i in 1:100) {
  codacore_to_plot_pretx[[i]] <-
    getTidyTable(codacore_replicates_pretx[[i]]) %>%
    select(-logRatioIndex, Features=Name) %>%
    mutate(`Model replicate`=i,
           testAUC=codacore_replicates_test_auc_pretx[[i]])
}


p_pretx <- 
  bind_rows(codacore_to_plot_pretx) %>% 
  left_join(annotated_long[,c("Features","compound_description","compound_class", "HMDB_ID", "merged_pathways", "pathway_name")], by = "Features" ) %>% 
  mutate(feature=ifelse(str_starts(Features, "s__"),
      Features,coalesce(compound_description,HMDB_ID)),
         Side=factor(Side, levels=c("Numerator","Denominator")),
         SUB_PATHWAY = if_else(
      str_starts(Features, "s__"),
      NA_character_,
    as.factor(pathway_name)
    )) %>% 
  ggplot(aes(x=`Model replicate`, y=feature, color=SUB_PATHWAY)) + 
  geom_point(aes(size=testAUC)) + 
  scale_size_continuous(range = c(1, 3)) +
  facet_grid(vars(Side), scales="free_y") +
  #scale_color_manual(values=c("#17154f", "#2f357c", "#6c5d9e", "#9d9cd5", "#b0799a", "#f6b3b0", "#e48171", "#bf3729", "#e69b00", "gray")) +
  theme(
    legend.key=element_rect(fill='white'),
    legend.text=element_text(size=10, color="black"),
    strip.background=element_blank(),
    strip.text=element_text(size=15, color="black"),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(fill="white"),
    panel.border=element_rect(colour="black", fill=NA, size=1),
    plot.title=element_text(size=15, color="black"),
    axis.title=element_text(size=15, color="black"),
    axis.ticks.length=unit(0.10,"cm"), 
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_text(size=12, color="black")) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol=1))

p_pretx

# Plot ROC

library(pROC)
library(ggplot2)
library(dplyr)
library(purrr)

# 1) Get the partial-model (logistic) linear predictor on the test set
partial_lp <- predict(
  partial_pretx, 
  newdata = dfTest_pretx, 
  type    = "link"         # gives log-odds
)

# 2) For each codacore replicate, predict with i = 1, isLogits = TRUE
roc_list <- map(
  codacore_replicates_pretx,
  ~{
    # predict codacore contribution
    cc_lp <- predict(
      object  = .x,
      newx = xTest_pretx,
      i       = 1,         # first (and only) learner
      isLogits = TRUE
    )
    # total logit = logistic offset + codacore log-ratio
    total_lp <- partial_lp + cc_lp

    # build ROC (use legacy axes: x = 1−spec, y = sens)
    roc(
      response   = yTest_pretx,
      predictor  = total_lp,
      levels     = rev(levels(yTest_pretx)),  # ensure correct pos/neg
      quiet      = TRUE
    )
  }
)

# 3) Extract AUCs
auc_vec_OTA_taxa_metab <- map_dbl(roc_list, auc)
summary(auc_vec_OTA_taxa_metab) 

mean(auc_vec_OTA_taxa_metab)
sd(auc_vec_OTA_taxa_metab)

# 4a) Plot all ROCs + mean curve
roc_df <- imap_dfr(
  roc_list,
  ~tibble(
    replicate = .y,
    fpr       = 1 - .x$specificities,
    tpr       = .x$sensitivities
  )
)

roc_summary_OTA_taxa_metab <- roc_df %>%
  group_by(fpr) %>%
  summarise(
    mean_tpr = mean(tpr, na.rm = TRUE),
    lo_tpr   = quantile(tpr, 0.025, na.rm = TRUE),
    hi_tpr   = quantile(tpr, 0.975, na.rm = TRUE)
  )

roc_summary_OTA_taxa_metab

# 4) Plot the average ROC with 95% ribbon
ggplot(roc_summary_EnnB_taxa, aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = lo_tpr, ymax = hi_tpr),
              fill = "#2b8cbe", alpha = 0.2) +
  geom_line(colour = "#2b8cbe", size = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  theme_classic() +
  labs(
    title = "Mean ROC Curve with 95% CI",
    x     = "1 – Specificity (FPR)",
    y     = "Sensitivity (TPR)"
  )
                                             
