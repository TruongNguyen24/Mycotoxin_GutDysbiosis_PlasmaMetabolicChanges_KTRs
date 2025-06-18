library(zCompositions)

microbiome_data # This is scaled microbiome data (row sum equal to 1)
phenotype_data # Phenotype data

ftbl <- microbiome_data

ftbl <- as.matrix(zCompositions::cmultRepl(ftbl, label=0, method="CZM", z.delete=F)) # Bayesian 0 imputation (this will generate a warning about features and samples with >=80% 0s. Here, we choose to keep those, hence z.delete=FALSE)

# Select SGBs for balance computing
# For OTA
numerator_taxa <- c("s__Intestinimonas_butyriciproducens",
                    "s__.Collinsella._massiliensis",
                    "s__Escherichia_coli",
                    "s__Bacteroides_thetaiotaomicron",
                    "s__Enorma_massiliensis",
                    "s__Coprobacter_fastidiosus",
                    "s__Phascolarctobacterium_faecium",
                    "s__Ruminococcaceae_bacterium_D5",
                    "s__Fusicatenibacter_saccharivorans",
                    "s__Firmicutes_bacterium_CAG_170")

numerator_taxa_index <- match(numerator_taxa, colnames(ftbl))

denominator_taxa <- c("s__Blautia_obeum",
                      "s__Lachnospira_pectinoschiza",
                      "s__Actinomyces_sp_HMSC035G02",
                      "s__Rothia_mucilaginosa",
                      "s__Enterococcus_faecium",
                      "s__Butyricimonas_virosa",
                      "s__Actinomyces_odontolyticus",
                      "s__Veillonella_atypica",
                      "s__Clostridium_saccharolyticum",
                      "s__Ruminococcus_lactaris",
                      "s__Blautia_obeum",
                      "s__Lachnospira_pectinoschiza",
                      "s__Eubacterium_rectale",
                      "s__Alistipes_inops")

denominator_taxa_index <- match(denominator_taxa, colnames(ftbl))

balance_df <- phenotype_data %>% mutate(balance_value=NA) 

# Compute balance
for(sample_i in rownames(ftbl)){
  balance_df[sample_i,]$balance_value <- log(exp(mean(log(ftbl[sample_i,numerator_taxa_index])))) - log(exp(mean(log(ftbl[sample_i,denominator_taxa_index]))))  
}

# Categorize balance based on median 
balance_df$balance_cat <- ifelse(balance_df$balance_value>quantile(balance_df$balance_value)[3], "high", "low")
balance_df$balance_cat <- as.factor(balance_df$balance_cat)
