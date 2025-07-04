# load dependencies

library(phyloseq)
library(tidyverse)
library(patchwork)
library(fido) # for Pibble model
library(wec) # for weighted sum contrasts
library(zCompositions) # for 0-imputation 
library(driver) # for a few CoDa functions; devtools::install_github("jsilve24/driver")
library(ggtext)


microbe # Scaled microbiome data ( samples (rows) x species (columns): row sums are equal to 1))
metadata_surv # phenotype data


# Formula

f <- c("sex", "age", "BMI", "OTA_encoded","EnnB_encoded","TeA_encoded", "KREA_BLO", "dum_calcineurin.inh", "dum_antibiotics","dum_ppi")
f <- reformulate(termlabels=f)

# Run Pibble

X <- t(model.matrix(f, data=metadata_surv))
dim(X)    
Y <- microbe
Y <- driver::miniclo(Y)
Y <- as.matrix(zCompositions::cmultRepl(Y, label=0, method="CZM", z.delete = F)) # multiplicative 0 replacement 
Y <- t(Y)

N <- ncol(Y) # number of samples
D <- nrow(Y) # number of features (categories)
    
# Specify uninformative (default) priors 
upsilon <- D+3
Omega <- diag(D)
G <- cbind(diag(D-1), -1)
Xi <- (upsilon-D)*G%*%Omega%*%t(G)
Theta <- matrix(0, D-1, nrow(X))
Gamma <- diag(nrow(X))
    
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi, n_samples=1000)
    
eta_init <- t(driver::alr(t(Y)))
eta_array <- array(eta_init, dim=c(nrow(eta_init), ncol(eta_init), 1000))
    
posterior <- uncollapsePibble(eta_array, priors$X, priors$Theta, priors$Gamma, priors$Xi, priors$upsilon, seed=2849)
posterior   
# Attach dimnames
dimnames(posterior$Lambda)[[2]] <- rownames(X)
dimnames(posterior$Lambda)[[1]] <- rownames(Y)[-length(rownames(Y))]
dimnames(posterior$Sigma)[[1]] <- dimnames(posterior$Sigma)[[2]] <- rownames(Y)[-length(rownames(Y))]

posterior <- pibblefit(D=D,
                           N=N,
                           Q=nrow(X),
                           coord_system="alr",
                           iter=1000L,
                           alr_base=D,
                           Eta=eta_array,
                           Lambda=posterior$Lambda,
                           Sigma=posterior$Sigma,
                           Y=Y,
                           X=X,
                           names_categories=rownames(Y),
                           names_samples=colnames(Y),
                           names_covariates=rownames(X),
                           upsilon=priors$upsilon,
                           Theta=priors$Theta,
                           Gamma=priors$Gamma,
                           Xi=priors$Xi)
# Change to CLR transform
    posterior_clr_taxa <- to_clr(posterior)
    posterior_clr_taxa
# Attached dimnames
    dimnames(posterior_clr_taxa$Lambda)[[2]] <- rownames(X)
    dimnames(posterior_clr_taxa$Lambda)[[1]] <- rownames(Y)
    dimnames(posterior_clr_taxa$Sigma)[[1]] <- dimnames(posterior_clr_taxa$Sigma)[[2]] <- rownames(Y)

