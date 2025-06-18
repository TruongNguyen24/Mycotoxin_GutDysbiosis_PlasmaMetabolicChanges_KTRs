microbe # Scaled microbiome data ( samples (rows) x species (columns): row sums are equal to 1))
metab # Scaled metabolomics data
metadata_surv # phenotype data

f <- c("sex", "age", "BMI" ,"KREA_BLO", "dum_calcineurin.inh", "dum_antibiotics","dum_ppi")

f <- reformulate(termlabels=f)
X <- t(model.matrix(f, data=metadata_surv))

X

Y <- microbe
Y <- driver::miniclo(Y)
Y <- as.matrix(zCompositions::cmultRepl(Y, label=0, method="CZM", z.delete = F))
       
Z <- metab
Z <- driver::miniclo(Z)
Z <- as.matrix(zCompositions::cmultRepl(Z, label=0, method="CZM", z.delete = F))
dim(Z)    

eta_init <- t(driver::alr(Y))
rownames(eta_init) <- colnames(Y)[-length(colnames(Y))]
Z_init <- t(driver::clr(Z)) 
    
eta_init_array <- array(eta_init, dim=c(nrow(eta_init), ncol(eta_init), 1000))
    
eta_Z <- t(cbind(t(eta_init), t(Z_init)))
    
eta_Z_array <- array(eta_Z, dim=c(nrow(eta_Z), ncol(eta_Z), 1000))
dimnames(eta_Z_array)[[1]] <- c(colnames(Y)[-length(colnames(Y))], colnames(Z))
dimnames(eta_Z_array)[[2]] <-  rownames(Y)
    
N <- ncol(t(Y))
P <- nrow(t(Z))
Q <- nrow(X)
D <- nrow(t(Y)) 
    
upsilon <- (D-1+P)+10
Xi <- diag(D-1+P)
GG <- cbind(diag(D-1), -1)
Xi[1:(D-1), 1:(D-1)] <- GG%*%diag(D) %*% t(GG)
Xi <- Xi * (upsilon-D-P) 
Gamma <- diag(Q)
Theta <- matrix(0, D-1+P, Q)
    
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi, n_samples=1000)

posterior <- uncollapsePibble(eta=eta_Z_array, X=priors$X, Theta=priors$Theta, Gamma=priors$Gamma, Xi=priors$Xi, upsilon=priors$upsilon, seed=2849)
    
posterior2 <- 
        orthusfit(
        D=D,
        N=N,
        P=P,
        Q=Q,
        coord_system="alr",
        iter=1000L,
        alr_base=D,
        Z=t(Z),
        Y=t(Y),
        Eta=eta_init_array,
        Lambda=posterior$Lambda,
        Sigma=posterior$Sigma,
        X=X,
        names_covariates=rownames(X),
        names_Zdimensions=rownames(t(Z)),
        names_categories=rownames(t(Y)),
        names_samples=colnames(t(Y))
      )
    
posterior_clr_orthus <- to_clr(posterior2)
    
    # Attached dimnames
dimnames(posterior_clr_orthus$Lambda)[[2]] <- rownames(X)
dimnames(posterior_clr_orthus$Lambda)[[1]] <- c(colnames(Y), colnames(Z))
dimnames(posterior_clr_orthus$Sigma)[[1]] <- dimnames(posterior_clr_orthus$Sigma)[[2]] <- c(colnames(Y), colnames(Z))

posterior_clr_orthus # This is the outcome for fido::orthus
posterior_clr_metab # This is the outcome for fido::pibble for metabolites
posterior_clr_taxa # # This is the outcome for fido::orthus for taxa

xcorTaxa <- posterior_clr_taxa$Sigma
dim(xcorTaxa)
xcorMetab <- posterior_clr_metab$Sigma
dim(xcorMetab)
xcorOrthus <- posterior_clr_orthus$Sigma[1:D, (D+1):((D-1+P)+1),]
dim(xcorOrthus)

D <- posterior_clr_orthus$D
P <- posterior_clr_orthus$P

xcorTaxa.mean <- apply(xcorTaxa, c(1,2), mean)
colnames(xcorTaxa.mean) <- rownames(xcorTaxa.mean) <- paste0("taxa_",colnames(xcorTaxa.mean))
  
xcorMetab.mean <- apply(xcorMetab, c(1,2), mean)
colnames(xcorMetab.mean) <- rownames(xcorMetab.mean) <- paste0("metab_",colnames(xcorMetab.mean))
  
xcorOrthus.mean <- apply(xcorOrthus, c(1,2), mean)
  
mega_V <- matrix(NA, nrow=D+P, ncol=D+P)
rownames(mega_V) <- colnames(mega_V) <- c(rownames(xcorTaxa.mean), rownames(xcorMetab.mean))
  
foo <- function(x) mega_V
mega_V_array <- sapply(1:1000, foo, simplify = 'array')
for(i in 1:1000){
    print(i)
    mega_V_array[1:D,1:D,i] <- xcorTaxa[,,i]
    mega_V_array[(D+1):(D+P),(D+1):(D+P),i] <- xcorMetab[,,i]
    mega_V_array[1:D,(D+1):(D+P),i] <- xcorOrthus[,,i]
    mega_V_array[(D+1):(D+P),1:D,i] <- t(xcorOrthus[,,i])
  }

dim(mega_V_array)

# compute proportionality

compute_proportionality_fast <- function(cov_array) {
  # Get dimensions
  n_features <- dim(cov_array)[1]
  n_samples <- dim(cov_array)[3]
  
  print(n_features)
  print(n_samples)
  # Initialize output array
  rho_array <- array(NA, dim = c(n_features, n_features, n_samples))
  
  # Compute proportionality for each posterior sample
  for (k in 1:n_samples) {
    print(k)
    cov_matrix <- cov_array[, , k]
    
    # Extract variances (diagonal elements)
    sigma_diag <- diag(cov_matrix)
    
    # Compute proportionality using vectorized operations
    denom <- outer(sigma_diag, sigma_diag, "+")  # Sum of variances
    rho_mat <- (2 * cov_matrix) / denom  # Compute full rho matrix
    
    # Extract lower triangle (including diagonal) and set upper triangle to NA
    rho_mat[upper.tri(rho_mat, diag = TRUE)] <- NA
    
    # Store result
    rho_array[, , k] <- rho_mat
  }
  
  return(rho_array)
}

rho_ij <- compute_proportionality_fast(mega_V_array)

rho_ij

# gather results
library(fido)
rho_summary <- driver::gather_array(rho_ij, rho, feature_i, feature_j, iter) %>%
    mutate(feature_i=c(rownames(xcorTaxa.mean), rownames(xcorMetab.mean))[feature_i], feature_j=c(rownames(xcorTaxa.mean), rownames(xcorMetab.mean))[feature_j]) %>% 
    drop_na() %>% 
    group_by(feature_i, feature_j) %>% 
    summarise_posterior(rho) %>% 
    arrange(mean)

rho_summary   %>% 
  filter(sign(p2.5)==sign(p97.5)) 
dim(rho_summary) 
