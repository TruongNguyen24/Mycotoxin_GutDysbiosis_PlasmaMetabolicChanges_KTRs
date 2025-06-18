library(tidyverse) 
library(vegan)     
library(ape)        
library(microbiome)
library(zCompositions)   # zero replacement
library(compositions)    # clr()
library(ggrepel)  

# Preprocessing
## Shorten taxonomy names to keep only the species names
rownames(microbiome_data) <- str_replace_all(rownames(microbiome_data), "^.*s__", "s__")

## Transpose abundance table: rows are samples and columns are species
microbiome_data <- t(microbiome_data)

## Rescale the table to ensure row sums are equal to 1
microbiome_data <- apply(microbiome_data, 1, function(x) x / sum(x)) %>% t() %>% as.data.frame()

rowSums(microbiome_data) ## Check if the row sums of 'microbiome_data' are all equal to 1

# Ecological parameters
## Alpha diversity
alpha_diversity <- data.frame(
  Shannon  = vegan::diversity(microbiome_data, index = "shannon"),
  Richness = rowSums(microbiome_data != 0)
)

wilcox.test(alpha_diversity$Shannon ~ phenotype_data$group)

data_plot <- cbind(phenotype_data,alpha_diversity)

comparisons <- list(c("Case", "Control"))
p <- ggplot(data_plot, aes(x = group, y = Shannon, color = group)) +
  geom_violin(trim = FALSE, fill = NA, linewidth = 1) +           # Only outline of violin plot
  geom_boxplot(width = 0.1, fill = NA, outlier.shape = NA, linewidth = 1) + # Only outline of boxplot
  geom_jitter(width = 0.15, size = 1, alpha = 0.2) +
  scale_color_manual(values = c("Control" = "#92c5de",  
                                "Case" = "#02818a")) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    aes(label    = paste0("10^", floor(log10(..p..)))),
    parse       = TRUE,
    tip.length  = 0.01
  ) +
  theme_classic() +
  labs(
       x = "",
       y = "Shannon diversity index") +
  theme(legend.position = "none")
print(p)

## Beta diversity

rowSums(microbiome_data) # To test if compositional data
tab_nz <- cmultRepl(microbiome_data, label = 0, method = "CZM", z.delete = FALSE) # zero replacement

clr_mat <- clr(tab_nz) # centred log-ratio transform

pca <- prcomp(clr_mat, scale. = FALSE) # PCA 

pca$x # prepare scores for ggplot

scores <- as_tibble(pca$x[, 1:2], rownames = "sample")   # PC1 & PC2 only
scores <- as.data.frame(scores)
rownames(scores) <- scores$sample
scores$group <- as.factor(data_plot$group)
eig    <- pca$sdev^2
varpad <- round(eig / sum(eig)*100, 1)                 # % variance

### Now build the plot of Aitchison distance between case and control
p <- ggplot(scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.5) +
  xlab(paste0("PC1 (", round(varpad[1], 2), "% variance)")) +
  ylab(paste0("PC2 (", round(varpad[2], 2), "% variance)")) +
  scale_color_manual(values = c("Donor" = "grey",  
                                "Survived" = "#02818a",  
                                "Deceased" = "#ce1256")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),                    # Remove grid lines
    axis.line = element_line(color = "black"),       # Add axis lines
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  coord_fixed(ratio = 1) +
  stat_ellipse(aes(group = group), type = "t", level = 0.95, linetype = 1) +
  # Add centroids for each group (here, shown as black crosses)
  geom_point(data = centroids, aes(x = centroid_axis1, y = centroid_axis2),
           shape = 21, size = 4, fill = c("Donor" = "grey",  
                                "Survived" = "#02818a",  
                                "Deceased" = "#ce1256"), color = "black", stroke = 1) 
  # Draw the circle for the control group average distance
  #geom_circle(data = circle_df, aes(x0 = x0, y0 = y0, r = r),
   #           color = "blue", linetype = "dashed", inherit.aes = FALSE)

print(p)

#### For the Control group itself, the distance will be 0.
centroids2 <- centroids %>% 
  mutate(dist_to_control = sqrt((centroid_axis1 - control_centroid$centroid_axis1)^2 +
                                (centroid_axis2 - control_centroid$centroid_axis2)^2))

print(centroids2)

### PERMANOVA
microbiome_dist_euclidean <- vegdist(microbiome_data_merged_sorted, method = "euclidean")

permanova_test <- adonis2(microbiome_dist_euclidean ~ Age + Sex + BMI + OTA_Presence + EnnB_Presence + TeA_Presence + dum_calcineurin.inh + dum_antibiotics +  dum_ppi, phenotype_data, permutations = 999, by = "term")
