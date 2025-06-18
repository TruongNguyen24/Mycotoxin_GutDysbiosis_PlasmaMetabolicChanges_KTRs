library(readxl)
library(dplyr)
library(patchwork)
library(ggbreak)
library(ggpubr)
library(rstatix)
library(ggplot2)

# Two-sided Wilcoxon test
wilcox_test(data, OTA ~ group)
wilcox_test(data, EnnB ~ group)
wilcox_test(data, CIT ~ group)

wilcox_effsize(data, OTA ~ group)

# Create the violin plot with testing results
p_OTA <- ggplot(data_compare, aes(x = group, y = OTA_log, fill = group)) +
  # Violin plot to show the overall density, including the spike at zero
  geom_violin(trim = FALSE, alpha = 0.5, color = "grey50") +
  # Boxplot overlay to highlight median and interquartile range (IQR)
  #geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  # Jittered points to display individual observations
  geom_jitter(aes(color = group), width = 0.15, alpha = 0.5, size = 1) +
  # Add the mean as a red point
  stat_summary(fun = mean, 
               geom = "point", 
               aes(shape = "Mean"), 
               size = 2, 
               color = "grey10", 
               show.legend = TRUE) +
  # Add median as a blue point and map its shape to "Median" for the legend
  stat_summary(fun = median, 
               geom = "point", 
               aes(shape = "Median"), 
               size = 2, 
               color = "grey10", 
               fill = "grey10", 
               show.legend = TRUE) +
  # Add 95% CI error bars for the mean and map their linetype to "CI" for the legend
  stat_summary(fun.data = mean_cl_normal, 
               geom = "errorbar",
               width = 0.1, 
               color = "grey10", 
               show.legend = TRUE) +
  scale_fill_manual(values = c("Donor" = "#ce1256", "KTR" = "#02818a")) +
  scale_color_manual(values = c("Donor" = "#ce1256", "KTR" = "#02818a")) +
  scale_y_break(c(0.4, 0.8), scales=0.25, ticklabels=c(0.4, 0.8), expand = TRUE) +
  annotate("text", x = 1.5, y = 0.25, label = sprintf("p < %.3f (*)", 0.000), vjust = 0) +
  labs(title = "OTA", y = "Concentration (µg/L)", x = "") +
  theme_minimal() + 
  theme(legend.box = "vertical")
