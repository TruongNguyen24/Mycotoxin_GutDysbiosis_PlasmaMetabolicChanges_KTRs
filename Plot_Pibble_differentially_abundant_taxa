

plot_pibble <- function(focus_var=focus_var, pibble_fit=pibble_fit, output_taxa=T, siglevel=siglevel) {
  
  if(output_taxa==T) {
    
    sig_decreasing <- 
      as.data.frame.table(pibble_fit$Lambda) %>%
      filter(Var2==focus_var) %>%
      mutate(idx=as.numeric(as.factor(Var1))) %>%
      select(idx, Freq, x=Var2, y=Var1) %>%
      group_by(idx) %>%
      ggdist::median_qi(Freq, .width = c(0.5, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97)) %>%
      mutate(species=rep(rownames(pibble_fit$Y),7),
             .width=factor(.width)) %>%
      arrange(species) %>%
      mutate(species=reorder(factor(species),Freq)) %>% 
      pivot_wider(species, names_from=.width, values_from=.upper) %>%
      select(species, p50=`0.5`, p75=`0.75`, p80=`0.8`, p85=`0.85`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) < 0) %>%
      mutate(species=factor(species)) %>% 
      select(species) %>% 
      pull() %>% 
      levels()
    
    sig_increasing <- 
      as.data.frame.table(pibble_fit$Lambda) %>%
      filter(Var2==focus_var) %>%
      mutate(idx=as.numeric(as.factor(Var1))) %>%
      select(idx, Freq, x=Var2, y=Var1) %>%
      group_by(idx) %>%
      ggdist::median_qi(Freq, .width = c(0.5, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97)) %>%
      mutate(species=rep(rownames(pibble_fit$Y),7),
             .width=factor(.width)) %>%
      arrange(species) %>%
      mutate(species=reorder(factor(species),Freq)) %>% 
      pivot_wider(species, names_from=.width, values_from=.lower) %>%
      select(species, p50=`0.5`, p75=`0.75`, p80=`0.8`, p85=`0.85`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) > 0) %>%
      mutate(species=factor(species)) %>% 
      select(species) %>% 
      pull() %>% 
      levels() %>% 
      rev()
    
    sig_taxa <- list(sig_decreasing=sig_decreasing, sig_increasing=sig_increasing)
    
    return(sig_taxa)
    
  } else {
    
    sig_decreasing <- 
      as.data.frame.table(pibble_fit$Lambda) %>%
      filter(Var2==focus_var) %>%
      mutate(idx=as.numeric(as.factor(Var1))) %>%
      select(idx, Freq, x=Var2, y=Var1) %>%
      group_by(idx) %>%
      ggdist::median_qi(Freq, .width = c(0.5, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97)) %>%
      mutate(species=rep(rownames(pibble_fit$Y),7),
             .width=factor(.width)) %>%
      arrange(species) %>%
      pivot_wider(species, names_from=.width, values_from=.upper) %>%
      select(species, p50=`0.5`, p75=`0.75`, p80=`0.8`, p85=`0.85`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) < 0) %>%
      select(species) %>% 
      pull()
    
    sig_increasing <- 
      as.data.frame.table(pibble_fit$Lambda) %>%
      filter(Var2==focus_var) %>%
      mutate(idx=as.numeric(as.factor(Var1))) %>%
      select(idx, Freq, x=Var2, y=Var1) %>%
      group_by(idx) %>%
      ggdist::median_qi(Freq, .width = c(0.5, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97)) %>%
      mutate(species=rep(rownames(pibble_fit$Y),7),
             .width=factor(.width)) %>%
      arrange(species) %>%
      pivot_wider(species, names_from=.width, values_from=.lower) %>%
      select(species, p50=`0.5`, p75=`0.75`, p80=`0.8`, p85=`0.85`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
      filter(!! rlang::sym(siglevel) > 0) %>%
      select(species) %>% 
      pull()
    
    p <- as.data.frame.table(pibble_fit$Lambda) %>% 
      filter(Var2==focus_var) %>% 
      mutate(idx=as.numeric(as.factor(Var1))) %>% 
      select(idx, Freq, x=Var2, y=Var1) %>% 
      group_by(idx) %>%
      ggdist::median_qi(Freq, .width = c(0.5, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97)) %>%
      mutate(species=rep(rownames(pibble_fit$Y),7)) %>% 
      filter(species %in% c(sig_decreasing, sig_increasing)) %>%
      
      ggplot(aes(y=reorder(factor(species),Freq), x=Freq, xmin=.lower, xmax=.upper)) +
      ggdist::geom_interval(aes(alpha=.width)) +
      
      scale_alpha_continuous("Credible interval", range=c(.7, .17), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
      geom_point(colour = "white",             # <— make points green
             size   = 1) +
      scale_color_brewer() +
      theme(#legend.position="none",
        legend.key=element_rect(fill='white'),
        legend.text=element_text(size=10, color="black"),
        strip.background=element_blank(),
        strip.text=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.ticks.length=unit(0.25,"cm"), 
        axis.text.x=element_text(size=10, color="black"),
        axis.text.y=element_text(size=10, color="black")) +
      labs(x="Log-Ratio Value", y=NULL, title=focus_var) +
      geom_vline(xintercept=0, linetype="dashed", color="darkgray")
    
    return(p)    
  }
}

# To plot the differentially abundant taxa between detectable and undetectable OTA
p1__taxa <- plot_pibble(focus_var = "OTA_encoded1", pibble_fit =   posterior_clr_taxa, siglevel = "p95", output_taxa = F)
p1_20250410_taxa
