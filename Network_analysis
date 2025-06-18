library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)

posterior_clr_orthus # This is the outcome for fido::orthus
posterior_clr_metab #This is the outcome for fido::pibble for metabolomics
posterior_clr_taxa # This is the outcome for fido::pibble for metagenomics

rho_summary # This is the result of joint modelling of taxa-metabolites

edge_attributes_faecal_mgs_plasma_mbs <- 
  rho_summary  %>% 
  mutate(feature_ii=feature_i, feature_jj=feature_j, 
         feature_i=str_remove_all(string = feature_i, paste(c("^taxa_","^metab_"), collapse = "|")),
         feature_j=str_remove_all(string = feature_j, paste(c("^taxa_","^metab_"), collapse = "|"))) %>% 
  filter(sign(p2.5)==sign(p97.5)) 

edge_attributes_faecal_mgs_plasma_mbs2 <- 
  edge_attributes_faecal_mgs_plasma_mbs %>% 
  mutate(pair=case_when((grepl(x = feature_ii, pattern = "^taxa_") & grepl(x = feature_jj, pattern = "^taxa_")) ~ "species-species", 
                        (grepl(x = feature_ii, pattern = "^metab_") & grepl(x = feature_jj, pattern = "^metab_")) ~ "metabolite-metabolite",
                         TRUE ~ "mixed")) %>% 
select(-feature_ii, -feature_jj)

mb_mb_pair_median <- edge_attributes_faecal_mgs_plasma_mbs2 %>% 
  filter(pair=="metabolite-metabolite") %>% pull(p50) %>% abs() %>% quantile(prob=0.75)

sp_sp_pair_median <- edge_attributes_faecal_mgs_plasma_mbs2 %>% 
  filter(pair=="species-species") %>% pull(p50) %>% abs() %>% quantile(prob=0.75)

mixed_pair_median <- edge_attributes_faecal_mgs_plasma_mbs2 %>% 
  filter(pair=="mixed") %>% pull(p50) %>% abs() %>% quantile(prob=0.75)


edge_attributes_faecal_mgs_plasma_mbs_mb_mb <- 
  edge_attributes_faecal_mgs_plasma_mbs2 %>% 
  filter(pair=="metabolite-metabolite" & abs(p50)>mb_mb_pair_median)

edge_attributes_faecal_mgs_plasma_mbs_sp_sp <- 
  edge_attributes_faecal_mgs_plasma_mbs2 %>% 
  filter(pair=="species-species" & abs(p50)>sp_sp_pair_median)

edge_attributes_faecal_mgs_plasma_mbs_sp_mb <- 
  edge_attributes_faecal_mgs_plasma_mbs2 %>% 
  filter(pair=="mixed" & abs(p50)>mixed_pair_median)

edge_attributes_faecal_mgs_plasma_mbs3 <- 
  bind_rows(edge_attributes_faecal_mgs_plasma_mbs_mb_mb,
            edge_attributes_faecal_mgs_plasma_mbs_sp_sp,
            edge_attributes_faecal_mgs_plasma_mbs_sp_mb)
            

edge_attributes <- edge_attributes_faecal_mgs_plasma_mbs3

edge_attributes$sign <- ifelse(edge_attributes$p50>=0, "positive", "negative")

edge_list <- edge_attributes[,c("feature_i","feature_j")]

g <- graph_from_data_frame(edge_list, directed=FALSE)
g
V(g)$group <- ifelse(startsWith(V(g)$name,"s__")==TRUE, "taxon", "metabolite") 
V(g)$shape <- ifelse(V(g)$group=="taxon", "circle", "square") 


g_edge_cols <- setNames(c("#e28686","#9494ea"),c("positive","negative"))
E(g)$color <- g_edge_cols[factor(edge_attributes$sign, levels=c("positive","negative"))]
E(g)$sign <- edge_attributes$sign

rescale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}

E(g)$transformed_weights <- 1 - ((edge_attributes$p50 + 1) / 2)
E(g)$transformed_weights[E(g)$transformed_weights==0] <- 0.000000001
E(g)$weight <- rescale(abs(edge_attributes$p50), 0, max(abs(edge_attributes$p50),na.rm=T), 1, 2)

# Modify node names conditionally
V(g)$name <- 
  ifelse(grepl("^f__", V(g)$name),  # Check if the name starts with "f__"
                    sub("^f__[^|]+ \\| ", "", V(g)$name),  # Remove the first part up to " | "
                    V(g)$name)  # Leave other names unchanged


# This is for species-species only
ll <- layout_with_kk(g, weights = E(g)$transformed_weights)
rownames(ll) <- V(g)$name

V(g)$size <- 2
plot(g, layout=ll, edge.color=E(g)$color, vertex.frame.color="black", vertex.label.dist=0.0, vertex.label.cex=0.00001,vertex.label.color="black")

g_components <- components(g, mode="strong")
biggest_cluster_id <- which.max(g_components$csize)
vert_ids <- V(g)[g_components$membership == biggest_cluster_id]
table(g_components$membership)

g_largest_components <- subgraph(g, vids = vert_ids)
springlass_cluster <- cluster_spinglass(g_largest_components, weights = E(g_largest_components)$weight)
unique(membership(springlass_cluster))
V(g_largest_components)$modules <- as.character(membership(springlass_cluster))

leiden_cluster <- cluster_leiden(g_largest_components,
                                 weights = E(g_largest_components)$weight_3,
                                 objective_function = "modularity",
                                 resolution = 1, n_iterations = 500)
unique(membership(leiden_cluster))
V(g_largest_components)$modules <- as.character(membership(leiden_cluster))
length(V(g_largest_components)$modules)

as_tbl_graph(g_largest_components) %>%
  ggraph(layout = "manual", x = ll[c(V(g_largest_components)$name),1], y = ll[c(V(g_largest_components)$name),2]) +
  geom_edge_link(aes(color = sign, width=weight), alpha = 0.5) +
  geom_node_point(aes(fill = group, shape=group), color="black", size=4) +
  geom_node_text(aes(label = name), size=2) +
  theme_graph() +
  scale_shape_manual(values=c(23,21)) +
  theme(legend.position = "none") +
  facet_nodes(~modules) +
  scale_edge_colour_manual(values=c("negative"="#92c5de", "positive"="#f4a582")) +
  scale_fill_manual(values=c("taxon"="#8c510a", "metabolite"="#01665e"))

as_tbl_graph(g_largest_components) %>%
  ggraph(layout = "manual", x = ll[c(V(g_largest_components)$name),1], y = ll[c(V(g_largest_components)$name),2]) +
  geom_edge_link(aes(color = sign, width=weight), alpha = 0.5) +
  geom_node_point(aes(fill = group, shape=group), color="black", size=4) +
  geom_node_text(aes(label = name), size=0.5) +
  theme_graph() +
  scale_shape_manual(values=c(23,21)) +
  theme(legend.position = "none") +
  facet_nodes(~modules) +
  scale_edge_colour_manual(values=c("negative"="#b8e186", "positive"="#f1b6da")) +
  scale_fill_manual(values=c("taxon"="#8c510a", "metabolite"="#01665e"))

# Sub modules
select_sub_graph <- V(g_largest_components)[V(g_largest_components)$modules=="1"]

sub_graph <- subgraph(g_largest_components, vids = select_sub_graph)
sub_graph

module <- as_tbl_graph(sub_graph) %>% 
  ggraph(layout = "manual", x = ll[c(V(sub_graph)$name),1], y = ll[c(V(sub_graph)$name),2]) +
  geom_edge_link(aes(color = sign, width=weight), alpha = 0.3) +
  geom_node_point(aes(fill = group, shape=group), color="black", size=6) +
 # geom_node_point(fill = "#1b7837", shape=23, color="black", size=6) +

  geom_node_label(aes(label = name), fill="gray", color="black", size=3, repel = T, label.padding = unit(0.2, "lines"), max.overlaps=100) +
  scale_shape_manual(values=c(23,21)) +
  scale_edge_colour_manual(values=c("negative"="#9494ea", "positive"="#e28686")) +
  scale_fill_manual(values=c("taxon"="#762a83", "metabolite"="#1b7837")) +
  theme_graph(fg_text_colour = 'white', base_family = 'Helvetica') +
  theme(legend.position = "none")
module 
