library(tidyverse)
library(readxl)
library(compositions)
library(robCompositions)
library(dendextend)
library(vegan)
library(ggvegan)
library(ggcorrplot)
library(patchwork)

#te data preparation --------

te_raw <- read_excel("data/LM_woda_statytyska.xlsx", sheet = "TE_1") %>% 
  mutate(ID = as.factor(ID),
         layer = as.factor(layer))

te_long <- te_raw %>% 
  pivot_longer(cols = c(Li:pH), names_to = "variable", values_to = "value")

te_boxplots <- ggplot(te_long, aes(y = value)) +
  geom_boxplot() +
  facet_wrap(.~variable, scales = "free")

te_clr_prep <- te_raw %>% 
  select(!c(ID, layer, pH))

te_clr <- clr(te_clr_prep) 

te_clr_long <- as_tibble(te_clr)%>% 
  mutate(ID = te_raw$ID) %>%
  pivot_longer(cols = c(Li:U), names_to = "variable", values_to = "value")

te_boxplot_clr <- ggplot(te_clr_long, aes(y = value)) +
  geom_boxplot() +
  facet_wrap(.~ variable, scales = "free")

#te pca--------

te_clr_z <- scale(te_clr)

te_pca <- rda(te_clr_z) 
screeplot(te_pca, bstick = TRUE) #first three PCs are signitficant according to "broken stick" model

te_fort <- fortify(te_pca, axes = c(1,2), scaling = "sites")

te_sites <- te_fort[te_fort$score %in% "sites",] %>% 
  mutate(ID = te_raw$ID)
te_sp <- te_fort[te_fort$score %in% "species",]

te_inertcomp <- inertcomp(te_pca) #contribution of each variable to the total inertia
te_sel_sp <- tibble(variable = rownames(te_inertcomp), inertcomp = te_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp),
         inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
te_sp_red <- te_sp[which(te_sp$label %in% te_sel_sp$variable), ] #leave 10 variables with the highest contrubution to the total inertia

te_ve_prep <- te_pca$CA$eig / te_pca$tot.chi * 100
(te_PC1_ve <- round(((te_ve_prep / sum(te_ve_prep))[c(1)]) * 100, digits = 1))#42.8%
(te_PC2_ve <- round(((te_ve_prep / sum(te_ve_prep))[c(2)]) * 100, digits = 1))#22.0%

te_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", te_PC2_ve, "%)", sep = ""), x = paste("PC1 (", te_PC1_ve, "%)", sep = ""),
       title = "PCA biplot for TE") +
  geom_segment(data = te_sp,
               color = "black", linewidth = 0.3,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = te_sites, aes(x = PC1, y = PC2, color = te_raw$pH, shape = te_raw$layer), size = 4) +
  scale_color_viridis_c() +
  ggrepel::geom_text_repel(data = te_sites, color = "black",
                           size = 4, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = ID)) +
  ggrepel::geom_text_repel(data = te_sp, color = "red",
                           size = 4, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6, linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6, linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.2)) +
#  guides(color = guide_legend(title = "pH")) +
  coord_equal() +
  theme(legend.position = "right") +
  labs(color = "pH",
       shape = "Layer")

ggsave(filename="figures/te.svg", 
       plot = te_pca_plot, 
       device = svg, 
       width = 6, 
       height = 4, 
       units = "in")

ggsave(filename="figures/te.jpeg", 
       plot = te_pca_plot, 
       device = jpeg, 
       width = 6, 
       height = 4, 
       units = "in")

#ree data preparation------
ree_raw <- read_excel("data/LM_woda_statytyska.xlsx", sheet = "REE_1") %>% 
  mutate(ID = as.factor(ID),
         layer = as.factor(layer))

ree_long <- ree_raw %>% 
  pivot_longer(cols = c(La:pH), names_to = "variable", values_to = "value")

ree_boxplots <- ggplot(ree_long, aes(y = value)) +
  geom_boxplot() +
  facet_wrap(.~variable, scales = "free")

ree_clr_prep <- ree_raw %>% 
  select(!c(ID, layer, pH))

ree_clr <- clr(ree_clr_prep) 

ree_clr_long <- as_tibble(ree_clr)%>% 
  mutate(ID = ree_raw$ID) %>%
  pivot_longer(cols = c(La:Lu), names_to = "variable", values_to = "value")

ree_boxplot_clr <- ggplot(ree_clr_long, aes(y = value)) +
  geom_boxplot() +
  facet_wrap(.~ variable, scales = "free")

#ree pca--------

ree_clr_z <- scale(ree_clr)

ree_pca <- rda(ree_clr_z) 
screeplot(ree_pca, bstick = TRUE) #first two PCs are signitficant according to "broken stick" model

ree_fort <- fortify(ree_pca, axes = c(1,2), scaling = "sites")

ree_sites <- ree_fort[ree_fort$score %in% "sites",] %>% 
  mutate(ID = ree_raw$ID)
ree_sp <- ree_fort[ree_fort$score %in% "species",]

ree_inertcomp <- inertcomp(ree_pca) #contribution of each variable to the total inertia
ree_sel_sp <- tibble(variable = rownames(ree_inertcomp), inertcomp = ree_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp),
         inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
ree_sp_red <- ree_sp[which(ree_sp$label %in% ree_sel_sp$variable), ] #leave 10 variables with the highest contrubution to the total inertia

ree_ve_prep <- ree_pca$CA$eig / ree_pca$tot.chi * 100
(ree_PC1_ve <- round(((ree_ve_prep / sum(ree_ve_prep))[c(1)]) * 100, digits = 1))#52.4%
(ree_PC2_ve <- round(((ree_ve_prep / sum(ree_ve_prep))[c(2)]) * 100, digits = 1))#21.0%

ree_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", ree_PC2_ve, "%)", sep = ""), x = paste("PC1 (", ree_PC1_ve, "%)", sep = ""),
       title = "PCA biplot for ree") +
  geom_segment(data = ree_sp,
               color = "black", linewidth = 0.3,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = ree_sites, aes(x = PC1, y = PC2, color = ree_raw$pH, shape = ree_raw$layer), size = 4) +
  scale_color_viridis_c() +
  ggrepel::geom_text_repel(data = ree_sites, color = "black",
                           size = 4, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = ID)) +
  ggrepel::geom_text_repel(data = ree_sp, color = "red",
                           size = 4, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6, linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6, linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.2)) +
  #  guides(color = guide_legend(title = "pH")) +
  coord_equal() +
  theme(legend.position = "right") +
  labs(color = "pH",
       shape = "Layer")

ggsave(filename="figures/ree.svg", 
       plot = ree_pca_plot, 
       device = svg, 
       width = 6, 
       height = 4, 
       units = "in")

ggsave(filename="figures/ree.jpeg", 
       plot = ree_pca_plot, 
       device = jpeg, 
       width = 6, 
       height = 4, 
       units = "in")
