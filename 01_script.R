library(cluster)
library(fpc)
library(tidyverse)
library(readxl)
library(compositions)
library(robCompositions)
library(dendextend)
ev_clust_function <- ("https://raw.githubusercontent.com/ethen8181/machine-learning/master/clustering_old/clustering/clustering_functions.R")
library(vegan)
library(ggvegan)
library(ggcorrplot)

#read data-------
LM_osady_statystyka <- read_excel("data/LM_osady_statystyka.xlsx", sheet = "1") %>% #sa 2 arkusze. w pierwszym ree sa rozbite, w drugim sa polaczone.
  mutate(LakeID = as.factor(paste("L", LakeID, sep = "")))

#box plots------
df_muzakow_long <- LM_osady_statystyka %>%  #LM_osady_statystyka bylo wczytane z excela; zwrocic uwage na nazwy
mutate(LakeID = as.factor(paste("L", LakeID, sep = "")),
       age = as.factor(age)) %>%
  pivot_longer(cols = c(Li:pH_bot, MD), names_to = "variable", values_to = "value")

muzakow_bp <- ggplot(df_muzakow_long, aes(value)) +
  geom_boxplot() +
  facet_wrap(.~ variable, scales = "free")

df_muzakow_clr_prep <- LM_osady_statystyka %>% 
  select(!c(LakeID, age))

df_muzakow_clr <- clr(df_muzakow_clr_prep) 

df_muzakow_clr_long <- as_tibble(df_muzakow_clr)%>% 
  mutate(ID = LM_osady_statystyka$LakeID) %>%
  pivot_longer(cols = c(Li:MD), names_to = "variable", values_to = "value")

muzakow_bp_clr <- ggplot(df_muzakow_clr_long, aes(value)) +
  geom_boxplot() +
  facet_wrap(.~ variable, scales = "free")

#clustering-------
df_muzakow_cluster_prep <- LM_osady_statystyka %>%  #LM_osady_statystyka bylo wczytane z excela; zwrocic uwage na nazwy
  select(!c(pH_surf, age, MD, LakeID)) #testowalem rozne warianty; tylko przy tym jeziora o spojnym wieku sie grupuja tak ladnie; wydaje sie to sensowne

df_muzakow_cluster <- as.data.frame(df_muzakow_cluster_prep)

#tranformation, i.e. first ilr and then standardisation to 0 mean and 1 stdv follows Templ et al. (2008): who wrote: Of all tested data transformations simple log-transformation and isometric logratio transformation delivered the most reliable results. Because geochemical data are always compositional data the isometric logratio transformation is preferable, it can, however, not be used for variable clustering because the direct relationship to the elements is lost. If the data show very different magnitude for the different variables (e.g., major, minor and trace elements mixed) the variables need to be standardised to mean 0 and variance 1. Following standardisation all variables will have the same influence on the results of cluster analysis.
df_muzakow_ilr <- ilr(df_muzakow_cluster)
df_muzakow_ilr_z <- scale(df_muzakow_ilr)

rownames(df_muzakow_ilr_z) <- paste(LM_osady_statystyka$LakeID, "(", LM_osady_statystyka$age, ")", sep = "")

df_muzakow_hclust <- hclust(dist(df_muzakow_ilr_z), "ward.D2")

plot(df_muzakow_hclust, hang = -1)
rect.hclust(df_muzakow_hclust, k = 4)

dend <- as.dendrogram(df_muzakow_hclust)

dend <- dend %>%
  set("branches_k_color", k = 4)

plot(dend)

criteria <- CHCriterion(data = df_muzakow_ilr_z, kmax = 10, clustermethod = "hclust", method = "ward.D2") #shows two popular indices that help selecting optimal number of clusters

criteria$plot #at the left plot showing Calinski-Harabasz Index we look for the largest number, while at the right we seek for "elbow". Right plot is not very helpful here. C-H index suggests solution with 6 clusters

#correlation matrix----

muzakow_cor <- round(cor(df_muzakow_cluster_prep, method = "spearman"), 1)

ggcorrplot(muzakow_cor, hc.order = TRUE, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"))

#PCA---------
#justification for clr from Filzmoser et al. (2018) (book on analyses of compositional data): Throughout the book, the analysis of compositional data is done with orthonormal coordinates, because they form the most natural way for a representation with respect to the Aitchison geometry. Traditionally, PCA is an exception with this respect in compositional data analysis, because clr coefficients are commonly applied instead. This dates back already to the work of Aitchison (1983). Here, a compromise is chosen: Although, due to methodological reasons, PCA will be developed in orthonormal coordinates, also the advantages of its computation and interpretation directly in clr coefficients will be shown. Namely, the latter approach turns out to be computationally much simpler and can thus be advantageously applied in practice

df_muzakow_clr <- clr(df_muzakow_cluster)
df_muzakow_clr_z <- scale(df_muzakow_clr)
  
muzakow_pca <- rda(df_muzakow_clr_z) 
screeplot(muzakow_pca, bstick = TRUE) #only first PC is signitficant according to "broken stick" model

muzakow_fort <- fortify(muzakow_pca, axes = c(1,2), scaling = "sites")

muzakow_sites <- muzakow_fort[muzakow_fort$score %in% "sites",] %>% 
  mutate(LakeID = LM_osady_statystyka$LakeID)
muzakow_sp <- muzakow_fort[muzakow_fort$score %in% "species",]

muzakow_inertcomp <- inertcomp(muzakow_pca) #contribution of each variable to the total inertia
muzakow_sel_sp <- tibble(variable = rownames(muzakow_inertcomp), inertcomp = muzakow_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp),
         inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
muzakow_sp_red <- muzakow_sp[which(muzakow_sp$label %in% muzakow_sel_sp$variable), ] #leave 10 variables with the highest contrubution to the total inertia

muzakow_ve_prep <- muzakow_pca$CA$eig / muzakow_pca$tot.chi * 100
(muzakow_PC1_ve <- round(((muzakow_ve_prep / sum(muzakow_ve_prep))[c(1)]) * 100, digits = 1))#52.2%
(muzakow_PC2_ve <- round(((muzakow_ve_prep / sum(muzakow_ve_prep))[c(2)]) * 100, digits = 1))#13.1%

muzakow_pca_plot_age <- ggplot() +
  labs(y = paste("PC2 (", muzakow_PC2_ve, "%)", sep = ""), x = paste("PC1 (", muzakow_PC1_ve, "%)", sep = "")) +
  geom_segment(data = muzakow_sp,
               color = "black", linewidth = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = muzakow_sites, aes(x = PC1, y = PC2, color = LM_osady_statystyka$age), size = 3) +
#  scale_color_brewer() +
  ggrepel::geom_text_repel(data = muzakow_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = LakeID)) +
  ggrepel::geom_text_repel(data = muzakow_sp, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6, linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6, linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.2)) +
  guides(color = guide_legend(title = "Age")) +
  theme(legend.position = "right")

muzakow_pca_plot_depth <- ggplot() +
    labs(y = paste("PC2 (", muzakow_PC2_ve, "%)", sep = ""), x = paste("PC1 (", muzakow_PC1_ve, "%)", sep = "")) +
    geom_segment(data = muzakow_sp,
                 color = "black", linewidth = 0.7,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = 
                   grid::arrow(length = grid::unit(0.25, "cm"))) +
    geom_point(data = muzakow_sites, aes(x = PC1, y = PC2, color = LM_osady_statystyka$MD), size = 3) +
    scale_color_viridis_c() +
    ggrepel::geom_text_repel(data = muzakow_sites, color = "black",
                             size = 2.5, segment.alpha = 0,
                             aes(x = PC1, y = PC2, 
                                 label = LakeID)) +
    ggrepel::geom_text_repel(data = muzakow_sp, color = "black",
                             size = 2.5, segment.alpha = 0,
                             aes(x = PC1, y = PC2, 
                                 label = label)) +
    geom_vline(xintercept = 0, color = 'black', linewidth = 0.6, linetype=2) + 
    geom_hline(yintercept = 0, color = 'black', linewidth = 0.6, linetype=2) +
    theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_line(color = "grey80", linewidth = 0.2)) +
    guides(color = guide_legend(title = "Age")) +
  theme(legend.position = "right")

muzakow_pca_plot_age + muzakow_pca_plot_depth
  
