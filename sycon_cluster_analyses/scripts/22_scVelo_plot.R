#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(SeuratExtend)
library(tidyverse)


##################################################
#     LOAD SEURAT DATA TO GET CLUSTER COLORS     #
##################################################

load("00_input/Sycon_Seuratv4.Rdata")

cluster_annotation <- read.table("07_notableGenes_clusterAnnotation/cluster_identity.tsv",
                                 header = TRUE) %>%
  mutate(cluster_ID = factor(cluster_ID))
cluster_annotation

# add cluster identity
Sycon[[]] <- Sycon[[]] %>%
  left_join(cluster_annotation,
            by = join_by("seurat_clusters" == "cluster_ID")) %>%
  rename("cluster_identity" = "value")

umap_clusters_raw <- DimPlot(Sycon, group.by = "seurat_clusters")
umap_clusters_raw

umap_clusters <- umap_clusters_raw@data %>%
  mutate(seurat_clusters = fct_infreq(as.factor(seurat_clusters))) %>%
  mutate(seurat_clusters = fct_rev(seurat_clusters)) %>%
  
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters, label = seurat_clusters)) +
  geom_point(size = 1) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 33,
                                                      sort = "diff"))
umap_clusters

cluster_color_list <- ggplot_build(umap_clusters)$data[[1]] %>%
  select(colour, label) %>%
  distinct() %>%
  arrange(desc(label)) %>%
  pull(colour)

blend_with_bg <- function(hex, alpha = 0.5, bg = "white") {
  fg_rgb <- col2rgb(hex) / 255
  bg_rgb <- as.vector(col2rgb(bg) / 255)  # drop dim -> plain length-3 vector, recycles per column
  blended_rgb <- alpha * fg_rgb + (1 - alpha) * bg_rgb
  rgb(blended_rgb[1,], blended_rgb[2,], blended_rgb[3,])
}

cluster_color_list_blended <- blend_with_bg(cluster_color_list, alpha = 0.5, bg = "white")
cluster_color_list_blended

annotation_color <- c("Unknown" = "#d9d9d9",
                      "Choanocytes" = "#67a64e",
                      "Uncertain" = "#999999",
                      "Metabolocyte-_and_pinacocyte-like_early_embryos" = "#747cc9",
                      "Myopeptidocyte-like_choanocytes" = "#ca603f",
                      "Oocytes/early_embryos" = "#48b1a7",
                      "Sclerocytes" = "#c8577b",
                      "Accessory_cells" = "#b49440")

annotation_color_list_blended <- blend_with_bg(annotation_color, alpha = 0.7, bg = "white")
names(annotation_color_list_blended) <- names(annotation_color)
annotation_color_list_blended


#############################
#     PREPARE ADATA OBJ     #
#############################

reticulate::use_condaenv("scvelo_env", required = TRUE)
scv <- reticulate::import("scvelo")
ad <- reticulate::import("anndata")

adata <- ad$read_h5ad("18_scvelo/4_Velocity_scVelo_processed_subset.h5ad")

# add cluster names to dataframe
obs_df <- adata$obs %>%
  rownames_to_column(var = "cell_id") %>%
  mutate(cluster_ID = factor(seurat_clusters)) %>%
  left_join(cluster_annotation, by = "cluster_ID") %>%
  rename(cluster_annotation = value) %>%
  select(-cluster_ID) %>%
  column_to_rownames(var = "cell_id") %>%
  mutate(cluster_annotation = str_replace(cluster_annotation, "_[0-9]+$", ""))
obs_df

# add back information to adata object
adata$obs <- obs_df

categories <- reticulate::py_to_r(adata$obs$cluster_annotation)

ordered_annotation_colors <- annotation_color_list_blended[categories]

adata$uns[["cluster_annotation_colors"]] <- reticulate::r_to_py(as.list(ordered_annotation_colors))


#############################
#     PLOT TRAJECTORIES     #
#############################

scv$pl$velocity_embedding_grid(adata, basis = "X_umap", color = "seurat_clusters",
                               arrow_length = 1.5, arrow_color = "black",
                               palette = cluster_color_list_blended, alpha = 1,
                               size = 20,
                               title = "", legend_loc = 'none',
                               save = "18_scvelo/figures/umap_trajectory_grid.png",
                               figsize = c(3.002, 3.118),  dpi = as.integer(800))

scv$pl$velocity_embedding_stream(adata, basis = "X_umap", color = "cluster_annotation",
                                 alpha = 1,
                                 density = 3, cutoff_perc = 0,
                                 arrow_color = "#222222", arrow_size = 0.5, arrow_style = "-|>",
                                 linewidth = 0.6,
                                 size = 20,
                                 title = "", legend_loc = 'none',
                                 save = "18_scvelo/figures/umap_trajectory_stream.svg",
                                 figsize = c(3.002, 3.118), dpi = as.integer(800))


################################
#     PLOT DYNAMICAL GENES     #
################################

# plot dynamic genes
top_genes <- adata$var['fit_likelihood'] %>%
  arrange(desc(fit_likelihood)) %>%
  head(n = 500) %>%
  rownames_to_column(var = "genes") %>%
  pull(genes)
top_genes

obs_df %>% 
  select(cluster_annotation) %>%
  distinct() %>%
  pull()

g <- scv$pl$heatmap(adata,
                    var_names = top_genes,
                    sortby = 'latent_time',
                    col_color = 'cluster_annotation',
                    n_convolve = as.integer(100),
                    yticklabels = FALSE,
                    figsize = c(2.43, 2.31),
                    show = FALSE)

g$savefig("18_scvelo/figures/heatmap_dynamical_genes.png",
          dpi = as.integer(600),
          bbox_inches = "tight")


#######################################
#     PLOT LENGTHS AND CONFIDENCE     #
#######################################

Velocity_scVelo_processed_subset_S.object <- sceasy::convertFormat("18_scvelo/4_Velocity_scVelo_processed_subset.h5ad",
                                                                   from = "anndata", to = "seurat")

boxplot_panels <- Velocity_scVelo_processed_subset_S.object[[]] %>%
  select(seurat_clusters, velocity_confidence, velocity_length) %>%
  pivot_longer(-seurat_clusters) %>%
  mutate(name = str_replace(name, "velocity_", "Velocity ")) %>%
  group_by(name) %>%
  mutate(p_low  = quantile(value, 0.05, na.rm = TRUE),
         p_high = quantile(value, 0.95, na.rm = TRUE),
         value_clipped = pmin(pmax(value, p_low), p_high),   # clip to 5th-95th percentile
         value_scaled  = scales::rescale(value_clipped, to = c(-1, 1))  # then rescale as before
         ) %>%
  ungroup() %>%
  
  ggplot(aes(x = seurat_clusters, y = value)) +
  geom_jitter(aes(colour = value_scaled), alpha = 0.6,
              size = 0.6, width = 0.1) +
  geom_boxplot(color = "grey10", outliers = FALSE, size = .6) +
  
  scale_color_gradient2(high = "#cc132e", mid = "#f3f3f3", low = "#3b4dbc",
                        midpoint = 0, guide = "none") +
  
  scale_x_discrete(breaks = seq(0, 32, 2),
                   minor_breaks = seq(0, 32)) +
  
  labs(x = "Cell clusters") +
  
  facet_wrap(~name, nrow = 2, scales = "free_y") +
  
  theme_bw(base_size = 18) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), 
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey90", lineend = "round", linewidth = .5),
    panel.grid.minor.x = element_line(color = "grey90", lineend = "round", linewidth = .5),
    panel.grid.minor.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "grey90", linewidth = .5),
    axis.ticks.length = unit(0.10, "cm"),
    axis.text.x = element_text(color = "black",
                               margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(color = "black",
                               margin = margin(t = 0, r = 4, b = 0, l = 0)),
    axis.title.y = element_blank(),
    axis.title.x = element_text(angle = 0, face = "bold",
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", angle = 0, hjust = 0),
    strip.clip = "off"
  )
boxplot_panels

ggsave("18_scvelo/fig2_scVelo_boxplots.pdf",
       boxplot_panels, device = cairo_pdf,
       height = 2.974*2, width = 3.02*2, units = "in", dpi = 300, bg = "white")
ggsave("18_scvelo/fig2_scVelo_boxplots.png",
       boxplot_panels, device = "png",
       height = 2.974*2, width = 3.02*2, units = "in", dpi = 300, bg = "white")
