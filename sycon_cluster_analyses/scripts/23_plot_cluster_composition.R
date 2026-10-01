#!/usr/bin/env Rscript

library(Seurat)
library(SeuratExtend)
library(tidyverse)

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

# load single cell data
load("00_input/Sycon_Seuratv4.Rdata")

DefaultAssay(Sycon) <- "RNA"

Sycon_not_integrated <- DietSeurat(Sycon, assay = "RNA")

DefaultAssay(Sycon) <- "SCT"

Sycon_not_integrated <- Sycon_not_integrated %>%
  SCTransform() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

umap_not_integrated_raw <- DimPlot(Sycon_not_integrated, split.by = "orig.ident")

umap_not_integrated <- umap_not_integrated_raw@data %>%
  ggplot(aes(umap_1, umap_2, col = orig.ident)) +
  geom_point() +
  
  scale_color_manual(values = SeuratExtend::color_pro(4, col.space = "bright"),
                     guide = "none") +
  facet_wrap(~orig.ident) +
  
  labs(col = "Libraries", title = "UMAP projection before integration",
       x = "UMAP 1", y = "UMAP 2") +
  
  theme_bw()
umap_not_integrated

umap_integrated_raw <- DimPlot(Sycon, split.by = "orig.ident")

umap_integrated <- umap_integrated_raw@data %>%
  ggplot(aes(UMAP_1, UMAP_2, col = orig.ident)) +
  geom_point() +
  
  scale_color_manual(values = SeuratExtend::color_pro(4, col.space = "bright"),
                     guide = "none") +
  facet_wrap(~orig.ident) +
  
  labs(col = "Seurat clusters", title = "UMAP projection after integration",
       x = "UMAP 1", y = "UMAP 2") +
  
  theme_bw()
umap_integrated

panel_umaps <- ggpubr::ggarrange(umap_not_integrated, umap_integrated,
                                 common.legend = TRUE, labels = "AUTO",
                                 legend = "right")
panel_umaps

colpot <- Sycon[[]] %>%
  count(seurat_clusters, orig.ident) %>%
  ggplot(aes(seurat_clusters, n, fill = orig.ident)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = SeuratExtend::color_pro(4, col.space = "bright")) +
  labs(fill = "Libraries", title = "Cell composition of clusters",
       x = "Seurat clusters", y = "Proportion of cells") +
  theme_bw()

panel_final <- ggpubr::ggarrange(panel_umaps, colpot,
                                 nrow = 2, labels = c("", "C"),
                                 heights = c(1, 0.5))
panel_final

ggsave("20_cluster_composition/cluster_composition.png",
       panel_final, device = "png",
       width = 12, height = 10, dpi = 300, units = "in", bg = "white")
ggsave("20_cluster_composition/cluster_composition.pdf",
       panel_final, device = cairo_pdf,
       width = 12, height = 10, dpi = 300, units = "in", bg = "white")
