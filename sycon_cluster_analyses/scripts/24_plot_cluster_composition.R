#!/usr/bin/env Rscript

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(scDblFinder)

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

# load single cell data
load("00_input/Sycon_Seuratv4.Rdata")


####################################
#     PLOT CLUSTER COMPOSITION     #
####################################

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
  dplyr::count(seurat_clusters, orig.ident) %>%
  ggplot(aes(seurat_clusters, n, fill = orig.ident)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = SeuratExtend::color_pro(4, col.space = "bright")) +
  labs(fill = "Libraries", title = "Cell composition of clusters",
       x = "Seurat clusters", y = "Proportion of cells") +
  theme_bw()


#############################
#     DOUBLET DETECTION     #
#############################

sycon_sce <- as.SingleCellExperiment(Sycon)

sycon_sce <- scDblFinder(sycon_sce, clusters = "seurat_clusters", samples = "orig.ident")

VlnPlot(sycon_seu, group.by = "orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2, pt.size = 0) + theme(legend.position = 'right')

barplot_doublets <- sycon_seu[[]] %>%
  dplyr::count(seurat_clusters, scDblFinder.class) %>%
  ggplot(aes(x = seurat_clusters, y = n, fill = scDblFinder.class)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = SeuratExtend::color_pro(n = 2, col.space = "bright", sort = "diff")) +
  labs(x = "Seurat clusters", y = "Proportion",
       title = "Doublet rates in each cluster",
       fill = "scDblFinder\nclassification") +
  theme_bw(base_size = 12)
barplot_doublets

panel_final <- ggpubr::ggarrange(panel_umaps,
                                 ggpubr::ggarrange(colpot, barplot_doublets,
                                                   align = "hv", nrow = 2,
                                                   labels = c("C", "D")),
                                 nrow = 2)
panel_final

ggsave("20_cluster_composition/cluster_composition.png",
       panel_final, device = "png",
       width = 12, height = 12, dpi = 300, units = "in", bg = "white")
ggsave("20_cluster_composition/cluster_composition.pdf",
       panel_final, device = cairo_pdf,
       width = 12, height = 12, dpi = 300, units = "in", bg = "white")

doublet_per_cluster <- table(sycon_seu$seurat_clusters, sycon_seu$scDblFinder.class)
doublet_per_cluster

# chi test with simulated p-values to account for low counts
chi <- chisq.test(doublet_per_cluster, simulate.p.value = TRUE, B = 1e5)
chi$expected

# cramer's V to test for effect size
n <- sum(tab)
V <- sqrt(chi$statistic / (n * (min(dim(tab)) - 1)))
V

round(chi$stdres, 2)

