#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggvenn)


##########################
#     PRE-PROCESSING     #
##########################

input_path <- "./00_input/scRNAseqs_rawTables/"

# Amphimedon queenslandica files
aque_files <- list(UMItable = paste0(input_path, "A_queenslandica/GSM3021561_Aque_Amphimedon_queenslandica_adult_UMI_table.txt.gz"),
                   metacell_assignments = paste0(input_path, "A_queenslandica/Amphimedon_adult_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "A_queenslandica/Aque_adult_metacell_annotation.tsv"))

# create the Seurat object from UMI table
aque_S.object <- read.table(file = gzfile(aque_files$UMItable),
                            header = TRUE, sep = "\t") %>%
  Seurat::CreateSeuratObject(project = "aque_adult", min.cells = 0, min.features = 0,
                             meta.data = read.table(aque_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(aque_files$metacell_annotation, header = TRUE, sep = "\t")) %>%
                               rename(cell_type = ID) %>%
                               column_to_rownames(var = "Cell"))

# plot feature statistics before filtering
vln_plot_nfeature_before_filtering <- VlnPlot(aque_S.object,
                                              features = "nFeature_RNA")
vln_plot_ncount_before_filtering <- VlnPlot(aque_S.object,
                                            features = "nCount_RNA")

FeatureScatter(aque_S.object,
               feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# filter based on features and counts
aque_S.object <- aque_S.object %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
           nCount_RNA > 200 & nCount_RNA < 11000)

# plot feature statistics after filtering
vln_plot_nfeature_after_filtering <- VlnPlot(aque_S.object,
                                             features = "nFeature_RNA")
vln_plot_ncount_after_filtering <- VlnPlot(aque_S.object,
                                           features = "nCount_RNA")

df_to_plot <- vln_plot_nfeature_after_filtering@data %>%
              rownames_to_column(var = "cells") %>%
  left_join(vln_plot_ncount_after_filtering@data %>%
              rownames_to_column(var = "cells"),
            by = c("cells", "ident")) %>%
  rename("nCount_RNA_after" = "nCount_RNA",
         "nFeature_RNA_after" = "nFeature_RNA") %>%
  pivot_longer(-c(cells, ident)) %>%
  mutate(name = str_replace(name, "_RNA", "RNA")) %>%
  separate(name, into = c("category", "stage"), sep = "_") %>%
  drop_na()
df_to_plot

violin_plots_faceted <- df_to_plot %>%
  mutate(stage = factor(stage, levels = c("before", "after"))) %>%
  ggplot(aes(x = ident, y = value)) +
  geom_jitter(col = "grey80", size = 0.1, width = 0.3, height = 0.1) +
  geom_violin(fill = alpha("white", 0)) +
  facet_wrap(~ category, scales = "free",
             labeller = labeller(
               category = c(nCountRNA = "nCount RNA",
                            nFeatureRNA = "nFeature RNA"))) +
  labs(y = "Counts", x = "Libraries") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = NA),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
violin_plots_faceted

# save the panel to file
ggsave("21_aque_plots/feature_plot_panel.pdf",
       violin_plots_faceted, device = cairo_pdf,
       width = 4, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("21_aque_plots/feature_plot_panel.png",
       violin_plots_faceted, device = "png",
       width = 4, height = 5, dpi = 300, units = "in", bg = "white")


######################
#     CLUSTERING     #
######################

aque_S.object <- aque_S.object %>%
  SCTransform(verbose = TRUE)

aque_S.object <- aque_S.object %>%
  RunPCA(verbose = TRUE)

ElbowPlot(aque_S.object, ndims = 50)

aque_S.object <- aque_S.object %>%
  RunUMAP(dims = 1:30, verbose = TRUE)

aque_S.object <- aque_S.object %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 1)

aque_S.object[[]] <- aque_S.object[[]] %>%
  mutate(Metacell = as_factor(Metacell))

umap <- DimPlot2(aque_S.object, features = "Metacell", cols = "bright",
                 theme = list(labs(title = expression(paste(bolditalic("Amphimedon queenslandica"),
                                                            bold(" cell atlas")))),
                              theme_classic(), theme_umap_arrows(),
                              theme(plot.title = element_text(hjust = 0.5))))
umap

final_panel <- ggpubr::ggarrange(violin_plots_faceted, umap, 
                                 nrow = 1, labels = "AUTO", widths = c(0.5, 1))
final_panel

ggsave("21_aque_plots/panel_amphimedon.pdf",
       final_panel, device = cairo_pdf,
       height = 4, width = 12, dpi = 300, units = ("in"), bg = 'white')
ggsave("21_aque_plots//panel_amphimedon.png",
       final_panel, device = "png",
       height = 4, width = 12, dpi = 300, units = ("in"), bg = 'white')
