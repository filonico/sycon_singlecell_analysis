#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

suppressMessages({library(tidyverse)
  library(hdWGCNA)
  library(SeuratExtend)})


######################
#     LOAD INPUT     #
######################

Sycon_blobOnly <- readRDS("13_recluster_blob/Sycon_blobOnly.Rds") %>%
  UpdateSeuratObject()

DimPlot(Sycon_blobOnly)


#####################
#     FUNCTIONS     #
#####################

plot_module_UMAP <- function(module_name) {
  
  # first get the coordinates of points in the umap
  umap <- Sycon_blobOnly@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    
    # then add information from meta data, including the score for each module
    left_join(Sycon_blobOnly@meta.data %>% rownames_to_column()) %>%
    select(rowname, umap_1, umap_2, all_of(module_name)) %>%
    
    # pivot data longer
    pivot_longer(-c(rowname, umap_1, umap_2), names_to = "module", values_to = "expr") %>%
    arrange(expr) %>%
    # group_by(module) %>%
    
    # plot data
    ggplot(aes(umap_1, umap_2)) +
    geom_point(aes(col = expr),
               # alpha = 0.5, 
               size = 0.7
    ) +
    
    # define color gradient rules
    # note that in the original hdWGCNA script they make the range symmetrical
    # is it legit? Isn't this affecting visualization?
    # here, we do not make the range symmetrical
    scale_color_gradient2(low = "grey75", mid = "grey95", high = module_name,
                          midpoint = 0) +
    
    ggtitle(module_name) +
    
    # make the plot squared
    coord_fixed() +
    
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
          panel.background = element_blank(),
          panel.grid = element_blank()
    )
  
  return(umap)
}


##########################################
#     hdWGCNA with standard pipeline     #
##########################################

DefaultAssay(Sycon_blobOnly) <- "SCT"
Sycon_blobOnly

# standard hdWGCNA pipeline
# pipeline with variable features is better than the other ("fraction")
# it retrieves better defined clusters
Sycon_blobOnly <- Sycon_blobOnly %>%
  SetupForWGCNA(gene_select = "variable",
                wgcna_name = "RNA_variableFeatures")

# Sycon <- SetupForWGCNA(Sycon,
#                        gene_select = "fraction",
#                        fraction = 0.05,
#                        wgcna_name = "RNA_fraction")

Sycon_blobOnly@misc$active_wgcna <- "RNA_variableFeatures"

# compute metacells by splitting both by samples and seurat clusters
# reduce number of cells to consider a cluster and the number of overlapping cells to 10
Sycon_blobOnly <- Sycon_blobOnly %>%
  MetacellsByGroups(group.by = c("orig.ident", "seurat_clusters_new"),
                    ident.group = "seurat_clusters_new",
                    min_cells = 10, k = 10,
                    assay = "RNA",
                    wgcna_name = "RNA_variableFeatures")

Sycon_blobOnly <- Sycon_blobOnly %>%
  NormalizeMetacells()

# it is important to use the RNA assay for DGE analysis, in order to use raw counts
# instead of normalized and/or integrated matrices
Sycon_blobOnly <- Sycon_blobOnly %>%
  SetDatExpr(group_name = levels(Sycon_blobOnly[[]]$seurat_clusters_new),
             group.by = "seurat_clusters_new",
             assay = "RNA",
             layer = "data")

# test different soft powers
Sycon_blobOnly <- Sycon_blobOnly %>%
  TestSoftPowers()

# # plot the results
# plotSoftPower_list <- PlotSoftPowers(Sycon)
# 
# # assemble with patchwork
# patchwork::wrap_plots(plotSoftPower_list, ncol = 2)

# construct co-expression network
Sycon_blobOnly <- Sycon_blobOnly %>%
  ConstructNetwork(tom_outdir = "13_recluster_blob/03_hdWGCNA/01_RNA_assay/",
                   tom_name = "RNA_variableFeatures", overwrite_tom = TRUE)

# PlotDendrogram(Sycon)

# this is necessary to include gorup by samples  in the next step
Sycon_blobOnly <- Sycon_blobOnly %>%
  ScaleData(features = VariableFeatures(Sycon_blobOnly))

# compute all MEs in the full single-cell dataset
Sycon_blobOnly <- Sycon_blobOnly %>%
  ModuleEigengenes(group.by.vars = "orig.ident")

# compute eigengene-based connectivity (kME):
Sycon_blobOnly <- Sycon_blobOnly %>%
  ModuleConnectivity(group.by = "seurat_clusters_new")

# make a featureplot of hMEs for each module
# here, plotted ranges are symmetrical (as hdWGCNA default settings)
RNA_featurePlot_list <- Sycon_blobOnly %>%
  ModuleFeaturePlot(features = "hMEs",
                    order = TRUE)

# plot patchworked plots
patchwork::wrap_plots(RNA_featurePlot_list,
                      ncol = 2)

# extract hME values and add them to meta data
hMEs <- Sycon_blobOnly %>%
  GetMEs(harmonized = TRUE)
Sycon_blobOnly@meta.data <- cbind(Sycon_blobOnly@meta.data, hMEs)

# create a custom UMAP for each module
# here, plotted ranges are not corrected, thus not symmetrical
RNA_umap_list <- list()
for (i in names(RNA_featurePlot_list)) {
  RNA_umap_list[[i]] <- plot_module_UMAP(i)
}

# plot patchworked plots
umaps <- print(patchwork::wrap_plots(RNA_umap_list, ncol = 2))
umaps

ggsave("13_recluster_blob/03_hdWGCNA/01_RNA_assay/gene_module_umaps.pdf",
       umaps,  device = cairo_pdf,
       dpi = 300, height = 6, width = 8, units = ("in"), bg = 'white')

ggsave("13_recluster_blob/03_hdWGCNA/01_RNA_assay/gene_module_umaps.png",
       umaps,  device = "png",
       dpi = 300, height = 6, width = 8, units = ("in"), bg = 'white')

# get the list of genes per module
gene_list <- Sycon_blobOnly@misc$RNA_variableFeatures$wgcna_modules %>%
  group_by(module) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  { setNames(.$genes, .$module) }

# write the list of genes per module to separate files
for (module in names(gene_list)) {
  file_path <- file.path("13_recluster_blob/03_hdWGCNA/01_RNA_assay/", paste0(module, "_RNA_module_genes.ls"))
  writeLines(gene_list[[module]],
             file_path)
}

saveRDS(Sycon_blobOnly, "13_recluster_blob/03_hdWGCNA/Sycon_blobOnly_hdWGCNA.Rds")

RNA_modules <- GetModules(Sycon_blobOnly)
RNA_mods <- levels(RNA_modules$module); RNA_mods <- RNA_mods[RNA_mods != 'grey']

RNA_dotplot <- Sycon_blobOnly %>%
  DotPlot(features = RNA_mods, group.by = 'seurat_clusters_new')

RNA_dotplot_custom <- RNA_dotplot$data %>%
  mutate(id = factor(id,
                     levels = c("2", "12",
                                "11", "9", "18",
                                "0", "4", "5", "7", "14", "17", "20",
                                "1", "3", "6", "8", "10", "13", "15", "16", "19")),
         features.plot = factor(features.plot,
                                levels = c("yellow", "brown", "turquoise", "blue"))) %>%
  filter(pct.exp > 0) %>%
  ggplot(aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
  
  scale_y_discrete(limits = rev) +
  scale_size_continuous(range = c(1,7)) +
  
  xlab("gene co-expression\nmodules") +
  ylab("Cell clusters") +
  labs(size = "Percent\nexpressed", fill = "Average\nexpression") +
  
  scale_fill_viridis_c(option = "plasma") +
  
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

RNA_dotplot_custom

ggsave("13_recluster_blob/03_hdWGCNA/01_RNA_assay/RNA_gene_module_dotplot.pdf",
       RNA_dotplot_custom,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

ggsave("13_recluster_blob/03_hdWGCNA/01_RNA_assay/RNA_gene_module_dotplot.png",
       RNA_dotplot_custom,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')
