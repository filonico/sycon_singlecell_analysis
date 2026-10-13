#!/usr/bin/env Rscript

library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(ggvenn)
library(grid)

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses")


#####################
#     FUNCTIONS     #
#####################

get_markers_and_save <- function(seurat_object, group, output_dir) {
  markers <- seurat_object %>%
    FindAllMarkers(group.by = group, only.pos = TRUE)
  
  markers_list <- markers %>%
    filter(avg_log2FC > 0,
           p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    # summarise(count = n())
    summarise(genes = list(gene), .groups = "drop") %>%
    { setNames(.$genes, .$cluster) }
  
  # write the list of upregulated genes per cluster to a file
  for (cluster_name in names(markers_list)) {
    file_path <- file.path(output_dir,
                           paste0("cluster", cluster_name, "_upregulatedGenes.ls"))
    writeLines(markers_list[[cluster_name]],
               file_path)
  }
  
  return(markers)
}

get_gene_universe <- function(seurat_object, out_filename){
  
  seurat_object@assays$RNA$counts %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "gene") %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    filter(sum > 0) %>%
    select(c("gene")) %>%
    
    write.table(file = out_filename,
                col.names = FALSE, row.names = FALSE, quote = FALSE)
}


###############################
#     DEFINE PLOT THEMES      #
###############################

theme_for_UMAPS <- theme(
  # aspect.ratio = 1,
  plot.background = element_blank(), 
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)


#################
#     INPUT     #
#################

load("00_input/Sycon_Seuratv4.Rdata")

Sycon <- Sycon %>%
  UpdateSeuratObject()


######################
#     RE-CLUSTER     #
######################

DefaultAssay(Sycon) <- "RNA"

# define list of lcusters to remove
clusters_to_remove <- c("15", "17", "19", "21", "24", "26", "27", "25", "30", "32")

# retain only RNA assay and remove cell belonging to above clusters
Sycon_blobOnly <- Sycon %>%
  UpdateSeuratObject() %>%
  DietSeurat(assays = "RNA") %>%
  subset(!(seurat_clusters %in% clusters_to_remove)) 

# split object back into individual batches, to perform integration
Sycon_blobOnly[["RNA"]] <- split(Sycon_blobOnly[["RNA"]],
                                 f = Sycon_blobOnly$orig.ident)

Sycon_blobOnly <- Sycon_blobOnly %>%
  SCTransform()

# Sycon_blobOnly <- Sycon_blobOnly %>%
#   RunPCA()
# 
# Sycon_blobOnly <- Sycon_blobOnly %>%
#   RunUMAP(dims = 1:30)
# 
# DimPlot(Sycon_blobOnly, group.by = "orig.ident")

Sycon_blobOnly <- Sycon_blobOnly %>%
  IntegrateLayers(method = HarmonyIntegration,
                  normalization.method = "SCT", verbose = TRUE)

ElbowPlot(Sycon_blobOnly, ndims = 50) %>%
  ggsave(plot = ., filename = "13_recluster_blob/elbowplot.pdf", device = cairo_pdf)

Sycon_blobOnly <- Sycon_blobOnly %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()

Sycon_blobOnly <- Sycon_blobOnly %>%
  RunUMAP(dims = 1:30, reduction = "harmony")

Sycon_blobOnly@meta.data <- Sycon_blobOnly@meta.data %>%
  rownames_to_column() %>%
  left_join(Sycon@meta.data %>%
              rownames_to_column() %>%
              select(c("rowname", "seurat_clusters")),
            by = c("rowname" = "rowname"),
            # keep = FALSE,
            suffix = c("_new", "_original")) %>%
  column_to_rownames(var = "rowname")

DimPlot(Sycon_blobOnly, group.by = c("orig.ident", "seurat_clusters_new", "seurat_clusters_original"))

Sycon_blobOnly <- PrepSCTFindMarkers(Sycon_blobOnly)

saveRDS(Sycon_blobOnly, file = "13_recluster_blob/Sycon_blobOnly.Rds")

scCustomize::as.anndata(x = Sycon_blobOnly, main_layer = "counts",
                        other_layers = NULL, file_path = "13_recluster_blob/04_SAMap/",
                        file_name = "Scil_blobOnly.h5ad")


#####################
#     DIM PLOTS     #
#####################

Sycon_blobOnly <- readRDS("13_recluster_blob/Sycon_blobOnly.Rds")

# dimplot with the original cluster names
dimplot_clusters_original_raw <- DimPlot(Sycon_blobOnly, group.by = "seurat_clusters_original")
  
dimplot_clusters_original <- dimplot_clusters_original_raw@data %>%
  mutate(seurat_clusters_original = fct_infreq(as.factor(seurat_clusters_original))) %>%
  mutate(seurat_clusters_original = fct_rev(seurat_clusters_original)) %>%
  ggplot(aes(x = umap_1, y = umap_2,
             colour = seurat_clusters_original, label = seurat_clusters_original)) +
  geom_point(size = 1) +
  
  geom_label(data = . %>%
               group_by(seurat_clusters_original) %>%
               summarise(umap_1 = median(umap_1),
                         umap_2 = median(umap_2),
                         .groups = "drop"),
             aes(label = seurat_clusters_original),
             size = 2.5, linewidth = 0.5,
             text.color = "black", fill = alpha("white", 0.8)) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 23,
                                                      sort = "diff")) + 
  
  labs(title = "Central blob (original clusters)",
       x = "UMAP 1", y = "UMAP 2",
       color = "Original clusters") +
  
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                        y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                        arrow = arrow(length = unit(2.5, "mm"),
                                                      ends = "last", type = "open"),
                                        gp = gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                        y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                        arrow = arrow(length = unit(2.5, "mm"),
                                                      ends = "last", type = "open"),
                                        gp = gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = textGrob(label = "UMAP 1",
                                    x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                    just = c(0, 1), gp = gpar(fontsize = 10))) +
  
  annotation_custom(grob = textGrob(label = "UMAP 2",
                                    x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                    just = c(0, 0), rot = 90, gp = gpar(fontsize = 10))) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 10, 8, "mm")) +
  theme_for_UMAPS
dimplot_clusters_original
  
# dimplot with the re-clustered cluster names
dimplot_clusters_new_raw <- DimPlot(Sycon_blobOnly, group.by = "seurat_clusters_new")

dimplot_clusters_new <- dimplot_clusters_new_raw@data %>%
  mutate(seurat_clusters_new = fct_infreq(as.factor(seurat_clusters_new))) %>%
  mutate(seurat_clusters_new = fct_rev(seurat_clusters_new)) %>%
  ggplot(aes(x = umap_1, y = umap_2,
             colour = seurat_clusters_new, label = seurat_clusters_new)) +
  geom_point(size = 1) +
  
  geom_label(data = . %>%
               group_by(seurat_clusters_new) %>%
               summarise(umap_1 = median(umap_1),
                         umap_2 = median(umap_2),
                         .groups = "drop"),
             aes(label = seurat_clusters_new),
             size = 2.5, linewidth = 0.5,
             text.color = "black", fill = alpha("white", 0.8)) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 17,
                                                      sort = "diff")) + 
  
  labs(title = "Central blob (reclustered)", x = "UMAP 1", y = "UMAP 2",
       color = "New clusters") +
  
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                        y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                        arrow = arrow(length = unit(2.5, "mm"),
                                                      ends = "last", type = "open"),
                                        gp = gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                        y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                        arrow = arrow(length = unit(2.5, "mm"),
                                                      ends = "last", type = "open"),
                                        gp = gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = textGrob(label = "UMAP 1",
                                    x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                    just = c(0, 1), gp = gpar(fontsize = 10))) +
  
  annotation_custom(grob = textGrob(label = "UMAP 2",
                                    x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                    just = c(0, 0), rot = 90, gp = gpar(fontsize = 10))) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 10, 8, "mm")) +
  theme_for_UMAPS
dimplot_clusters_new

ggsave(filename = "13_recluster_blob/dimplot_blobOnly_newClusters.pdf",
       dimplot_clusters_new, device = cairo_pdf,
       dpi = 300, height = 8/1.5, width = 8/1.5, units = ("in"), bg = 'white')
ggsave(filename = "13_recluster_blob/dimplot_blobOnly_newClusters.png",
       dimplot_clusters_new, device = "png",
       dpi = 300, height = 8/1.5, width = 8/1.5, units = ("in"), bg = 'white')

# get the panel
panel_clusters <- ggpubr::ggarrange(dimplot_clusters_original, dimplot_clusters_new,
                  ncol = 2)

panel_clusters

ggsave("13_recluster_blob/dimplot_blobOnly_clusters_panel.pdf",
       panel_clusters, device = cairo_pdf,
       dpi = 300, height = 8/1.5, width = 16/1.5, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/dimplot_blobOnly_clusters_panel.png",
       panel_clusters, device = "png",
       dpi = 300, height = 8/1.5, width = 16/1.5, units = ("in"), bg = 'white')


################################
#      CLUSTER COMPOSITION     #
################################

# get the new cluster composition relative to the original ones
new_cluster_composition <- Sycon_blobOnly[[]] %>%
  select(seurat_clusters_new, seurat_clusters_original) %>%
  group_by(seurat_clusters_new, seurat_clusters_original) %>%
  count()

# plot the cluster composition
cluster_composition_seurat_style <- ClusterDistrBar(Sycon_blobOnly$seurat_clusters_new,
                                                    Sycon_blobOnly$seurat_clusters_original,
                                                    percent = FALSE)
cluster_composition_seurat_style

color_palette <- ggplot_build(cluster_composition_seurat_style)$data[[1]]$fill %>%
  unique()

new_cluster_composition_col <- new_cluster_composition %>%
  ggplot(aes(x = n, y = seurat_clusters_new, fill = seurat_clusters_original)) +
  geom_col(position = "stack", colour = "white") +
  
  labs(#title = "New cluster composition",
       fill = "Original clusters") +
  ylab("New clusters") +
  xlab("Number of cells") +
  
  scale_fill_manual(values = c(color_palette)) +
  scale_x_continuous(expand = c(0, 0)) +
  
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

new_cluster_composition_col

# ggsave("13_recluster_blob/new_cluster_composition_barplot.pdf",
#        new_cluster_composition_col, device = cairo_pdf,
#        dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

new_cluster_composition_tile <- table(Sycon_blobOnly$seurat_clusters_original, Sycon_blobOnly$seurat_clusters_new) %>%
  as.data.frame() %>%
  rename(original_clusters = Var1, new_clusters = Var2, count = Freq) %>%
  group_by(new_clusters) %>%
  mutate(normalized = count / sum(count)) %>%  # proportion of cluster
  ungroup() %>%
  filter(count != 0) %>%
  
  ggplot(aes(new_clusters, original_clusters, fill = normalized)) +
  geom_tile() +
  geom_text(aes(label = count), color = "black", size = 3) +
  
  scale_fill_gradient(low = alpha("white", 0), high = "#37d461", limits = c(0,1),
                      breaks = seq(0, 1, 0.25), labels = format(seq(0, 1, 0.25), nsmall = 2)) +
  labs(#title = "New cluster composition",
       fill = "Normalized\nvalues") +
  xlab("New clusters") +
  ylab("Original clusters") +
  coord_flip() +
  
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

new_cluster_composition_tile

# ggsave("13_recluster_blob/new_cluster_composition_tile.pdf",
#        new_cluster_composition_tile, device = cairo_pdf,
#        dpi = 300, height = 7, width = 9, units = ("in"), bg = 'white')

cluster_composition_panel <- ggpubr::ggarrange(new_cluster_composition_col, new_cluster_composition_tile,
                  ncol = 2, align = "hv", widths = c(2,1.5)) %>%
  ggpubr::annotate_figure(top = ggpubr::text_grob("Cluster composition", face = "bold"))

ggsave("13_recluster_blob/new_cluster_composition_panel.pdf",
       cluster_composition_panel, device = cairo_pdf,
       dpi = 300, height = 7, width = 18, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/new_cluster_composition_panel.png",
       cluster_composition_panel, device = "png",
       dpi = 300, height = 7, width = 18, units = ("in"), bg = 'white')


##################################
#     GET MARKERS PER CLUSTER    #
##################################

# all_markers <- Sycon_blobOnly %>%
#   FindAllMarkers(group.by = "seurat_clusters_new")
# 
# write.table(all_markers, "13_recluster_blob/02_onlyBlob_newClusters/all_markers_blobOnly.tsv",
#             quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

markers_blobOnly_originalClusters <- get_markers_and_save(Sycon_blobOnly,
                                                          "seurat_clusters_original",
                                                          "13_recluster_blob/01_onlyBlob_originalClusters/")

markers_blobOnly_newClusters <- get_markers_and_save(Sycon_blobOnly,
                                                     "seurat_clusters_new",
                                                     "13_recluster_blob/02_onlyBlob_newClusters/")

markers_venn <- ggvenn::ggvenn(list("Original cluster\nmarkers" = markers_blobOnly_originalClusters$gene,
                                    "New cluster\nmarkers" = markers_blobOnly_newClusters$gene),
                               auto_scale = TRUE)
markers_venn

ggsave("13_recluster_blob/venn_markers.pdf",
       markers_venn, device = cairo_pdf,
       dpi = 300, height = 6, width = 6, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/venn_markers.png",
       markers_venn, device = "png",
       dpi = 300, height = 6, width = 6, units = ("in"), bg = 'white')


#############################
#     GET GENE UNIVERSE     #
#############################

get_gene_universe(Sycon_blobOnly, "13_recluster_blob/sycon_onlyBlob_geneUniverse.ls")
