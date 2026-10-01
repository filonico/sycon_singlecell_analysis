#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(SeuratExtend)
library(tidyverse)


#######################
#     LOAD INPUTS     #
#######################

load("00_input/Sycon_Seuratv4.Rdata")

sycon_allMarkers_default <- FindAllMarkers(Sycon,
                                           only.pos = TRUE)


###############################
#     DEFINE PLOT THEMES      #
###############################

theme_for_UMAPS <- theme(
  # aspect.ratio = 1,
  plot.background = element_blank(), 
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 10, face = "bold"),
  plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)


#######################
#     GET MARKERS     #
#######################

# DimPlot2(Sycon, features = c("g8720", "g6836", "g2683"))

sycon_allMarkers_default %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

heatmap <- DoHeatmap(Sycon, features = top5$gene) + NoLegend() +
  theme(axis.text = element_text(size = 5))
heatmap

DotPlot2(Sycon, features = top5$gene)

top5 %>%
  write.table(file = "07_notableGenes_clusterAnnotation/top5_markers_perCluster.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


ggsave("07_notableGenes_clusterAnnotation/top5_marker_heatmap.pdf",
       heatmap, device = cairo_pdf,
       height = 8*2, width = 10*2, dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/top5_marker_heatmap.png",
       heatmap, device = "png",
       height = 8*2, width = 10*2, dpi = 300, bg = "white")

sycon_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  write.table(file = "07_notableGenes_clusterAnnotation/all_markers_perCluster.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

  

###############################
#     PLOT NOTABLE GENES      #
###############################

notable_genes <- list(
  "Sox" = c("g11067","g5646","g5684","g1695","g11976","g8388"),
  "Tgf-B\nligands" = c("g4396","g13016","g8882","g8564","g2075","g13405","g6992","g2931","g12087","g13365"),
  "Wnt" = c("g7201","g8624","g13150","g4399","g6515","g1564","g13140","g2924","g13308","g123","g6996","g12657","g13684","g7152","g11232","g1969"),
  "Pax/Six/\nEya" = c("g7040","g8356","g3655"),
  "Spicule\nformation" = c("g782","g11209","g10087","g8503","g7906","g9905", "g9130", "g6914", "g5755", "g5244", "g2825"),
  "Fzd\n(Wnt\npath)" = c("g4648","g3042"),
  "Dvl\n(Wnt\npath)" = c("g949"),
  "Lrp\n(Wnt\npath)" = c("g13476"),
  "Tcf\n(Wnt\npath)" = c("g9","g11378"),
  "Beta-cat\n(Wnt path)" = c("g13719","g4975"),
  "Smad" = c("g9893","g11243","g4903","g1511","g9154","g11177"),
  "Bra/Tbox" = c("g3032","g3118","g9948","g10519"),
  "Totipotency/\ngermline/\noocytes" = c("g3867","g9458", "g12083", "g11831", "g3684", "g12115", "g6724","g12152", "g3613"),
  "Elav" = c("g7460"),
  "Gata" = c("g2248"),
  "Gli" = c("g8905"),
  "HD\ngenes" = c("g8932","g8162","g9066","g9084","g8884","g6644","g8884","g6644")
)

p <- Sycon %>% SeuratExtend::DotPlot2(features = notable_genes, flip = TRUE)

dotplot_notable_genes <- p$data %>%
  mutate(Var2 = factor(Var2,
                       levels = c("1", "7", "11", "12", "13", "16", "23", "21", "26", "29",
                                  "15", "25",
                                  "32", "27", "24", "19", "17",
                                  "22", "20", "28", "30", "31",
                                  "0", "2", "3", "4", "5", "6", "8", "9", "10", "14", "18"))) %>%
  filter(pct > 0) %>%
  semi_join(sycon_allMarkers_default,
            by = c("Var1" = "gene", "Var2" = "cluster")) %>%
  # left_join(notable_gene_names, by = join_by("Var1" == "gene_ID")) %>%
  
  ggplot(aes(Var1, Var2)) +
  geom_point(aes(size = pct, fill = zscore), shape = 21) +
  
  scale_fill_viridis_c(option = "plasma") +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(2,9)) +
  
  xlab("Notable genes") +
  ylab("Cell clusters") +
  labs(size = "Percent\nexpressed") +
  
  facet_grid(cols = vars(FeatureGroup), scales = "free_x", space = "free") +
  
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 0, size = 7),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

dotplot_notable_genes

ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes.pdf",
       dotplot_notable_genes, device = cairo_pdf,
       dpi = 300, height = 8, width = 20, units = ("in"), bg = 'white')
ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes.png",
       dotplot_notable_genes, device = "png",
       dpi = 300, height = 8, width = 20, units = ("in"), bg = 'white')

dotplot_notable_genes_sig <- p$data %>%
  mutate(Var2 = factor(Var2,
                       levels = c("1", "7", "11", "12", "13", "16", "23", "21", "26", "29",
                                  "15", "25",
                                  "32", "27", "24", "19", "17",
                                  "22", "20", "28", "30", "31",
                                  "0", "2", "3", "4", "5", "6", "8", "9", "10", "14", "18"))) %>%
  filter(pct > 0) %>%
  semi_join(sycon_allMarkers_default %>%
              filter(p_val_adj >= 0.05),
            by = c("Var1" = "gene", "Var2" = "cluster")) %>%
  
  ggplot(aes(Var1, Var2)) +
  geom_point(aes(size = pct, fill = zscore), shape = 21) +
  
  scale_fill_viridis_c(option = "plasma") +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(2,9)) +
  
  xlab("Notable genes") +
  ylab("Cell clusters") +
  labs(size = "Percent\nexpressed") +
  
  facet_grid(cols = vars(FeatureGroup), scales = "free_x", space = "free") +
  
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 0, size = 7),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

dotplot_notable_genes_sig

ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.pdf",
       dotplot_notable_genes_sig, device = cairo_pdf,
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')
ggsave("07_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.png",
       dotplot_notable_genes_sig, device = "png",
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')


##################################
#     PLOT ANNOTATED_MARKERS     #
##################################

# if a cluster significantly expresses known markers, we give it a name
# if a cluster maps with a spongilla cluster, we append the name 
# if a cluster does not express any known marker, we call it "unknown"
# if a cluster expresses known markers (but not significantly), we call it "uncertain"
# if a cluster significantly expresses known markers but with multiple identities, we call it "uncertain"
cluster_identity <- c("0" = "Uncertain_0", # uncertain identity (choanos + development)
                      "1" = "Choanocytes_1",
                      "2" = "Unknown_2",
                      "3" = "Unknown_3",
                      "4" = "Unknown_4",
                      "5" = "Uncertain_5", # uncertain identity (choanos + development)
                      "6" = "Unknown_6",
                      "7" = "Choanocytes_7",
                      "8" = "Unknown_8",
                      "9" = "Unknown_9",
                      "10" = "Uncertain_10", # uncertain identity (choanos + development)
                      "11" = "Choanocytes_11",
                      "12" = "Unknown_12",
                      "13" = "Unknown_13",
                      "14" = "Unknown_14",
                      "15" = "Metabolocyte-_and_pinacocyte-like_early_embryos",
                      "16" = "Choanocytes_16",
                      "17" = "Oocytes/early_embryos",
                      "18" = "Unknown_18",
                      "19" = "Oocytes/early_embryos",
                      "20" = "Unknown_20",
                      "21" = "Myopeptidocyte-like_choanocytes",
                      "22" = "Unknown_22",
                      "23" = "Choanocytes_23",
                      "24" = "Oocytes/early_embryos",
                      "25" = "Sclerocytes",
                      "26" = "Uncertain_26", # uncertain identity (choanos + development)
                      "27" = "Oocytes/early_embryos",
                      "28" = "Unknown_28",
                      "29" = "Unknown_29",
                      "30" = "Accessory_cells",
                      "31" = "Unknown_31",
                      "32" = "Oocytes/early_embryos")
                      

cluster_identity_df <- as_tibble(cluster_identity) %>%
  rownames_to_column(var = "cluster_ID") %>%
  mutate(cluster_ID = as.character(as.numeric(cluster_ID) - 1))

write.table(cluster_identity_df, "07_notableGenes_clusterAnnotation/cluster_identity.tsv", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

Sycon[[]] <- Sycon[[]] %>%
  left_join(cluster_identity_df,
            by = join_by("seurat_clusters" == "cluster_ID")) %>%
  rename("cluster_identity" = "value")

umap_clusters_raw <- DimPlot(Sycon, group.by = "seurat_clusters")
umap_clusters_raw

umap_clusters <- umap_clusters_raw@data %>%
  mutate(seurat_clusters = fct_infreq(as.factor(seurat_clusters))) %>%
  mutate(seurat_clusters = fct_rev(seurat_clusters)) %>%
  
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters, label = seurat_clusters)) +
  geom_point(size = 1) +
  
  geom_label(data = . %>%
               group_by(seurat_clusters) %>%
               summarise(UMAP_1 = median(UMAP_1),
                         UMAP_2 = median(UMAP_2),
                         .groups = "drop"),
             aes(label = seurat_clusters),
             size = 2.5, linewidth = 0.5,
             text.color = "black", fill = alpha("white", 0.8)) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 33,
                                                      sort = "diff")) + 
  
  labs(x = "UMAP 1", y = "UMAP 2",
       color = "Cell clusters") +
  
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"),
                                                            ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"),
                                                            ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))) +
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))) +
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 10, 8, "mm")) +
  theme_for_UMAPS
umap_clusters

umap_annotated_raw <- DimPlot(Sycon, group.by = "cluster_identity")

umap_annotated <- umap_annotated_raw@data %>%
  mutate(cluster_identity = str_replace_all(cluster_identity, "_", " "),
         cluster_identity = str_replace_all(cluster_identity, " [0-9]+$", ""),
         cluster_identity = str_replace(cluster_identity, "-like ", "-like\n"),
         cluster_identity = str_replace(cluster_identity, "/", "/\n"),
         cluster_identity = str_replace(cluster_identity, "and ", "and \n")) %>%
  
  ggplot(aes(UMAP_1, UMAP_2)) +
  
  geom_point(data = . %>%
               subset(grepl("Unknown", cluster_identity)),
             aes(color = "Unknown"), size = 1) +
  geom_point(data = . %>%
               subset(grepl("Uncertain", cluster_identity)),
             aes(colour = "Uncertain"), size = 1) +
  geom_point(data = . %>%
               subset(!grepl("Unknown|Uncertain", cluster_identity)),
             aes(colour = as.factor(cluster_identity)), size = 1) +
  
  scale_color_manual(values = c(SeuratExtend::color_iwh("default", n = 7)[c(1,seq(3,7))],
                                "grey60", "grey85")) + 
  
  guides(color = guide_legend(override.aes = list(size = 2),
                              byrow = TRUE, position = "right",
                              keyheight = 0.8, default.unit = "cm")) +
  
  labs(colour = "Main cell types") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS
umap_annotated

panel_umaps <- ggpubr::ggarrange(umap_clusters, umap_annotated, ncol = 2,
                           widths = c(0.8, 1))
panel_umaps

ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.pdf",
       panel_umaps, device = cairo_pdf,
       width = 18/1.5, height = 8/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/cluster_annotation_umaps.png",
       panel_umaps, device = "png",
       width = 18/1.5, height = 8/1.5, units = "in", dpi = 300, bg = "white")


##############################################
#     PLOT NOTABLE GENES FOR PUBLICATION     #
##############################################

notable_gene_names <- data.frame(gene_name = c("AE-like-1","βcat-A","βcat-B","Bra-1","Bra-2","CA-1","CA-2",
                                               "Cdx","Diactinin","Dvl-A","Dvl-B","Elav","Eya",
                                               "Fzd-A","Fzd-B","Fzd-C","Fzd-D","Gata","Gli","Hedgling",
                                               "Hex","Hmx","Lrp-5/6","Msi-A","Msi-B","Msx",
                                               "Nanos","NCBT-like-1","NK-A","NK-B","Pax-B","Pl10-A","Pl10-B",
                                               "Six-A","Smad-1/5","Smad-4","Smad-Ra",
                                               "Sox-6","Sox-7","Sox-B","Sox-C","Sox-E","Sox-F1","Sox-F2","Spiculin",
                                               "Tbox-A","Tbox-B","Tbox-C","Tcf-A","Tcf-B",
                                               "Tgfβ-B","Tgfβ-D","Tgfβ-F","Tgfβ-G","Tgfβ-N","Triactinin","Vasa-A","Vasa-B",
                                               "Wnt-A","Wnt-C","Wnt-D","Wnt-F","Wnt-G","Wnt-I","Wnt-J","Wnt-L","Wnt-N","Wnt-Q","Wnt-R","Wnt-S","Wnt-T"),
                                 gene_ID = c("g8503","g13719","g4975","g3032","g3118","g782","g11209","g8932",
                                             NA,"g949",NA,"g7460","g8356",NA,NA,"g4648","g3042","g2248","g8905",
                                             NA,"g8162","g9066","g13476","g11831","g12083","g9084","g3613","g10087",
                                             "g8884","g6644","g3655","g6724","g12152","g7040","g9893",
                                             "g11243","g9154","g11067","g5646","g5684","g1695","g11976","g8388",NA,
                                             "g7906","g9948",NA,"g10519","g9","g11378","g13016",NA,"g8564","g2075","g6992",
                                             "g9905","g3867","g9458","g7201","g13150","g4399",NA,NA,"g13140",NA,
                                             "g13308","g6996","g13684","g7152","g11232","g1969")) %>% as_tibble()


p <- Sycon %>% SeuratExtend::DotPlot2(features = notable_gene_names$gene_ID)

dotplot_markers_publication <- p$data %>%
  
  # sort y axis values to match cell types
  mutate(Var2 = factor(Var2, levels = c(
    "1","7","11","16","23",
    "21",
    "30",
    "25",
    "17","19","24","27", "32",
    "15",
    "0","5","10","26",
    "2","3","4","6","8","9","12","13","14","18","20","22","28","29","31"
  ))) %>%
  
  # remove genes with 0 expression in cells
  filter(pct > 0) %>%
  
  # add gene names
  left_join(notable_gene_names, by = join_by("Var1" == "gene_ID")) %>%
  filter(!is.na(gene_name)) %>%
  
  # remove genes that are not significant in any cell
  filter(! gene_name %in% c("Sox-C", "Sox-F1", "Lrp-5/6", "Wnt-D", "Wnt-S", "Msi-A", "Gata", "Tcf-A", "βcat-A", "Bra-2", "Pl10-B")) %>%
  
  # add cell cluster names
  left_join(cluster_identity_df, by = join_by("group" == "cluster_ID")) %>%
  mutate(value = str_replace_all(value, "_", " "),
         value = str_replace_all(value, " [0-9]+$", ""),
         value = str_replace(value, "-like ", "-like\n"),
         value = str_replace(value, "/", "/\n"),
         value = str_replace(value, "and ", "and\n"),
         value = factor(value, levels = c("Choanocytes", "Myopeptidocyte-like\nchoanocytes",
                                          "Accessory cells", "Sclerocytes",
                                          "Oocytes/\nearly embryos", "Metabolocyte- and\npinacocyte-like\nearly embryos",
                                          "Uncertain", "Unknown"))) %>%
  rename("cell_identity" = "value") %>%
  
  # add p-values according to FindAllMarkers
  left_join(sycon_allMarkers_default %>%
              select(c(p_val_adj, cluster,gene)), by = join_by("group" == "cluster", "Var1" == "gene"),
            relationship = "many-to-many") %>%

  mutate(gene_name = factor(gene_name, levels = c("Sox-E","Pl10-A","Elav","Wnt-Q","Vasa-A","Vasa-B","Smad-1/5","Sox-6","Sox-7","Sox-B","Sox-C","Sox-F1",
                                                  "Tgfβ-B","Tgfβ-F","Tgfβ-G","Tgfβ-N","Wnt-A","Wnt-C","Wnt-D","Wnt-I","Wnt-L","Wnt-N","Wnt-R","Wnt-S",
                                                  "Wnt-T","Eya","Smad-Ra","Msi-B","Gata","Fzd-C","Fzd-D","Dvl-A","Tcf-B","Cdx","Gli","Pax-B",
                                                  "AE-like-1","CA-1","CA-2","NCBT-like-1","Spiculin","Triactinin","βcat-B","Tbox-C","NK-A","NK-B",
                                                  "Nanos","Hex","Smad-4","Tbox-A","Hmx","Six-A","Bra-1","Msi-A","Lrp-5/6","Tcf-A","βcat-A","Bra-2",
                                                  "Pl10B")),
         sig = case_when(p_val_adj < 0.05 ~ "*",
                         TRUE ~ NA)
         ) %>%
  arrange(zscore) %>%
  
  
  ggplot(aes(gene_name, as_factor(Var2))) +
  
  geom_rect(data = . %>% distinct(cell_identity),
            aes(fill = cell_identity), alpha = 0.15, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = Inf, color = "#4f4f4f") +
  scale_fill_manual(values = c("Accessory cells" = "#B49440",
                               "Archaeocyte-like\nstem cells" = "#B65CBF",
                               "Choanocytes" = "#67A64E",
                               "Metabolocyte- and\npinacocyte-like\nearly embryos" = "#747CC9",
                               "Myopeptidocyte-like\nchoanocytes" = "#CA5F3E",
                               "Oocytes/\nearly embryos" = "#48B1A7",
                               "Sclerocyte-like\npinacocytes" = "#C8577B",
                               "Uncertain" = "grey60",
                               "Unknown" = NA),
                    guide = "none") +
  ggnewscale::new_scale_fill() +
  
  geom_point(aes(size = pct, fill = zscore), shape = 21, color = "#08306B") +
  geom_text(aes(label = sig), size = 3, color = "white", fontface = "bold", vjust = 0.65) +
  
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_color_gradient(high = "#08306B", low = "#DEEBF7") +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1,7)) +
  
  labs(x = "Marker genes", y = "Cell clusters",
       fill = "z-score", size = "Percent expressed") +
  
  facet_wrap(. ~ cell_identity, space = "free_y", scales = "free_y", strip.position = "right") +
  
  coord_cartesian(clip = "off") +
  guides(size = guide_legend(label.position = "bottom", override.aes = list(color = "#08306B")),
         color = guide_colorbar(barheight = 1.4)) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA), 
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", lineend = "round", linewidth = .5),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(.4, "cm"),
    legend.key.height = unit(.4, "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size = 8),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "grey90", linewidth = .5),
    axis.ticks.length = unit(0.10, "cm"),
    axis.text.x = element_text(color = "black", angle = 45, face = "italic", hjust = 0,
                               margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(color = "black",
                               margin = margin(t = 0, r = 4, b = 0, l = 0)),
    axis.title.y = element_text(angle = 90, size = 10, face = "bold",
                                margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(angle = 0, size = 10, face = "bold",
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", size = 8, angle = 0, hjust = 0),
    strip.clip = "off",
  )
dotplot_markers_publication

panel_umaps_noLegend <- ggpubr::ggarrange(umap_clusters,
                                          umap_annotated +
                                            ggrepel::geom_text_repel(data = data.frame(x = -2, y = -5.5, label = "Uncertain"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = -0.5, nudge_y = -2, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 1, y = -3.4, label = "Unknown"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = 2.5, nudge_y = -2, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 7.3, y = -10.4, label = "Myopeptidocyte-like\nchoanocytes"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = -1, nudge_y = 2.5, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 9.4, y = -3.1, label = "Metabolocyte- and\npinacocyte-like\nearly-embryos"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = -4, nudge_y = 0.2, max.overlaps = Inf, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 11.5, y = 0, label = "Sclerocyte-like\npinacocytes"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = -4, nudge_y = 0.1, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 10.1, y = 1.9, label = "Accessory cells"),
                                                                     aes(x = x, y = y, label = label), nudge_x = -1.8, nudge_y = 1.2) +
                                            # ggrepel::geom_text_repel(data = data.frame(x = 11.9, y = 5.2, label = "Archaeocyte-like\n stem cells"),
                                            #                          aes(x = x, y = y, label = label),
                                            #                          nudge_x = -4, nudge_y = 0.4, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 9.6, y = 9.6, label = "Oocytes/early embryos"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = -3.5, nudge_y = -2, lineheight = 0.9) +
                                            ggrepel::geom_text_repel(data = data.frame(x = 0.9, y = 4.1, label = "Choanocytes"),
                                                                     aes(x = x, y = y, label = label),
                                                                     nudge_x = 2.5, nudge_y = 1.4, lineheight = 0.9) +
                                            NoLegend(), ncol = 2)
panel_umaps_noLegend

panel_fig1 <- ggpubr::ggarrange(panel_umaps_noLegend, dotplot_markers_publication,
                                nrow = 2, labels = "AUTO", heights = c(1,2))
panel_fig1

ggsave("07_notableGenes_clusterAnnotation/fig1_panel.pdf",
       panel_fig1, device = cairo_pdf,
       width = 20/1.5, height = 22/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/fig1_panel.png",
       panel_fig1, device = "png",
       width = 20/1.5, height = 22/1.5, units = "in", dpi = 300, bg = "white")


#############################################
#     PLOT TOP5 HEATMAP FOR PUBLICATION     #
#############################################

top5_gene_names <- data.frame(gene_name = c("Gs","g9227","Ntf2","g4161","g5797","Fibrinogen-like_1","Dscam-like","Hmx","g13838","Reeler_1","EGF/Laminin_1","g6974","vWA_1",
                                            "vWA_2","Ser_protease_inhibitor_1","Collagen","EGF/Laminin_2","g10751","g7467","g7466","Fibrinogen-like_2","Atp5pb","g12902","Profilin","g10750","Fibrinogen-like_3",
                                            "Litaf","g13598","g7090","g10478","g9945","Actin_1","Vatg","g972","Cox6B1","Trx","Granulin","Grx","g7110","g7109",
                                            "Efh-like_1","g9244","g4644","g5695","g4274","g7689","Actin_2","g11077","Ser_protease_inhibitor_2","Adf_1","g7145","g13220","g11260","g13222",
                                            "g13219","g9786","Rho_GTPase","PK","Efh-like_2","HMG-box","g11803","Ubiquitin-like","g12430","g12426","g12411","Adf_2","Ferritin_1","Actin_3",
                                            "Rpp1C","Thy-β","g5815","Reeler_2","g258","g7014","g3337","g2945","g13794","Cbr/Clec-78","EGF/Laminin_3","Lrrk2","g4687","g11675","g7941","g7939",
                                            "Hsp40/DnaJ","Hsp70_1","Hsp70_2","g10859","g10860","g10798","g11387","Banf2","g6738","DNA-polI","Limr","g13350","g6099","Tsp1","g13311","Clec3a",
                                            "g2120","g2025","g13339","g13346","g13347","g13351","g13184","Col6α3","Spondin","Ferritin_2","Ferritin_3","Armc1","Histone_H2A","g6217","g13237",
                                            "g10075","g7921","Ependymin_1","Fas1","Ependymin_2","g7108","g2122","g1771","g7688","g7949","Cda9","Cystatin_1","Cystatin_2","g13348","g13341","Cystatin_3","g6860",
                                            "g12872","g12871","g4666","g11926","g4241","g5159","g7075","g13456","g4172"),
                              gene_ID = c("g2318","g9227","g849","g4161","g5797","g12641","g4153","g9066","g13838","g609","g7679","g6974","g7611",
                                          "g6674","g11550","g3083","g7469","g10751","g7467","g7466","g7468","g2894","g12902","g13119","g10750","g7428",
                                          "g6842","g13598","g7090","g10478","g9945","g6719","g5874","g972","g838","g4517","g11545","g307","g7110","g7109",
                                          "g9957","g9244","g4644","g5695","g4274","g7689","g4792","g11077","g11549","g10508","g7145","g13220","g11260","g13222",
                                          "g13219","g9786","g7883","g2323","g1061","g10741","g11803","g3232","g12430","g12426","g12411","g1554","g7031","g4223",
                                          "g9721","g5447","g5815","g5618","g258","g7014","g3337","g2945","g13794","g2336","g10457","g11295","g4687","g11675","g7941","g7939",
                                          "g3385","g8652","g9627","g10859","g10860","g10798","g11387","g7923","g6738","g3745","g13498","g13350","g6099","g12405","g13311","g11398",
                                          "g2120","g2025","g13339","g13346","g13347","g13351","g13184","g11565","g5471","g7113","g7027","g10528","g5592","g6217","g13237",
                                          "g10075","g7921","g8206","g1039","g8207","g7108","g2122","g1771","g7688","g7949","g662","g7363","g13344","g13348","g13341","g13349","g6860",
                                          "g12872","g12871","g4666","g11926","g4241","g5159","g7075","g13456","g4172")) %>% as_tibble()

p_heatmap <- Sycon %>% SeuratExtend::DotPlot2(features = heatmap@data$Feature)

top5markers_dotplot <- p_heatmap@data %>%
  mutate(group = "Top-5 expressed genes per cluster") %>%
  mutate(Var2 = as_factor(Var2)) %>%
  
  # remove genes with 0 expression in cells
  filter(pct > 0) %>%
  
  # add gene names
  left_join(top5_gene_names, by = join_by("Var1" == "gene_ID")) %>%
  
  left_join(sycon_allMarkers_default %>%
              select(c(p_val_adj, cluster, gene)), by = join_by("Var2" == "cluster", "Var1" == "gene"),
            relationship = "many-to-many") %>%
  
  mutate(Var1 = factor(Var1, levels = levels(heatmap@data$Feature)),
         gene_name = coalesce(gene_name, Var1),
         gene_name = factor(gene_name, levels = c("Gs","g9227","Ntf2","g4161","g5797","Fibrinogen-like_1","Dscam-like","Hmx","g13838","Reeler_1","EGF/Laminin_1","g6974","vWA_1",
                                                  "vWA_2","Ser_protease_inhibitor_1","Collagen","EGF/Laminin_2","g10751","g7467","g7466","Fibrinogen-like_2","Atp5pb","g12902","Profilin","g10750","Fibrinogen-like_3",
                                                  "Litaf","g13598","g7090","g10478","g9945","Actin_1","Vatg","g972","Cox6B1","Trx","Granulin","Grx","g7110","g7109",
                                                  "Efh-like_1","g9244","g4644","g5695","g4274","g7689","Actin_2","g11077","Ser_protease_inhibitor_2","Adf_1","g7145","g13220","g11260","g13222",
                                                  "g13219","g9786","Rho_GTPase","PK","Efh-like_2","HMG-box","g11803","Ubiquitin-like","g12430","g12426","g12411","Adf_2","Ferritin_1","Actin_3",
                                                  "Rpp1C","Thy-β","g5815","Reeler_2","g258","g7014","g3337","g2945","g13794","Cbr/Clec-78","EGF/Laminin_3","Lrrk2","g4687","g11675","g7941","g7939",
                                                  "Hsp40/DnaJ","Hsp70_1","Hsp70_2","g10859","g10860","g10798","g11387","Banf2","g6738","DNA-polI","Limr","g13350","g6099","Tsp1","g13311","Clec3a",
                                                  "g2120","g2025","g13339","g13346","g13347","g13351","g13184","Col6α3","Spondin","Ferritin_2","Ferritin_3","Armc1","Histone_H2A","g6217","g13237",
                                                  "g10075","g7921","Ependymin_1","Fas1","Ependymin_2","g7108","g2122","g1771","g7688","g7949","Cda9","Cystatin_1","Cystatin_2","g13348","g13341","Cystatin_3","g6860",
                                                  "g12872","g12871","g4666","g11926","g4241","g5159","g7075","g13456","g4172")),
         sig = case_when(p_val_adj < 0.05 ~ "*",
                         TRUE ~ NA)) %>%

  arrange(zscore) %>%
  
  ggplot(aes(gene_name, as_factor(Var2))) +
  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = NA, color = "#4f4f4f", linewidth = 0.4) +
  
  # annotate("rect", xmin = seq(0,160,10), xmax = seq(5,165,10), ymin = -Inf, ymax = Inf,
  #          fill = "grey60", color = NA, alpha = 0.2) +
  
  geom_point(aes(size = pct, color = zscore), shape = 19) +
  geom_text(aes(label = sig), size = 3, color = "white", fontface = "bold", vjust = 0.65) +
  
  scale_color_distiller(palette = "Blues", direction = 1) +
  
  scale_y_discrete(limits = rev, position = "left", sec.axis = sec_axis(transform = identity)) +
  scale_x_discrete(limits = rev, position = "top") +
  scale_size_continuous(range = c(1,7)) +
  
  labs(x = "Top-5 markers per cluster", y = "Cell clusters",
       color = "z-score", size = "Percent expressed") +
  
  coord_cartesian(clip = "off") +
  guides(size = guide_legend(label.position = "bottom", override.aes = list(color = "#77b3d6")),
         color = guide_colorbar(barheight = 1.4, barwidth = 5)) +
  
  theme_bw(base_size = 14) +
  theme(
    plot.margin = margin(0, 5, 0, 5, "mm"),
    plot.background = element_rect(fill = "transparent", colour = NA), 
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", lineend = "round", linewidth = .5),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(.4, "cm"),
    legend.key.height = unit(.4, "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size = 10),
    axis.line = element_blank(),
    axis.ticks.x = element_line(colour = "grey90", linewidth = .5),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.10, "cm"),
    axis.text.x = element_text(color = "black", angle = 45, face = "italic", hjust = 0, size = 7,
                               margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(color = "black",
                               margin = margin(t = 0, r = 4, b = 0, l = 0)),
    axis.title.y = element_text(angle = 90, size = 12, face = "bold",
                                margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(angle = 0, size = 12, face = "bold",
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 12, angle = 0),
    strip.clip = "off"
  )
top5markers_dotplot

ggsave("07_notableGenes_clusterAnnotation/top5_marker_dotplot.pdf",
       top5markers_dotplot, device = cairo_pdf,
       width = 31/1.5, height = 12/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/top5_marker_dotplot.png",
       top5markers_dotplot, device = "png",
       width = 31/1.5, height = 12/1.5, units = "in", dpi = 300, bg = "white")
