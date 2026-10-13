#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)


######################
#     LOAD INPUT     #
######################

load("00_input/Sycon_Seuratv4.Rdata")

sycon_allMarkers_default <- FindAllMarkers(Sycon,
                                           only.pos = TRUE)

sycon_allMarkers_default %>%
  filter(cluster %in% c(1, 7, 11, 16, 23)) %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  select(gene) %>%
  distinct() %>%
  pull(gene)

choanocyte_markers <- c("g7939", "g13347", "g13351", "g13339", "g10860", "g11926", "g10859", "g4666", "g12872", "g3385", "g12871", "g6860", "g307", "g9627", "g8652")

SeuratExtend::DotPlot2(Sycon, features = choanocyte_markers)

SeuratExtend::DimPlot2(Sycon, features = choanocyte_markers)


###############################
#     DEFINE PLOT THEMES      #
###############################

umap_arrows <- list(
  
  annotation_custom(grob = grid::segmentsGrob(
    x0 = unit(0, "mm"), x1 = unit(12, "mm"),
    y0 = unit(0, "mm"), y1 = unit(0, "mm"),
    arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
    gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  
  annotation_custom(grob = grid::segmentsGrob(
    x0 = unit(0, "mm"), x1 = unit(0, "mm"),
    y0 = unit(0, "mm"), y1 = unit(12, "mm"),
    arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
    gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  
  annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  
  coord_cartesian(clip = "off"))

umap_coord_x <- list(annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                                             x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                                             just = c(0, 1), gp = grid::gpar(fontsize = 10))),
                     coord_cartesian(clip = "off"))

umap_coord_y <- list(annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                                             x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                                             just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
                     coord_cartesian(clip = "off"))

theme_for_UMAPS <- theme(
  # aspect.ratio = 1,
  plot.background = element_blank(), 
  panel.border = element_rect(color = "#545454"),
  panel.grid = element_blank(),
  legend.position = "none",
  plot.title = element_text(size = 13, hjust = 0, vjust = 1.75, face = "bold.italic"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_text(color = "#545454"),
  axis.title = element_blank()
)


#####################
#     FUNCTIONS     #
#####################

plot_dimplot_gene <- function(gene_id, gene_name) {
  
  dimplot_raw <- Sycon %>%
    Seurat::FeaturePlot(features = gene_id)
  
  dimplot_def <- dimplot_raw@data %>%
    arrange(.data[[gene_id]]) %>%
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_vline(xintercept = seq(-5, 10, 5), col = "grey95") +
    geom_hline(yintercept = seq(-10, 10, 5), col = "grey95") +
    geom_point(aes(col = .data[[gene_id]]), size = 1) +
    
    # scale_color_distiller(palette = "Blues", direction = 1) +
    # scale_color_gradient(high = "#1c57d5", low = "#e5ecfb") +
    scale_color_gradient(high = "#d51c3b", low = "#e5ecfb") +
    
    scale_y_continuous(breaks = seq(-10, 10, 5)) +
    
    annotate("text", x = -7.3, y = 9.9, hjust = 0,
             label = paste0("italic('", gene_name, "')~'(", gene_id, ")'"),
             parse = TRUE) +
    
    # labs(title = paste0(gene_name, " (", gene_id, ")")) +
    
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme_for_UMAPS
  
  return(dimplot_def)
}


#################################
#     PLOT HCR MARKERS UMAP     #
#################################

hcr_markers <- c(
  "Sox-E" = "g11976",
  "Mpeg-1" = "g4312",
  "Laminin" = "g7679",
  "Spiculin" = "g7906",
  "Granulin" = "g11545",
  "Wnt-N" = "g6996"
  # "Vasa-B" = "g9458",
  # "Elav" = "g7460",
  # "Calcineurin" = "g2493",
  # "Actin" = "g13004",
  # "TgfB-F" = "g8564",
  # "Reeler" = "g5618",
  # "ApoL" = "g7154",
  # "Triactinin" = "g9905",
  # "Tbox-A" = "g9948",
  )


dimplots_hcr_markers <- Map(plot_dimplot_gene, unname(hcr_markers), names(hcr_markers))

dimplots_hcr_markers$g11976 <- dimplots_hcr_markers$g11976 +
  scale_x_continuous(breaks = seq(0, 10, 5)) +
  scale_y_continuous(breaks = seq(-5, 5, 5)) +
  umap_coord_x + umap_coord_y

panel <- patchwork::wrap_plots(dimplots_hcr_markers, nrow = 1) +
  # patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.location = "plot")
panel

ggsave("07_notableGenes_clusterAnnotation/fig2_hcrMarkers_new.png",
       panel, device = "png",
       height = 4.5/1.5, width = 25/1.5, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/fig2_hcrMarkers_new.pdf",
       panel, device = cairo_pdf,
       height = 4.5/1.5, width = 25/1.5, units = "in", dpi = 300, bg = "white")


####################################
#     PLOT HCR MARKERS DOTPLOT     #
####################################

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
                      "15" = "Developing_embryos",
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

hcr_markers <- c(
  "Sox-E" = "g11976",
  "Mpeg-1" = "g4312",
  "Laminin" = "g7679",
  "Spiculin" = "g7906",
  "Granulin" = "g11545",
  "Wnt-N" = "g6996")

hcr_markers_df <- data.frame(gene_name = names(hcr_markers),
                             gene_ID = unname(hcr_markers))

p <- Sycon %>% SeuratExtend::DotPlot2(features = hcr_markers_df$gene_ID)
p

hcr_markers_dotplot <- p$data %>%
  filter(Var2 %in% c("1","7","11","16","23","17","19","24","27","32","15","25","21","30")) %>%
  
  # sort y axis values to match cell types
  mutate(Var2 = factor(Var2, levels = c(
    "1","7","11","16","23",
    "17","19","24","27",
    "32",
    "15",
    "25",
    "21",
    "30"
    # "0","5","10","26",
    # "2","3","4","6","8","9","12","13","14","18","20","22","28","29","31"
  ))) %>%
  
  # remove genes with 0 expression in cells
  filter(pct > 0) %>%
  
  # add gene names
  left_join(hcr_markers_df, by = join_by("Var1" == "gene_ID")) %>%
  
  # add cell cluster names
  left_join(cluster_identity_df, by = join_by("group" == "cluster_ID")) %>%
  mutate(value = str_replace_all(value, "_", " "),
         value = str_replace_all(value, " [0-9]+$", ""),
         value = str_replace(value, "-like ", "-like\n"),
         value = str_replace(value, "/", "/\n"),
         value = str_replace(value, "and ", "and\n"),
         value = str_replace(value, "Developing ", "Developing\n"),
         value = factor(value, levels = c("Choanocytes", "Myopeptidocyte-like\nchoanocytes",
                                          "Accessory cells", "Sclerocytes",
                                          "Oocytes/\nearly embryos", "Developing\nembryos"
                                          ))) %>%
  rename("cell_identity" = "value") %>%
  
  mutate(gene_name = factor(gene_name, levels = c("Sox-E", "Mpeg-1", "Laminin",
                                                  "Spiculin", "Granulin", "Wnt-N"))) %>%
  arrange(zscore) %>%
  
  
  ggplot(aes(as_factor(Var2), gene_name)) +
  
  geom_rect(data = . %>% distinct(cell_identity),
            aes(fill = cell_identity), alpha = 0.15, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_hline(yintercept = -Inf, color = "#4f4f4f") +
  scale_fill_manual(values = c("Accessory cells" = "#B49440",
                               "Choanocytes" = "#67A64E",
                               "Developing\nembryos" = "#747CC9",
                               "Myopeptidocyte-like\nchoanocytes" = "#CA5F3E",
                               "Oocytes/\nearly embryos" = "#48B1A7",
                               "Sclerocytes" = "#C8577B"),
                    guide = "none") +
  ggnewscale::new_scale_fill() +
  
  geom_point(aes(size = pct, fill = zscore), shape = 21, color = "#08306B") +
  
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_color_gradient(high = "#08306B", low = "#DEEBF7") +
  
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(range = c(5,10)) +
  
  labs(y = "Marker genes", x = "Cell clusters",
       fill = "z-score", size = "Percent\nexpressed") +
  
  facet_wrap(. ~ cell_identity, space = "free_x", scales = "free_x", strip.position = "bottom") +
  
  coord_cartesian(clip = "off") +
  guides(size = guide_legend(label.position = "right", override.aes = list(color = "#08306B")),
         color = guide_colorbar(barheight = 1.4)) +
  
  theme_bw(base_size = 18) +
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
    legend.position = "right",
    legend.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size = 8),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "grey90", linewidth = .5),
    axis.ticks.length = unit(0.10, "cm"),
    axis.text.x = element_text(color = "black",
                               margin = margin(t = 4, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(color = "black", face = "italic",
                               margin = margin(t = 0, r = 4, b = 0, l = 0)),
    axis.title.y = element_text(angle = 90, face = "bold",
                                margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(angle = 0, face = "bold",
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", size = 10, angle = 0, hjust = 0),
    strip.clip = "off",
  )
hcr_markers_dotplot

ggsave("07_notableGenes_clusterAnnotation/fig2_hcrMarkers_dotplot.png",
       hcr_markers_dotplot, device = "png",
       height = 2.0*2, width = 6.334*2, units = "in", dpi = 300, bg = "white")
ggsave("07_notableGenes_clusterAnnotation/fig2_hcrMarkers_dotplot.png",
       hcr_markers_dotplot, device = "png",
       height = 2.0*2, width = 6.334*2, units = "in", dpi = 300, bg = "white")
