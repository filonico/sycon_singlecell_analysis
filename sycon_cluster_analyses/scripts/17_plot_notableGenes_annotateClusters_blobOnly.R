#!/usr/bin/env Rscript

setwd("/lustre/alice3/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(SeuratExtend)
library(DESeq2)
library(tidyverse)


#######################
#     LOAD INPUTS     #
#######################

Sycon_blobOnly <- readRDS("13_recluster_blob/Sycon_blobOnly.Rds")
load("00_input/Sycon_Seuratv4.Rdata")

sycon_blobOnly_allMarkers_default <- FindAllMarkers(Sycon_blobOnly,
                                                    group.by = "seurat_clusters_new",
                                                    only.pos = TRUE)

DimPlot(Sycon_blobOnly, group.by = "seurat_clusters_new")


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

sycon_blobOnly_allMarkers_default %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

heatmap <- DoHeatmap(Sycon_blobOnly, features = top5$gene) + NoLegend() +
  theme(axis.text = element_text(size = 5))
heatmap

top5 %>%
  write.table(file = "13_recluster_blob/06_notableGenes_clusterAnnotation/top5_markers_perCluster.tsv",
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/top5_marker_heatmap.pdf",
       heatmap, device = cairo_pdf,
       height = 8*2, width = 10*2, dpi = 300, bg = "white")
ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/top5_marker_heatmap.png",
       heatmap, device = "png",
       height = 8*2, width = 10*2, dpi = 300, bg = "white")

sycon_blobOnly_allMarkers_default %>%
  filter(p_val_adj < 0.05) %>%
  write.table(file = "13_recluster_blob/06_notableGenes_clusterAnnotation/all_markers_perCluster.tsv",
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
  "Gli" = c("g8905")
)

p <- Sycon_blobOnly %>% SeuratExtend::DotPlot2(features = notable_genes, flip = TRUE)

dotplot_notable_genes <- p$data %>%
  mutate(Var2 = factor(Var2,
                       levels = c("1", "7", "11", "12", "13", "16", "23", "21", "26", "29",
                                  "15", "25",
                                  "32", "27", "24", "19", "17",
                                  "22", "20", "28", "30", "31",
                                  "0", "2", "3", "4", "5", "6", "8", "9", "10", "14", "18"))) %>%
  filter(pct > 0) %>%
  
  semi_join(sycon_blobOnly_allMarkers_default,
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
        strip.text = element_text(face = "bold", size = 6),
        strip.clip = "off",
        axis.text.x = element_text(angle = 45, hjust = 0, size = 5),
        panel.border = element_rect(color = "#4f4f4f", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        panel.grid = element_line(color = "#dbdbdb"))

dotplot_notable_genes

ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes.pdf",
       dotplot_notable_genes, device = cairo_pdf,
       dpi = 300, height = 4*1.2, width = 6*1.2, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes.png",
       dotplot_notable_genes, device = "png",
       dpi = 300, height = 4*1.2, width = 6*1.2, units = ("in"), bg = 'white')

dotplot_notable_genes_sig <- p$data %>%
  mutate(Var2 = factor(Var2,
                       levels = c("1", "7", "11", "12", "13", "16", "23", "21", "26", "29",
                                  "15", "25",
                                  "32", "27", "24", "19", "17",
                                  "22", "20", "28", "30", "31",
                                  "0", "2", "3", "4", "5", "6", "8", "9", "10", "14", "18"))) %>%
  filter(pct > 0) %>%
  semi_join(sycon_blobOnly_allMarkers_default %>%
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

ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.pdf",
       dotplot_notable_genes_sig, device = cairo_pdf,
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')
ggsave("13_recluster_blob/06_notableGenes_clusterAnnotation/dotplot_notableGenes_sig.png",
       dotplot_notable_genes_sig, device = "png",
       dpi = 300, height = 8, width = 15, units = ("in"), bg = 'white')


######################################
#     PLOT UMAPS FOR PUBLICATION     #
######################################

dimplot_clusters_new_raw <- DimPlot(Sycon_blobOnly, group.by = "seurat_clusters_new")

dimplot_clusters_new <- dimplot_clusters_new_raw@data %>%
  mutate(seurat_clusters_new = fct_infreq(as.factor(seurat_clusters_new))) %>%
  mutate(seurat_clusters_new = fct_rev(seurat_clusters_new)) %>%
  ggplot(aes(x = umap_1, y = umap_2,
             colour = seurat_clusters_new, label = seurat_clusters_new)) +
  geom_point(size = 2) +
  
  geom_label(data = . %>%
               group_by(seurat_clusters_new) %>%
               summarise(umap_1 = median(umap_1),
                         umap_2 = median(umap_2),
                         .groups = "drop"),
             aes(label = seurat_clusters_new),
             size = 4.5, linewidth = 0.5,
             text.color = "black", fill = alpha("white", 0.8)) +
  
  scale_color_manual(values = SeuratExtend::color_pro("default", n = 17,
                                                      sort = "diff")) + 
  
  labs(x = "UMAP 1", y = "UMAP 2",
       color = "New clusters") +
  
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
dimplot_clusters_new

inset_umap_clusters_raw <- DimPlot(Sycon, group.by = "seurat_clusters")

inset_umap_clusters <- inset_umap_clusters_raw@data %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = . %>%
               subset(! seurat_clusters %in% c("14","8","31","12","16","23","4","28","7","1","11","2","13","20","3","29","6","18","9","10","5","0","22")),
             color = "grey85", size = 0.3) +
  geom_point(data = . %>%
               subset(seurat_clusters %in% c("14","8","31","12","16","23","4","28","7","1","11","2","13","20","3","29","6","18","9","10","5","0","22")),
             color = "grey40", size = 0.3) +
  annotate("text", x = -6, y = 9,
           label = "i", fontface = "bold") +
  
  theme_bw(base_size = 12) +
  theme_for_UMAPS +
  theme(plot.background = element_blank(), 
        panel.border = element_rect(color = "#4f4f4f"),
        panel.background = element_blank())
inset_umap_clusters

umaps_with_inset <- cowplot::ggdraw() +
  cowplot::draw_plot(inset_umap_clusters, x = 0, y = 0.74, width = 0.22, height = 0.22) +
  cowplot::draw_plot(dimplot_clusters_new)
umaps_with_inset


##############################################
#     PLOT NOTABLE GENES FOR PUBLICATION     #
##############################################

notable_gene_names <- data.frame(gene_name = c("AE-like1","βcatA","βcatB","Bra1","Bra2","CA1","CA2",
                                               "Cdx","Diactinin","DvlA","DvlB","Elav","Eya",
                                               "FzdA","FzdB","FzdC","FzdD","Gata","Gli","Hedgling",
                                               "Hex","Hmx","Lrp5/6","MsiA","MsiB","Msx",
                                               "Nanos","Ncbt-like1","NkA","NkB","PaxB","Pl10A","Pl10B",
                                               "SixA","Smad1/5","Smad4","SmadRa",
                                               "Sox6","Sox7","SoxB","SoxC","SoxE","SoxF1","SoxF2","Spiculin",
                                               "TboxA","TboxB","TboxC","TcfA","TcfB",
                                               "TgfβB","TgfβD","TgfβF","TgfβG","TgfβN","Triactinin","VasaA","VasaB",
                                               "WntA","WntC","WntD","WntF","WntG","WntI","WntJ","WntL","WntN","WntQ","WntR","WntS","WntT"),
                                 gene_ID = c("g8503","g13719","g4975","g3032","g3118","g782","g11209","g8932",
                                             NA,"g949",NA,"g7460","g8356",NA,NA,"g4648","g3042","g2248","g8905",
                                             NA,"g8162","g9066","g13476","g11831","g12083","g9084","g3613","g10087",
                                             "g8884","g6644","g3655","g6724","g12152","g7040","g9893",
                                             "g11243","g9154","g11067","g5646","g5684","g1695","g11976","g8388",NA,
                                             "g7906","g9948",NA,"g10519","g9","g11378","g13016",NA,"g8564","g2075","g6992",
                                             "g9905","g3867","g9458","g7201","g13150","g4399",NA,NA,"g13140",NA,
                                             "g13308","g6996","g13684","g7152","g11232","g1969")) %>% as_tibble()

top5_gene_names <- data.frame(gene_name = c("Dbh-like","g11926","g4666","g12872","g12871","MacPf","g13341","g13347","g13348","Cystatin_1","g13339","g7108","Baf",
                                            "g11393","Cystatin_2","Itln-1","Cystatin_3","vWD","Histone_H1/5","Histone_H1","g13237","Histone_H2A","Tmem14DP","Cystatin_4","EGF/Laminin_1","g7939",
                                            "g7941","Hsp70_1","Hsp70_2","g10860","Ferritin_1","Ferritin_2","g7090","g13184","vWA","g2025","g2120","EGF/Laminin_2","g13311",
                                            "Tsp1","g6738","PK-like","Tmem121","g10798","Pb1","Hsp40","Tmsβ4X","g12714","DD3-3","Emp","Actin","Tsp-1",
                                            "EGF/Laminin_3","g4274","g5695","g4644","Cfap141","Dm10","g13219","g13222","g1429","Collagen","g13220","g12411","g12430",
                                            "g12426","g11803","g12419"),
                              gene_ID = c("g7011","g11926","g4666","g12872","g12871","g6860","g13341","g13347","g13348","g13351","g13339","g7108","g7921",
                                          "g11393","g13349","g12641","g13325","g3841","g10075","g2502","g13237","g5592","g6217","g13344","g10457","g7939",
                                          "g7941","g8652","g9627","g10860","g7027","g7113","g7090","g13184","g11565","g2025","g2120","g11398","g13311",
                                          "g12405","g6738","g3910","g3985","g10798","g10859","g3385","g1131","g12714","g6294","g5159","g13524","g2945",
                                          "g11306","g4274","g5695","g4644","g4754","g9957","g13219","g13222","g1429","g7145","g13220","g12411","g12430",
                                          "g12426","g11803","g12419")) %>% as_tibble()


p <- Sycon_blobOnly %>% SeuratExtend::DotPlot2(features = notable_gene_names$gene_ID)
p_heatmap <- Sycon_blobOnly %>% SeuratExtend::DotPlot2(features = heatmap@data$Feature)

blobOnly_combined_dotplot <- p@data %>%
  mutate(group = "Marker genes") %>%
  add_row(p_heatmap@data %>%
            mutate(group = "Top-5 expressed genes per cluster")) %>%
  
  mutate(Var2 = as_factor(Var2)) %>%
  
  # remove genes with 0 expression in cells
  filter(pct > 0) %>%
  
  # add gene names
  left_join(notable_gene_names %>%
              add_row(top5_gene_names), by = join_by("Var1" == "gene_ID")) %>%
  
  # remove genes that are not significant in any cell
  filter(! gene_name %in% c("SoxC", "SoxF1", "Lrp5/6", "WntD", "WntS", "MsiA", "Gata", "TcfA", "βcatA", "Bra2", "Pl10B")) %>%
  
  # add p-values according to FindAllMarkers
  left_join(sycon_blobOnly_allMarkers_default %>%
              select(c(p_val_adj, cluster, gene)), by = join_by("Var2" == "cluster", "Var1" == "gene"),
            relationship = "many-to-many") %>%
  
  mutate(gene_name = coalesce(gene_name, Var1),
         gene_name = factor(gene_name, levels = c("SoxE","Pl10A","Elav","WntQ","VasaA","VasaB","Smad1/5","Sox6","Sox7","SoxB","SoxC","SoxF1",
                                                  "TgfβB","TgfβF","TgfβG","TgfβN","WntA","WntC","WntD","WntI","WntL","WntN","WntR","WntS",
                                                  "WntT","Eya","SmadRa","MsiB","Gata","FzdC","FzdD","DvlA","TcfB","Cdx","Gli","PaxB",
                                                  "AE-like1","CA1","CA2","NCBT-like1","Spiculin","Triactinin","βcatB","TboxC","NKA","NKB",
                                                  "Nanos","Hex","Smad4","TboxA","Hmx","SixA","Bra1","MsiA","Lrp5/6","TcfA","βcatA","Bra2","Pl10B",
                                                  "Dbh-like","g11926","g4666","g12872","g12871","MacPf","g13341","g13347","g13348","Cystatin_1","g13339","g7108","Baf",
                                                  "g11393","Cystatin_2","Itln-1","Cystatin_3","vWD","Histone_H1/5","Histone_H1","g13237","Histone_H2A","Tmem14DP","Cystatin_4","EGF/Laminin_1","g7939",
                                                  "g7941","Hsp70_1","Hsp70_2","g10860","Ferritin_1","Ferritin_2","g7090","g13184","vWA","g2025","g2120","EGF/Laminin_2","g13311",
                                                  "Tsp1","g6738","PK-like","Tmem121","g10798","Pb1","Hsp40","Tmsβ4X","g12714","DD3-3","Emp","Actin","Tsp-1",
                                                  "EGF/Laminin_3","g4274","g5695","g4644","Cfap141","Dm10","g13219","g13222","g1429","Collagen","g13220","g12411","g12430",
                                                  "g12426","g11803","g12419")),
         sig = case_when(p_val_adj < 0.05 ~ "*",
                         TRUE ~ NA)) %>%
  
  arrange(zscore) %>%
  
  ggplot(aes(gene_name, as_factor(Var2))) +
  
  geom_hline(yintercept = -Inf, color = "#4f4f4f") +
  
  geom_point(aes(size = pct, color = zscore), shape = 19) +
  geom_text(aes(label = sig), size = 3, color = "white", fontface = "bold", vjust = 0.65) +
  
  scale_color_distiller(palette = "Blues", direction = 1) +
  
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1,7)) +
  
  labs(x = "", y = "Cell clusters",
       color = "z-score", size = "Percent expressed") +
  
  facet_wrap(. ~ group, ncol = 2, space = "free_x", scales = "free_x", strip.position = "bottom") +
  
  coord_cartesian(clip = "off") +
  guides(size = guide_legend(label.position = "bottom", override.aes = list(color = "#77b3d6")),
         color = guide_colorbar(barheight = 1.4, barwidth = 5)) +
  
  theme_bw(base_size = 14) +
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
blobOnly_combined_dotplot

panel_fig2 <- ggpubr::ggarrange(umaps_with_inset, blobOnly_combined_dotplot,
                  ncol = 2, labels = "AUTO", widths = c(0.8,2))
panel_fig2

ggsave("13_recluster_blob/fig2_panel.pdf",
       panel_fig2, device = cairo_pdf,
       width = 30/1.5, height = 10/1.5, units = "in", dpi = 300, bg = "white")
ggsave("13_recluster_blob/fig2_panel.png",
       panel_fig2, device = "png",
       width = 30/1.5, height = 10/1.5, units = "in", dpi = 300, bg = "white")
