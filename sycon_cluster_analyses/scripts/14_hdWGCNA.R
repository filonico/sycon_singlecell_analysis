#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

suppressMessages({library(tidyverse)
  library(hdWGCNA)
  library(SeuratExtend)})
  # library(smacof)})


######################
#     LOAD INPUT     #
######################

# load("00_input/Sycon_Seuratv4.Rdata")
# 
# Sycon <- UpdateSeuratObject(Sycon)

Sycon <- readRDS("12_hdWGCNA/Sycon_hdWGCNA.Rds")


#####################
#     FUNCTIONS     #
#####################

# plot UMAP of selected modules
plot_module_UMAP <- function(module_name, color_name) {
  
  umap_df <- Sycon@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column() %>%
    left_join(Sycon@meta.data %>% rownames_to_column()) %>%
    select(rowname, UMAP_1, UMAP_2, all_of(module_name)) %>%
    pivot_longer(-c(rowname, UMAP_1, UMAP_2),
                 names_to = "module",
                 values_to = "expr")
  
  # enforce symmetric range
  max_abs <- max(abs(umap_df$expr), na.rm = TRUE)
  
  umap_df <- umap_df %>%
    mutate(expr = pmax(pmin(expr, max_abs), -max_abs)) %>%
    arrange(expr)
  
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_vline(xintercept = seq(-5, 10, 5), col = "grey95") +
    geom_hline(yintercept = seq(-10, 10, 5), col = "grey95") +
    geom_point(aes(col = expr), size = 1.5) +
  
    scale_color_gradient(low = "#f5f8ff", high = color_name)
  
    # scale_color_gradient2(low = "white", mid = "#e5ecfb", high = color_name,
    #                       midpoint = 0)
    # 
  return(p)
}

# compute MDS of go term semantics
get_GOterm_semantic_mds <- function(goterm_enrich_file) {
  goerms_bp_elim <- read.table(goterm_enrich_file, sep = "\t", header = TRUE, quote = "") %>%
    filter(classicFisher < 0.05)
  
  annoDb_sycon <- read.table("09_gene_annotation/GOterms_OMA.tsv", sep = "\t", skip = 4) %>%
    select(V2, V5, V9) %>%
    mutate(V9 = case_when(V9 == "C" ~ "CC",
                          V9 == "P" ~ "BP",
                          V9 == "F" ~ "MF",
                          TRUE ~ "")) %>%
    rename(c("V2" = "GENE", "V5" = "GO", "V9" = "ONTOLOGY"))
  annoDb_sycon
  
  semantic_data <- godata(annoDb = annoDb_sycon, ont = "BP", computeIC = FALSE)
  
  semantic_similarity <- mgoSim(goerms_bp_elim$GO.ID, goerms_bp_elim$GO.ID,
                                semData = semantic_data, measure = "Wang", combine = NULL)
  
  semantic_dissimilarity <- sim2diss(semantic_similarity, method = "confusion", to.dist = TRUE)
  
  mds <- mds(semantic_dissimilarity)
  
  tibble_to_plot <- mds$conf %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "GO") %>%
    left_join(goerms_bp_elim, by = join_by(GO == GO.ID)) %>%
    arrange(desc(classicFisher))
  
  return(tibble_to_plot)
}

# plot MDS GO term
plot_semantic_mds <- function(file, labels_list) {
  df <- get_GOterm_semantic_mds(file)
  
  p <- ggplot(df, aes(D1, D2)) +
    geom_vline(xintercept = seq(-0.5, 0.5, 0.5), col = "grey95") +
    geom_hline(yintercept = seq(-0.5, 0.5, 0.5), col = "grey95") +
    geom_point(aes(size = Significant, fill = -log(classicFisher)),
               shape = 21, color = "#7e7bc7", alpha = 0.7)
  
  # add GO labels
  for (lbl in labels_list) {
    p <- p + geom_label(
      data = df %>% filter(GO %in% lbl$GO) %>% summarise(D1 = mean(D1), D2 = mean(D2), .groups = "drop"),
      label = lbl$label,
      fill = alpha("white", 0.8), size = 3,
      label.padding = unit(0.3, "lines")
    )
  }
  
  # Add common scales, guides, and themes
  p <- p +
    labs(x = "MD1", y = "MD2", fill = "-log(p-val)", size = "Sig. genes") +
    scale_size_continuous(range = c(1,14), breaks = c(5,10)) +
    scale_fill_gradient(high = "darkblue", low = "lightblue") +
    guides(size = guide_legend(override.aes = list(color = "#706db8")),
           color = guide_colorbar(barheight = 4, barwidth = 1))
  
  return(p)
}


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
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))))

umap_coord_x <- list(annotation_custom(grob = grid::textGrob(label = "UMAP 1",
                                                             x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                                             just = c(0, 1), gp = grid::gpar(fontsize = 18))))

umap_coord_y <- list(annotation_custom(grob = grid::textGrob(label = "UMAP 2",
                                                             x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                                             just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 18))))

md_arrows <- list(
  
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
  
  annotation_custom(grob = grid::textGrob(label = "MD 1",
                                          x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  
  annotation_custom(grob = grid::textGrob(label = "MD 2",
                                          x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  
  coord_cartesian(clip = "off"))

md_coord_x <- list(annotation_custom(grob = grid::textGrob(label = "MD 1",
                                                           x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                                           just = c(0, 1), gp = grid::gpar(fontsize = 10))),
                   coord_cartesian(clip = "off"))

md_coord_y <- list(annotation_custom(grob = grid::textGrob(label = "MD 2",
                                                           x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                                           just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
                   coord_cartesian(clip = "off"))


##########################################
#     hdWGCNA with standard pipeline     #
##########################################

DefaultAssay(Sycon) <- "RNA"
Sycon

# standard hdWGCNA pipeline
# pipeline with variable features is better than the other ("fraction")
# it retrieves better defined clusters
Sycon <- SetupForWGCNA(Sycon,
                       gene_select = "variable",
                       wgcna_name = "RNA_variableFeatures")

# Sycon <- SetupForWGCNA(Sycon,
#                        gene_select = "fraction",
#                        fraction = 0.05,
#                        wgcna_name = "RNA_fraction")

Sycon@misc$active_wgcna <- "RNA_variableFeatures"

# compute metacells by splitting both by samples and seurat clusters
# reduce number of cells to consider a cluster and the number of overlapping cells to 10
# otherwise it get rids of cluster 28..32
Sycon <- MetacellsByGroups(seurat_obj = Sycon,
                           group.by = c("orig.ident", "seurat_clusters"),
                           ident.group = "seurat_clusters",
                           min_cells = 10, k = 10,
                           assay = "RNA",
                           wgcna_name = "RNA_variableFeatures")

Sycon <- NormalizeMetacells(Sycon)

# it is important to use the RNA assay for DGE analysis, in order to use raw counts
# instead of normalized and/or integrated matrices
Sycon <- SetDatExpr(Sycon,
                    group_name = levels(Sycon[[]]$seurat_clusters),
                    group.by = "seurat_clusters",
                    assay = "RNA",
                    layer = "data")

# test different soft powers
Sycon <- TestSoftPowers(Sycon)

# # plot the results
# plotSoftPower_list <- PlotSoftPowers(Sycon)
# 
# # assemble with patchwork
# patchwork::wrap_plots(plotSoftPower_list, ncol = 2)

# construct co-expression network
Sycon <- ConstructNetwork(Sycon,
                          tom_outdir = "12_hdWGCNA/01_RNA_assay/",
                          tom_name = "RNA_variableFeatures", overwrite_tom = TRUE)

# PlotDendrogram(Sycon)

# this is necessary to include gorup by samples  in the next step
Sycon <- ScaleData(Sycon,
                   features = VariableFeatures(Sycon))

# compute all MEs in the full single-cell dataset
Sycon <- ModuleEigengenes(Sycon,
                          group.by.vars = "orig.ident")

# compute eigengene-based connectivity (kME)
Sycon <- ModuleConnectivity(Sycon,
                            group.by = "seurat_clusters")


########################
#     PLOT MODULES     #
########################

# make a featureplot of hMEs for each module
# here, plotted ranges are symmetrical (as hdWGCNA default settings)
RNA_featurePlot_list <- ModuleFeaturePlot(Sycon,
                                          features = "hMEs",
                                          order = TRUE)

# plot patchworked plots
patchwork::wrap_plots(RNA_featurePlot_list,
                      ncol = 3)

# extract hME values and add them to meta data
hMEs <- Sycon %>%
  GetMEs(harmonized = TRUE)
Sycon@meta.data <- cbind(Sycon@meta.data, hMEs)

# create a custom UMAP for each module
# here, plotted ranges are not corrected, thus not symmetrical
RNA_umap_list <- list()
for (i in names(RNA_featurePlot_list)) {
  RNA_umap_list[[i]] <- plot_module_UMAP(i, i) + theme_bw()
}

# plot patchworked plots
umaps <- print(patchwork::wrap_plots(RNA_umap_list, ncol = 3))
umaps

ggsave("12_hdWGCNA/01_RNA_assay/gene_module_umaps.pdf",
       umaps,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

ggsave("12_hdWGCNA/01_RNA_assay/gene_module_umaps.png",
       umaps,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')


############################################
#     EXPORT LISTS OF GENE PER MODULES     #
############################################

# get the list of genes per module
gene_list <- Sycon@misc$RNA_variableFeatures$wgcna_modules %>%
  group_by(module) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  { setNames(.$genes, .$module) }

# write the list of genes per module to separate files
for (module in names(gene_list)) {
  file_path <- file.path("12_hdWGCNA/01_RNA_assay/", paste0(module, "_RNA_module_genes.ls"))
  writeLines(gene_list[[module]],
             file_path)
}

saveRDS(Sycon, "12_hdWGCNA/Sycon_hdWGCNA.Rds")

RNA_modules <- GetModules(Sycon)
RNA_mods <- levels(RNA_modules$module); mods <- mods[mods != 'grey']

RNA_dotplot <- Sycon %>%
  DotPlot(features = RNA_mods, group.by = 'seurat_clusters')


RNA_dotplot_custom <- RNA_dotplot$data %>%
  mutate(id = factor(id,
                     levels = c("32", "27", "24", "19", "17",
                                "30", "25", "15",
                                "26", "21",
                                "23", "18", "16", "14", "11",
                                "31", "29", "28", "22", "20", "13", "12",
                                seq(10, 0, -1))),
         features.plot = factor(features.plot,
                                levels = c("yellow", "turquoise", "brown", "green",
                                           "blue", "red"))) %>%
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

ggsave("12_hdWGCNA/01_RNA_assay/RNA_gene_module_dotplot.pdf",
       RNA_dotplot_custom,  device = cairo_pdf,
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')

ggsave("12_hdWGCNA/01_RNA_assay/RNA_gene_module_dotplot.png",
       RNA_dotplot_custom,  device = "png",
       dpi = 300, height = 6, width = 12, units = ("in"), bg = 'white')


#################################
#     UMAPS FOR PUBLICATION     #
#################################

SeuratExtend:::color_iwh("intense", n = 5)
#"#996742" "#5e824b" "#5d7880" "#7169a7" "#9f516c"
#"#c7956d" "#8fae6c" "#6db7ad" "#879dcf" "#cc88ad"
# "#b88340" "#9749ba" "#659d56" "#b64a56" "#757bad"

modules <- c("yellow","turquoise","blue","green","brown","red")
# cols <- c("#e1c23a","#0faef2","#0f4ff2","#1fde22","#aa2831","#e0080b") # these colors to keep the correspondence with hdWGCNA scheme
cols <- c("#48B1A7","#747CC9","#67A64E","#CA5F3E","#C8577B","#67A64E") # these colors to keep the correspondence with cluster annotations


module_umaps <- Map(plot_module_UMAP, modules, cols)
module_umaps <- lapply(module_umaps, function(p)
  p +
    scale_y_continuous(breaks = seq(-10, 10, 5)) +
    # coord_cartesian(clip = "off") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none",
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
          axis.title.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x = element_text(face = "bold", angle = 0, hjust = 0),
          strip.clip = "off")
)


module_umaps$red <- module_umaps$red +
  annotate("text", x = -7.3, y = 9.9, hjust = 0,
           label = "Module 1")

module_umaps$blue <- module_umaps$blue +
  annotate("text", x = -7.3, y = 9.9, hjust = 0,
           label = "Module 2")

module_umaps$green <- module_umaps$green +
  annotate("text", x = -7.3, y = 9.9, hjust = 0,
           label = "Module 3")

module_umaps$brown <- module_umaps$brown +
  annotate("text", x = -7.3, y = 9.9, hjust = 0,
           label = "Module 4") +
  scale_x_continuous(breaks = seq(0, 10, 5)) +
  scale_y_continuous(breaks = seq(-5, 5, 5)) +
  umap_coord_x + umap_coord_y

module_umaps$yellow <- module_umaps$yellow +
  annotate("text", x = -7.3, y = 9.9, hjust = 0,
           label = "Module 5")

module_umaps$turquoise <- module_umaps$turquoise +
  annotate("text", x = -7.3, y = 9.9, hjust = 0,
           label = "Module 6")
  
  
# module_umaps[[1]] <- module_umaps[[1]] +
#   annotate("text", x = -7.3, y = 11.5, hjust = 0,
#            label = "Oocytes/early embryos")
# 
# module_umaps[[2]] <- module_umaps[[2]] +
#   annotate("text", x = -7.3, y = 11.5, hjust = 0,
#            label = "Embryos")
# 
# module_umaps[[3]] <- module_umaps[[3]] +
#   annotate("text", x = -7.3, y = 11.5, hjust = 0,
#            label = "Choanocytes")
# 
# module_umaps[[4]] <- module_umaps[[4]] +
#   annotate("text", x = -7.3, y = 11.5, hjust = 0,
#            label = "Sclerocytes")
# 
# module_umaps[[5]] <- module_umaps[[5]] +
#   annotate("text", x = -7.3, y = 11.5, hjust = 0,
#            label = "Choanocytes") +
#   scale_x_continuous(breaks = seq(0, 10, 5)) +
#   scale_y_continuous(breaks = seq(-5, 5, 5)) +
#   umap_coord_x + umap_coord_y

panel_def <-
  module_umaps$red +
  module_umaps$blue +
  module_umaps$green +
  module_umaps$brown +
  module_umaps$yellow + 
  module_umaps$turquoise +
  patchwork::plot_layout(nrow = 2)
panel_def

ggsave("12_hdWGCNA/01_RNA_assay/fig2_hdWGCNA_new.png",
       panel_def, device = "png",
       height = 5.201*2, width = 7.600*2, units = "in", dpi = 300, bg = "white")
ggsave("12_hdWGCNA/01_RNA_assay/fig2_hdWGCNA_new.pdf",
       panel_def, device = cairo_pdf,
       height = 5.201*2, width = 7.600*2, units = "in", dpi = 300, bg = "white")


#######################################
#     ENRICHMENTS FOR PUBLICATION     #
#######################################

# list of GO term files
module_enrich_files <- c("12_hdWGCNA/01_RNA_assay/yellow_RNA_module_genes_GOterms_topGO_BP_elim.txt",
                         "12_hdWGCNA/01_RNA_assay/turquoise_RNA_module_genes_GOterms_topGO_BP_elim.txt",
                         "12_hdWGCNA/01_RNA_assay/green_RNA_module_genes_GOterms_topGO_BP_elim.txt",
                         "12_hdWGCNA/01_RNA_assay/brown_RNA_module_genes_GOterms_topGO_BP_elim.txt",
                         "12_hdWGCNA/01_RNA_assay/red_RNA_module_genes_GOterms_topGO_BP_elim.txt")

# define GO labels for each module
module_GO_labels <- list(
  
  # Module 1
  list(list(GO = c("GO:0030833", "GO:0007015"), label = "Actin filament\norganization"),
       list(GO = c("GO:0030334", "GO:2000145", "GO:0040012", "GO:0040011"), label = "Cell\nmigration"),
       list(GO = c("GO:0007131", "GO:0035825", "GO:0140527", "GO:1903046",
                   "GO:0140013", "GO:0061982", "GO:0007127", "GO:0051321"), label = "Meiotic\nprocesses")),
  
  # Module 2
  list(list(GO = c("GO:0046942", "GO:0015711", "GO:0050877"), label = "Organic molecule\nmetabolism"),
       list(GO = c("GO:0032836","GO:0021554","GO:0097065","GO:0033333","GO:0033338","GO:0060039",
                   "GO:0048589","GO:0048468","GO:0048736","GO:0010720","GO:0048638","GO:0051094",
                   "GO:0072163","GO:0072164","GO:0001823","GO:0001657","GO:0060537","GO:0007399",
                   "GO:0010975","GO:0002064","GO:0072006","GO:0048731","GO:0007517","GO:0043588",
                   "GO:0050793","GO:0048856"), label = "Anatomical\ndevelopment"),
       list(GO = c("GO:0007160", "GO:0031589", "GO:0042221", "GO:0007155",
                   "GO:0030198", "GO:0007156", "GO:0007166"), label = "Cell adhesion\nand proliferation")),
  
  # Module 3
  list(list(GO = c("GO:0007009", "GO:0071554", "GO:0006366"), label = "Cell membrane\norganization"),
       list(GO = c("GO:0007229", "GO:0007155"), label = "Cell adhesion"),
       list(GO = c("GO:0009913", "GO:0008544", "GO:0060429"), label = "Epithelium\ndevelopment")),
  
  # Module 4
  list(# list(GO = c("GO:0050801", "GO:0055080", "GO:0003169"), label = "Ion\nhomeostasis"),
       list(GO = c("GO:0034220", "GO:0008656", "GO:0002600", "GO:0006820", "GO:0005698"), label = "Ion transmembrane\ntransport"),
       list(GO = c("GO:0098609", "GO:0034446", "GO:0007160", "GO:0033627", "GO:0007156"), label = "Cell adhesion"),
       list(GO = c("GO:0009718", "GO:0071895"), label = "Mineralized tissue\ndevelopment"),
       list(GO = c("GO:0016055", "GO:0051493"), label = "Cytoskeleton\norganization")),
  
  # Module 5
  list(list(GO = c("GO:0006915", "GO:0007166", "GO:0043122", "GO:0071356"), label = "Cell adhesion\nand signalling"))
)

# generate all plots in one go
module_mds <- map2(module_enrich_files, module_GO_labels, plot_semantic_mds)
names(module_mds) <- names(module_umaps)
module_mds <- lapply(module_mds, function(p)
  p +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(legend.position = "none",
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(color = "#545454"),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "#545454"),
          plot.margin = margin(2, 0, 3, 8, "mm")
          )
)

module_mds[[2]] <- module_mds[[2]] +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5))

module_mds[[3]] <- module_mds[[3]] +
  scale_y_continuous(breaks = c(-0.5, 0))

module_mds[[4]] <- module_mds[[4]] +
  scale_y_continuous(breaks = c(-0.5, 0.5))

module_mds[[5]] <- module_mds[[5]] +
  scale_x_continuous(breaks = c(0, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.5)) +
  md_coord_x + md_coord_y


#################
#     PANEL     #
#################

panel_def <-
  module_umaps$yellow + module_mds$yellow +
  module_umaps$turquoise + module_mds$turquoise +
  module_umaps$green + module_mds$green +
  module_umaps$brown + module_mds$brown +
  module_umaps$red + module_mds$red +
  patchwork::plot_layout(ncol = 2, widths = c(0.8, 1)) +
  patchwork::plot_annotation(tag_levels = list(c("A", "B",
                                                 "D", "E",
                                                 "G", "H",
                                                 "J", "K",
                                                 "M", "N"))) &
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.location = "plot")

# panel_def

ggsave("fig2.png",
       panel_def, device = "png",
       width = 10/1.5, height = 22/1.5, units = "in", dpi = 300, bg = "white")