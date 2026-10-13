#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(ggpubr)
library(sceasy)

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")


##################
#     INPUTS     #
##################

# directory where cellranger results are stored
input_directory <- getwd()

# pre-processed porifera data
scil_data <- "00_input/Sycon_Seuratv4.Rdata"
# aque_data <- "04_preprocessed_scRNAseqs/Aque_cellFiltered.h5ad"
# slac_data <- "04_preprocessed_scRNAseqs/Slac_cellFiltered.h5ad"

# files with cell-cycle markers annotation
cellcycle_markers <- list.files(path = "00_input/cell_cycle_markers", pattern = "withNCBI.csv", full.names = TRUE)

# cellcycle markers per reference orthogroup
cellcycle_orthogroupsPerReference_file <- "16_cellcycle_scoring/03_orthofinder/cellcycle_orthogroupsPerReference.tsv"

# cellcycle markers per porifera orthogroup
cellcycle_orthogroupsPerSycon_file <- "16_cellcycle_scoring/03_orthofinder/cellcycle_orthogroupsPerSycon.tsv"
# cellcycle_orthogroupsPeraque_file <- "16_cellcycle_scoring/03_orthofinder/cellcycle_orthogroupsPerAmphimedon.tsv"
# cellcycle_orthogroupsPerslac_file <- "16_cellcycle_scoring/03_orthofinder/cellcycle_orthogroupsPerSpongilla.tsv"


###################
#     OUTPUTS     #
###################

PCA_cellcycle_scoring_file_scil <- "16_cellcycle_scoring/04_scoring_results/PCA_cellcycle_scoring_scil"
PCA_cellcycle_scoring_file_aque <- "16_cellcycle_scoring/04_scoring_results/PCA_cellcycle_scoring_aquemedon.pdf"
PCA_cellcycle_scoring_file_slac <- "16_cellcycle_scoring/04_scoring_results/PCA_cellcycle_scoring_slacilla.pdf"


#####################
#     FUNCTIONS     #
#####################

import_cellCycle_annotation_tables <- function(tsv_to_read, colnames) {
  
  tibble <- read.table(tsv_to_read, sep = "\t", na.strings = "NA") %>%
    `colnames<-`(c(colnames[1], colnames[2])) %>%
    tidyr::separate_rows(colnames[2], sep = ",") %>%
    mutate_if(is.character, ~na_if(., '')) %>%
    drop_na()
  
  return(tibble)
}

annotate_slace_cellCycle_markers <- function(reference_markers, slace_markers) {
  
  tibble <- reference_cellCycle_markers %>%
    select(-c(modified,geneID)) %>%
    left_join(reference_markers) %>%
    left_join(slace_markers, by = "OGs", relationship = "many-to-many") %>%
    drop_na() %>%
    group_by(geneID) %>%
    reframe(phase) %>%
    distinct() %>%
    arrange(phase)
  
  return(tibble)
}

write_annotation_to_file <- function(input_table, output_filename) {
  
  write.table(x = input_table, file = output_filename,
              col.names = TRUE, row.names = FALSE,
              sep = "\t", quote = FALSE)
  
}

extract_phase_markers <- function(list_of_markers, cell_phase) {
  
  genes <- list_of_markers %>%
    filter(phase == cell_phase) %>%
    pull(geneID)
  
  return(genes)
  
}

score_cellcycle <- function(seurat_object, s_genes, g2m_genes) {
  
  seurat_object <- seurat_object %>%
    CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE) %>%
    RunPCA(features = c(s_genes, g2m_genes))
  
  return(seurat_object)
}

plot_PCA <- function(seurat_object, plot_title) {
  
  plot <- DimPlot(seurat_object, reduction = "pca") +
    ggtitle(plot_title) +
    labs(colour="Cell cycle phases") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))
  
  return(plot)
}


###############################################################
#     IMPORT AND PROCESS DATA ABOUT CELL CYCLE ANNOTATION     #
###############################################################

# load orthogroup annotations for reference species
cellcycle_orthogroupsPerReference <- import_cellCycle_annotation_tables(cellcycle_orthogroupsPerReference_file,
                                                                        c("OGs", "NCBI"))

# load orthogroup annotations for slace species
cellcycle_orthogroupsPerSycon <- import_cellCycle_annotation_tables(cellcycle_orthogroupsPerSycon_file,
                                                                    c("OGs", "geneID"))
# cellcycle_orthogroupsPerslac <- import_cellCycle_annotation_tables(cellcycle_orthogroupsPerslac_file,
#                                                                     c("OGs", "geneID"))
# cellcycle_orthogroupsPeraque <- import_cellCycle_annotation_tables(cellcycle_orthogroupsPeraque_file,
#                                                                     c("OGs", "geneID"))

# create a tibble out of cell-cycle markers of reference species
reference_cellCycle_markers <- do.call(rbind, lapply(cellcycle_markers, read.csv))

# annotate scil cell-cycle marker genes
cellcycle_scil <- annotate_slace_cellCycle_markers(cellcycle_orthogroupsPerReference,
                                                     cellcycle_orthogroupsPerSycon)
# cellcycle_aque <- annotate_slace_cellCycle_markers(cellcycle_orthogroupsPerReference,
#                                                      cellcycle_orthogroupsPeraque)
# cellcycle_slac <- annotate_slace_cellCycle_markers(cellcycle_orthogroupsPerReference,
#                                                      cellcycle_orthogroupsPerslac)

# write annotation to files
write_annotation_to_file(cellcycle_scil, "16_cellcycle_scoring/03_orthofinder/cellcycle_markers_scil.tsv")
# write_annotation_to_file(cellcycle_aque, "16_cellcycle_scoring/03_orthofinder/cellcycle_markers_aquemedon.tsv")
# write_annotation_to_file(cellcycle_slac, "16_cellcycle_scoring/03_orthofinder/cellcycle_markers_slacilla.tsv")

# extract S phase markers
s_genes_scil <- extract_phase_markers(cellcycle_scil, "S")
# s_genes_aque <- extract_phase_markers(cellcycle_aque, "S")
# s_genes_slac <- extract_phase_markers(cellcycle_slac, "S")

# extract G2/M phase markers
g2m_genes_scil <- extract_phase_markers(cellcycle_scil, "G2/M")
# g2m_genes_aque <- extract_phase_markers(cellcycle_aque, "G2/M")
# g2m_genes_slac <- extract_phase_markers(cellcycle_slac, "G2/M")


###################################
#     IMPORT SINGLE CELL DATA     #
###################################

# load single cell data
load(scil_data)

# aque <- sceasy::convertFormat(aque_data, from = "anndata", to = "seurat") %>%
#   NormalizeData() %>%
#   FindVariableFeatures(selection.method = "vst") %>%
#   ScaleData(features = rownames(.))
# 
# slac <- sceasy::convertFormat(slac_data, from = "anndata", to = "seurat") %>%
#   NormalizeData() %>%
#   FindVariableFeatures(selection.method = "vst") %>%
#   ScaleData(features = rownames(.))

# score cells based on cycle phases
Sycon <- score_cellcycle(Sycon, s_genes_scil, g2m_genes_scil)
# aque <- score_cellcycle(aque, s_genes_aque, g2m_genes_aque)
# slac <- score_cellcycle(slac, s_genes_slac, g2m_genes_slac)

# plot PCs scored by cell-cycle phases before regression
PCAplot_before_cellCycle_scaling_scil <- plot_PCA(Sycon, "Cell cycle scoring before regression")
# PCAplot_before_cellCycle_scaling_aque <- plot_PCA(aque, "Cell cycle scoring before regression (aquemedon)")
# PCAplot_before_cellCycle_scaling_slac <- plot_PCA(slac, "Cell cycle scoring before regression (slacilla)")

# # plot umap before regression
# UMAPplot_before_cellCycle_scaling <- DimPlot(Sycon, reduction = "umap") +
#   # ggtitle("Cell cycle scoring before regression on UMAP") +
#   labs(colour="Cell cycle phases") +
#   theme_classic()
# 
# UMAPplot_before_cellCycle_scaling

# regress out cell cycle phases
Sycon <- Sycon %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
  RunPCA(features = c(s_genes_scil, g2m_genes_scil))

# aque <- aque %>%
#   ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
#   RunPCA(features = c(s_genes_aque, g2m_genes_aque))
# 
# slac <- slac %>%
#   ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(.)) %>%
#   RunPCA(features = c(s_genes_slac, g2m_genes_slac))

# run PCA again
Sycon <- Sycon %>%
  RunPCA(features = c(s_genes_scil, g2m_genes_scil))

# plot PCs scored by cell-cycle phases after regression
PCAplot_after_cellCycle_scaling_scil <- plot_PCA(Sycon, "Cell cycle scoring\nafter regression")
# PCAplot_after_cellCycle_scaling_aque <- plot_PCA(aque, "Cell cycle scoring after regression (aquemedon)")
# PCAplot_after_cellCycle_scaling_slac <- plot_PCA(slac, "Cell cycle scoring after regression (slacilla)")

# # plot umap after regression
# UMAPplot_after_cellCycle_scaling <- DimPlot2(Sycon, reduction = "umap") +
#   # ggtitle("Cell cycle scoring before regression on UMAP") +
#   labs(colour="Cell cycle phases") +
#   theme_classic()
# 
# UMAPplot_after_cellCycle_scaling

# plot panel
panel_scil <- ggpubr::ggarrange(PCAplot_before_cellCycle_scaling_scil + ggtitle("Cell cycle scoring\nbefore regression"),
                                PCAplot_after_cellCycle_scaling_scil,
                                # UMAPplot_before_cellCycle_scaling,
                                # UMAPplot_after_cellCycle_scaling,
                                common.legend = TRUE, legend = "right", labels = "AUTO")
panel_scil
# panel_aque <- ggpubr::ggarrange(PCAplot_before_cellCycle_scaling_aque, PCAplot_after_cellCycle_scaling_aque,
#                                 ncol = 2, align = "hv", common.legend = TRUE, legend = "right", labels = "AUTO")
# panel_slac <- ggpubr::ggarrange(PCAplot_before_cellCycle_scaling_slac, PCAplot_after_cellCycle_scaling_slac,
#                                 ncol = 2, align = "hv", common.legend = TRUE, legend = "right", labels = "AUTO")



########################
#     SAVE OUTPUTS     #
########################

ggsave(paste0(PCA_cellcycle_scoring_file_scil, ".pdf"),
       plot = panel_scil, device = "pdf",
       dpi = 300, height = 4, width = 9, units = ("in"), bg = 'white')
ggsave(paste0(PCA_cellcycle_scoring_file_scil, ".png"),
       plot = panel_scil, device = "png",
       dpi = 300, height = 4, width = 9, units = ("in"), bg = 'white')

ggsave(PCA_cellcycle_scoring_file_aque,
       plot = panel_aque, device = "pdf",
       dpi = 300, height = 6, width = 14, units = ("in"), bg = 'white')

ggsave(PCA_cellcycle_scoring_file_slac,
       plot = panel_slac, device = "pdf",
       dpi = 300, height = 6, width = 14, units = ("in"), bg = 'white')
