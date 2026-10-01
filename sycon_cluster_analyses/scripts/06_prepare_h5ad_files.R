#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(Seurat)
library(reticulate)
library(scCustomize)
library(tidyverse)
library(metacell)
# library(sceasy)

# initialize python from the conda env
reticulate::use_condaenv("RSeurat_env")

# import the anndata module to converto to/from .h5ad
ad <- reticulate::import("anndata", convert = FALSE)


#####################
#     I/O FILES     #
#####################

output_path <- "./04_preprocessed_scRNAseqs/"
input_path <- "./00_input/scRNAseqs_rawTables/"

# Sycon ciliatum files
scil_files <- list(Rdata = "00_input/Sycon_Seuratv4.Rdata",
                   out_h5ad = "Scil_cellFiltered.h5ad")

# Amphimedon queenslandica files
aque_files <- list(UMItable = paste0(input_path, "A_queenslandica/GSM3021561_Aque_Amphimedon_queenslandica_adult_UMI_table.txt.gz"),
                   metacell_assignments = paste0(input_path, "A_queenslandica/Amphimedon_adult_metacell_assignments.tsv"),
                   metacell_annotation = paste0(input_path, "A_queenslandica/Aque_adult_metacell_annotation.tsv"),
                   out_h5ad = "Aque_cellFiltered.h5ad")

# Spongilla lacustris files
# SPONGILLA LACUSTRIS FILE ARE NOW IN ../spongilla_remapping/
# WE REMAPPED THE SINGLE-CELL DATASET AGAINST THE GENOME
#slac_files <- list(UMItable = paste0(input_path, "S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix_edited.txt.gz"),
#                   cell_metadata = paste0(input_path, "S_lacustris/spongilla_cell_metadata_newNames.tsv"),
#                   out_h5ad = "Slac_cellFiltered.h5ad")


#####################
#     FUNCTIONS     #
#####################

plot_featuresVScounts <- function(seuratObject,
                                  nCount_min, nCount_max,
                                  nFeature_min, nFeature_max) {
  
  # if none of the filtering parameter is set, then plot non-filtered data
  if (missing(nCount_min) && missing(nCount_max) && missing(nFeature_min) && missing(nFeature_max)) {
    seuratObject_toPlot <- seuratObject
  }
  
  # else, filter low-quality cell and then plot to data
  else {
    seuratObject_toPlot <- seuratObject %>%
      subset(subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max &
               nCount_RNA > nCount_min & nCount_RNA < nCount_max)
  }
  
  # produce violin plot
  vlnplot <- VlnPlot(seuratObject_toPlot,
                     features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident")
  
  # produce scatter plot
  scatterplot <- FeatureScatter(seuratObject_toPlot,
                                feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #+
    # scale_x_log10() +
    # scale_y_log10()

  # put the plots in a panel
  panel <- ggpubr::ggarrange(vlnplot, scatterplot)
  
  if (missing(nCount_min) && missing(nCount_max) && missing(nFeature_min) && missing(nFeature_max)) {
    return(panel)
  }
  else {
    return(list(S.object = seuratObject_toPlot, panel = panel))
  }
}

# function to convert from Seurat to AnnData object and save to file
fromSeurat_toAnndata <- function(s.object, outdir, outfile) {
  
  scCustomize::as.anndata(x = s.object, main_layer = "counts",
                          other_layers = NULL, file_path = outdir,
                          file_name = outfile)
}


#######################
#     S. CILIATUM     #
#######################

load(scil_files$Rdata)

scil_S.object <- Seurat::CreateSeuratObject(counts = Seurat::GetAssayData(object = Sycon,
                                                                          assay = "RNA",
                                                                          layer = 'counts'),
                                  project = "scil", min.cells = 0, min.features = 0) %>%
  Seurat::AddMetaData(metadata = Sycon@meta.data)

# sceasy::convertFormat(scil_files_reAnno$in_h5ad, from = "anndata", to = "seurat", outFile = "Sycon.ReAnno.rds")
# scil_S.object_reAnno <- readRDS("Sycon.ReAnno.rds")

# convert SeuratObject to anndata file
fromSeurat_toAnndata(scil_S.object,
                       output_path, scil_files$out_h5ad)


############################
#     A. QUEENSLANDICA     #
############################

# create the Seurat object from UMI table
aque_S.object <- read.table(file = gzfile(aque_files$UMItable),
                            header = TRUE, sep = "\t") %>%
  Seurat::CreateSeuratObject(project = "aque_adult", min.cells = 0, min.features = 0,
                             meta.data = read.table(aque_files$metacell_assignments, header = TRUE, sep = "\t") %>%
                               left_join(read.table(aque_files$metacell_annotation, header = TRUE, sep = "\t")) %>%
                               rename(cell_type = ID) %>%
                               column_to_rownames(var = "Cell"))

# plot features and counts before cell filtering
aque_beforeFiltering <- plot_featuresVScounts(aque_S.object)

# plot features and counts after filtering
aque_afterFiltering <- plot_featuresVScounts(aque_S.object,
                                             nCount_min = 200, nCount_max = 11000,
                                             nFeature_min = 200, nFeature_max = 4000)

# convert SeuratObject to anndata file
fromSeurat_toAnndata(aque_afterFiltering$S.object,
                     output_path, aque_files$out_h5ad)

# SCTransform data
aque_afterFiltering$S.object <- aque_afterFiltering$S.object %>%
  SCTransform(verbose = TRUE)

# dim reductions
aque_afterFiltering$S.object <- aque_afterFiltering$S.object %>%
  RunPCA(verbose = TRUE)
aque_afterFiltering$S.object <- aque_afterFiltering$S.object %>%
  RunUMAP(dims = 1:30, verbose = TRUE)

fromSeurat_toAnndata(aque_afterFiltering$S.object,
                     output_path, "Aque_cellFiltered_SCT_PCA_UMAP.h5ad")


########################
#     S. LACUSTRIS     #
########################

#  OLD CODE TO PROCESS THE ALREADY-MAPPED SC-RNASEQ DATASET

# create the Seurat object
# slac_S.object <- read.table(file = gzfile(slac_files$UMItable),
#                             header = TRUE, sep = "\t") %>%
#   Seurat::CreateSeuratObject(project = "slac", min.cells = 0, min.features = 0,
#                             
#                              # add meta data with pre-computed cell clusters and annotations
#                              meta.data = read.table(slac_files$cell_metadata, header = TRUE, sep = "\t") %>%
#                                # select(cell, clusterID, cell_type, cell_type_newName) %>%
#                                mutate(cell = str_replace_all(cell, c('-' = '.'))) %>%
#                                column_to_rownames(var = "cell"))

# # plot features and counts before cell filtering
# slac_beforeFiltering <- plot_featuresVScounts(slac_S.object)

# # plot features and counts after cell filtering
# slac_afterFiltering <- plot_featuresVScounts(slac_S.object,
#                                              nCount_min = 200, nCount_max = 40000,
#                                              nFeature_min = 200, nFeature_max = 6000)

# # convert SeuratObject to anndata file
# fromSeurat_toAnndata(slac_afterFiltering$S.object,
#                      output_path, slac_files$out_h5ad)

# # SCTransform data
# slac_afterFiltering$S.object <- slac_afterFiltering$S.object %>%
#   SCTransform(verbose = TRUE)

# # dim reductions
# slac_afterFiltering$S.object <- slac_afterFiltering$S.object %>%
#   RunPCA(verbose = TRUE)
# slac_afterFiltering$S.object <- slac_afterFiltering$S.object %>%
#   RunUMAP(dims = 1:30, verbose = TRUE)

# fromSeurat_toAnndata(slac_afterFiltering$S.object,
#                      output_path, "Slac_cellFiltered_SCT_PCA_UMAP.h5ad")