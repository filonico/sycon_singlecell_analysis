### SYCON PRE-PROCESSING

# Load in necessary libraries
library(Seurat)
library(tidyverse)
library(harmony)
library(dittoSeq)
library(RColorBrewer)

# Load in sample data
Sycon.S231 <- Read10X("raw_data/matrices/S231/")
Sycon.S232 <- Read10X("raw_data/matrices/S232/")
Sycon.S251 <- Read10X("raw_data/matrices/S251/")
Sycon.S252 <- Read10X("raw_data/matrices/S252/")

Sycon.S231 <- CreateSeuratObject(counts = Sycon.S231,
                                 project = "Replicate 1",
                                 min.cells = 0, min.features = 0)
# Convert the count matrix to a data frame for manipulation
Sycon.S231_matrix <- Sycon.S231[["RNA"]]$counts
Sycon.S231_df <- as.data.frame(Sycon.S231_matrix)

# Update row names with best annotations
row.names(Sycon.S231_df) <- Sycon.S231_df$best_annot

# Remove the 'best_annot' column before recreating the Seurat object if it's not needed as a feature
Sycon.S231_df <- Sycon.S231_df[, -which(colnames(Sycon.S231_df) == "best_annot")]

# Recreate the Seurat object
Sycon.S231 <- CreateSeuratObject(counts = as.matrix(Sycon.S231_df),
                                 project = "Replicate 1",
                                 min.cells = 0, min.features = 0)
rownames(Sycon.S231)

Sycon.S232 <- CreateSeuratObject(counts = Sycon.S232,
                                 project = "Replicate 2",
                                 min.cells = 0, min.features = 0)
# Convert the count matrix to a data frame for manipulation
Sycon.S232_matrix <- Sycon.S232[["RNA"]]$counts
Sycon.S232_df <- as.data.frame(Sycon.S232_matrix)

# Update row names with best annotations
row.names(Sycon.S232_df) <- Sycon.S232_df$best_annot

# Remove the 'best_annot' column before recreating the Seurat object if it's not needed as a feature
Sycon.S232_df <- Sycon.S232_df[, -which(colnames(Sycon.S232_df) == "best_annot")]

# Recreate the Seurat object
Sycon.S232 <- CreateSeuratObject(counts = as.matrix(Sycon.S232_df),
                                 project = "Replicate 2",
                                 min.cells = 0, min.features = 0)
rownames(Sycon.S232)

Sycon.S251 <- CreateSeuratObject(counts = Sycon.S251,
                                 project = "Replicate 3",
                                 min.cells = 0, min.features = 0)
# Convert the count matrix to a data frame for manipulation
Sycon.S251_matrix <- Sycon.S251[["RNA"]]$counts
Sycon.S251_df <- as.data.frame(Sycon.S251_matrix)

# Update row names with best annotations
row.names(Sycon.S251_df) <- Sycon.S251_df$best_annot

# Remove the 'best_annot' column before recreating the Seurat object if it's not needed as a feature
Sycon.S251_df <- Sycon.S251_df[, -which(colnames(Sycon.S251_df) == "best_annot")]

# Recreate the Seurat object
Sycon.S251 <- CreateSeuratObject(counts = as.matrix(Sycon.S251_df),
                                 project = "Replicate 3",
                                 min.cells = 0, min.features = 0)
rownames(Sycon.S251)

Sycon.S252 <- CreateSeuratObject(counts = Sycon.S252,
                                 project = "Replicate 4",
                                 min.cells = 0, min.features = 0)
# Convert the count matrix to a data frame for manipulation
Sycon.S252_matrix <- Sycon.S252[["RNA"]]$counts
Sycon.S252_df <- as.data.frame(Sycon.S252_matrix)

# Update row names with best annotations
row.names(Sycon.S252_df) <- Sycon.S252_df$best_annot

# Remove the 'best_annot' column before recreating the Seurat object if it's not needed as a feature
Sycon.S252_df <- Sycon.S252_df[, -which(colnames(Sycon.S252_df) == "best_annot")]

# Recreate the Seurat object
Sycon.S252 <- CreateSeuratObject(counts = as.matrix(Sycon.S252_df),
                                 project = "Replicate 4",
                                 min.cells = 0, min.features = 0)
rownames(Sycon.S252)

# Normalise the data
Sycon.S231 <- SCTransform(Sycon.S231)
Sycon.S232 <- SCTransform(Sycon.S232)
Sycon.S251 <- SCTransform(Sycon.S251)
Sycon.S252 <- SCTransform(Sycon.S252)

# Merge replicates together
sycon_merge <- merge(Sycon.S231, c(Sycon.S232, Sycon.S251, Sycon.S252))

# Filter cells
p <- multi_dittoPlot(sycon_merge,
                     var = c("nFeature_RNA", "nCount_RNA"),
                     group.by = "orig.ident",
                     plots = c("jitter", "vlnplot"),
                     legend.show = TRUE,
                     ncol = 1,
                     color.panel = c("orangered1", "purple1", "olivedrab3", "deepskyblue"))
ggsave(filename = "results/plots/1_Basic_metrics.pdf",
       plot = p, dpi = 600, height = 6, width = 5)

p <- dittoScatterPlot(sycon_merge,
                      x.var = "nCount_RNA",
                      y.var = "nFeature_RNA",
                      split.by = "orig.ident",
                      color.panel = c("orangered1", "purple1", "olivedrab3", "deepskyblue"),
                      color.var = "orig.ident")
ggsave(filename = "results/plots/2_Basic_metrics_scatterplot.pdf",
       plot = p, dpi = 600, height = 4, width = 5)

# Filtering
sycon_merge <- subset(sycon_merge,
                      subset = nFeature_RNA > 200 & nCount_RNA < 10000 & nCount_RNA > 1500)

# Normalise the data
sycon_merge <- SCTransform(sycon_merge, return.only.var.genes = FALSE)
sycon_merge <- ScaleData(sycon_merge)

# Dimensionality reduction
## Merged datset
sycon_merge <- RunPCA(sycon_merge, npcs = 100)
ElbowPlot(sycon_merge, ndims = 100)
sycon_merge <- FindNeighbors(sycon_merge, dims = 1:30)
sycon_merge <- FindClusters(sycon_merge, resolution = 1.2)
sycon_merge <- RunUMAP(sycon_merge, dims = 1:30)
p <- dittoDimPlot(sycon_merge,
                  var = "orig.ident",
                  reduction = "umap",
                  do.label = T,
                  color.panel = c("orangered1", "purple1", "olivedrab3", "deepskyblue"))
print(p)
ggsave(filename = "results/plots/3_UMAP_merged.pdf",
       plot = p, dpi = 600, height = 4, width = 5)


saveRDS(sycon_merge, file = "results/sycon_merged.rds")

# Harmony Integration
sycon_integrated <- RunHarmony(sycon_merge, group.by.vars = c("orig.ident"))

sycon_integrated <- FindNeighbors(sycon_integrated, reduction = "harmony", dims = 1:30)
sycon_integrated <- FindClusters(sycon_integrated, resolution = 1.2)
sycon_integrated <- RunUMAP(sycon_integrated, dims = 1:30, reduction = "harmony")
sycon_integrated <- RunTSNE(sycon_integrated, dims = 1:30, reduction = "harmony")

p <- dittoDimPlot(sycon_integrated,
                  var = "orig.ident",
                  reduction = "umap",
                  do.label = T,
                  color.panel = c("orangered1", "purple1", "olivedrab3", "deepskyblue"))
print(p)
ggsave(filename = "results/plots/4_UMAP_harmonyintegrated.pdf", plot = p, dpi = 600, height = 4, width = 5)

p <- dittoDimPlot(sycon_integrated,
                  var = "orig.ident",
                  reduction = "tsne",
                  do.label = T,
                  color.panel = c("orangered1", "purple1", "olivedrab3", "deepskyblue"))
ggsave(filename = "results/plots/4_TSNE_harmonyintegrated.pdf",
       plot = p, dpi = 600, height = 4, width = 5)

# Function to generate and randomize an extended color palette
extend_and_randomize_palette <- function(colors, num) {
  reps <- ceiling(num / length(colors))
  pal <- rep(colors, reps)
  randomized_pal <- sample(pal[1:num])  # Randomize the colors
  return(randomized_pal)
}

# Combining multiple palettes
palette1 <- brewer.pal(12, "Set3")
palette2 <- brewer.pal(12, "Paired")
palette3 <- brewer.pal(8, "Dark2")  # This combination adds up to 32 colors exactly

# Combining and randomizing the palettes
full_palette <- c(palette1, palette2, palette3)
extended_random_palette <- extend_and_randomize_palette(full_palette, 32)

p <- dittoDimPlot(sycon_integrated,
                  var = "seurat_clusters",
                  do.label = T,
                  reduction = "umap",
                  color.panel = as.vector(extended_random_palette))
print(p)
ggsave(filename = "results/plots/6_UMAP_harmonyintegrated_clusters.pdf",
       plot = p, dpi = 600, height = 4, width = 5)

p <- dittoDimPlot(sycon_integrated,
                  var = "seurat_clusters",
                  do.label = T,
                  reduction = "tsne",
                  color.panel = as.vector(extended_random_palette))
ggsave(filename = "results/plots/7_TSNE_harmonyintegrated_clusters.pdf",
       plot = p, dpi = 600, height = 4, width = 5)

saveRDS(sycon_integrated, file = "results/sycon_integrated.rds")