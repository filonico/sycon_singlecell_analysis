#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses")

library(Seurat)
library(tidyverse)
library(DESeq2)


# read in the list of connecting gene pairs
gene_pairs_table <- read.table("05NEW_SAMap_porifera/02_gene_pairs/AqueScilSlac_leiden3Clusters_costum_04threhsold_100topCells_allCluster_samapGenePairs.tsv",
                               fill = TRUE, na.strings = "", header = TRUE, sep = "\t") %>%
  
  rename_with(~ paste0(., "_gene_pair"), .cols = !matches("_pval[12]$")) %>%
  
  pivot_longer(cols = everything(),
               names_to = c("mapping_clusters", ".value"),
               names_pattern = "^(.*)_(gene_pair|pval1|pval2)$") %>%
  drop_na() %>%
  mutate(mapping_clusters = str_replace(mapping_clusters, "\\.Scil", "|Scil"),
         mapping_clusters = str_replace(mapping_clusters, "\\.Slac", "|Slac"),
         mapping_clusters = str_replace(mapping_clusters, "\\.Aque", "|Aque")) %>%
  separate(gene_pair,
           into = c("source_gene", "target_gene"),
           sep = ";") %>% 
  separate(mapping_clusters,
           into = c("source_cluster", "target_cluster"),
           sep = "\\|") %>% 
  
  # keep only pairs that include Scil (discards Aque vs Slac)
  filter(str_starts(source_cluster, "Scil") |
           str_starts(target_cluster, "Scil")) %>%
  
  # flag rows where Scil is currently in column 2 and swap them
  mutate(swap = str_starts(target_cluster, "Scil"),
         cluster_a = if_else(swap, target_cluster, source_cluster),
         cluster_b = if_else(swap, source_cluster, target_cluster),
         gene_a = if_else(swap, target_gene, source_gene),
         gene_b = if_else(swap, source_gene, target_gene),
         pval_a = if_else(swap, pval2, pval1),
         pval_b = if_else(swap, pval1, pval2)) %>%
  
  select(source_cluster = cluster_a, target_cluster = cluster_b,
         source_gene = gene_a, target_gene = gene_b,
         pval1 = pval_a, pval2 = pval_b)
  
gene_pairs_table
  
# read in mapping scores
mapping_scores <- read.table("05NEW_SAMap_porifera/01_mapping_scores/AqueScilSlac_leiden3Clusters_costum_100topCells_samapMappingTable.tsv",
                             header = TRUE, sep = "\t") %>%
  
  # create a long-format dataframe
  pivot_longer(-X, names_to = "target_ID", values_to = "mapp_score") %>%
  rename("source_ID" = "X") %>% 
  
  # filter out rows with mapping scores below the threshold
  filter(mapp_score > 0.4,
         str_detect(source_ID, "^Scil"))

# create a vector of mapping cluster pairs
significant_mapping_clusters <- mapping_scores %>% 
  mutate(mapping_cluster_pairs = paste0(source_ID, "|", target_ID)) %>%
  pull(mapping_cluster_pairs)
significant_mapping_clusters
  
# extract genes connecting the significant mapping clusters (> 0.4)
significant_gene_pairs <- gene_pairs_table %>% 
  mutate(mapping_cluster_pairs = paste0(source_cluster, "|", target_cluster)) %>%
  select(mapping_cluster_pairs, source_gene, target_gene) %>%
  filter(mapping_cluster_pairs %in% significant_mapping_clusters) %>%
  distinct() %>%
  mutate(source_gene = str_remove(source_gene, "^Scil_"),
         target_gene = str_remove(target_gene, "^Aque_"))

# create lists of connecting gene pairs according to their contribution
samap_genes_choanocytes <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(1, 7, 11, 16, 23))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_myopeptido <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(21))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_sclero <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(25))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_myopeptido <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(21))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_oocytes_embryos <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(17, 19, 24, 27, 32))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_embryos <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(15))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_uncertain <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(0, 5, 10, 26))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

samap_genes_unknown <- significant_gene_pairs %>%
  separate(mapping_cluster_pairs, into = c("Scil", "Other"),
           sep = "\\|") %>%
  filter(Scil %in% paste0("Scil_", c(2, 3, 4, 6, 8, 9, 12, 13, 14, 18, 20, 22, 28, 29, 31))) %>%
  distinct(source_gene) %>%
  pull(source_gene)

# ggvenn::ggvenn(
#   list("choanocytes" = samap_genes_choanocytes,
#        # "myopeptidocytes" = samap_genes_myopeptido,
#        # "sclerocytes" = samap_genes_sclero,
#        "oocytes_embryos" = samap_genes_oocytes_embryos,
#        "embryos" = samap_genes_embryos,
#        # "uncertain" = samap_genes_uncertain,
#        "unknown" = samap_genes_unknown),
#   text_size = 3, stroke_size = 0.8)

# create a list of lists
samap_gene_lists <- list("samap_choanocytes" = samap_genes_choanocytes,
                         "samap_myopeptidocytes" = samap_genes_myopeptido,
                         "samap_sclerocytes" = samap_genes_sclero,
                         "samap_oocytes_embryos" = samap_genes_oocytes_embryos,
                         "samap_embryos" = samap_genes_embryos,
                         "samap_uncertain" = samap_genes_uncertain,
                         "samap_unknown" = samap_genes_unknown)

# write to files
for (cluster_name in names(samap_gene_lists)) {
  file_path <- file.path("05NEW_SAMap_porifera/02_gene_pairs/01_lists_for_enrichments",
                         paste0(cluster_name, "_connectingGenes.ls"))
  writeLines(samap_gene_lists[[cluster_name]],
             file_path#,
             # row.names = FALSE,
             # col.names = FALSE,
             # quote = FALSE
  )
}




samap <- schard::h5ad2seurat("05NEW_SAMap_porifera/ScilSlac_leiden3Clusters_samap.h5ad")

identity <- read.table("07_notableGenes_clusterAnnotation/cluster_identity.tsv",
                       sep = "\t", header = TRUE) 

spongilla_identity <- c("0" = "Archaeocytes_1",
                      "1" = "Archaeocytes_2",
                      "2" = "Archaeocytes_3",
                      "3" = "Archaeocytes_4",
                      "4" = "Archaeocytes_5",
                      "5" = "Archaeocytes_4",
                      "6" = "Cluster_6",
                      "7" = "Cluster_7",
                      "8" = "Cluster_8",
                      "9" = "Myopeptidocytes",
                      "10" = "Cluster_10",
                      "11" = "Cluster_11",
                      "12" = "Cluster_12",
                      "13" = "Cluster_13",
                      "14" = "Metabolocytes",
                      "15" = "Myopeptidocytes",
                      "16" = "Cluster_16",
                      "17" = "Pinacocytes",
                      "18" = "Cluster_18",
                      "19" = "Cluster_19",
                      "20" = "Cluster_20",
                      "21" = "Archaeocytes_6",
                      "22" = "Cluster_22",
                      "23" = "Archaeocytes-like",
                      "24" = "Myopeptidocytes",
                      "25" = "Choanocytes/-blasts",
                      "26" = "Cluster_26",
                      "27" = "Cluster_27",
                      "28" = "Cluster_28",
                      "29" = "Cluster_29",
                      "30" = "Sclerocytes",
                      "31" = "Amoebocytes/Neuroid",
                      "32" = "Cluster_32",
                      "33" = "Cluster_33",
                      "34" = "Cluster_34",
                      "35" = "Basopinacocytes",
                      "36" = "Mesocytes",
                      "37" = "Granulocytes-like",
                      "38" = "Cluster_38",
                      "39" = "Cluster_39",
                      "40" = "Cluster_40",
                      "41" = "Mesocytes",
                      "42" = "Cluster_42") %>%
  as_tibble() %>%
  rownames_to_column() %>%
  mutate(rowname = as.integer(rowname),
         rowname = rowname - 1)
spongilla_identity
  
samap[[]] %>% 
  filter(species == "Slac") %>%
  count(Slac_leiden_clusters, Slac_seurat_clusters) %>%
  mutate(Slac_seurat_clusters = as.integer(Slac_seurat_clusters)) %>%
  left_join(spongilla_identity, by = join_by("Slac_seurat_clusters" == "rowname")) %>%
  group_by(Slac_leiden_clusters) %>%
  slice_max(order_by = n,
            n = 5) %>% View
