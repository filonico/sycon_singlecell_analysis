#!/usr/bin/env Rscript

setwd("/data/evassvis/fn76/sycon/sycon_clusterAnnotation/sycon_cluster_analyses/")

library(tidyverse)
library(Seurat)
library(SeuratExtend)
library(VennDiagram)
library(patchwork)

###############################
#     DEFINE PLOT THEMES      #
###############################

theme_for_plots <- theme(
  plot.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey90", lineend = "round"),
  panel.border = element_rect(colour = "black", linewidth = .6),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  legend.key.width = unit(.4, "cm"),
  legend.key.height = unit(.4, "cm"),
  legend.position = "right",
  legend.title = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 10),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = "black", linewidth = .4),
  axis.ticks.length = unit(0.10, "cm"),
  axis.text.x = element_text(color = "black", margin = margin(t = 4, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(color = "black", margin = margin(t = 0, r = 4, b = 0, l = 0)),
  axis.title.y = element_text(angle = 90, size = 13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(angle = 0, size = 13, margin = margin(t = 10, r = 0, b = 0, l = 0)),
  strip.text = element_text(color = "black", face = "bold", hjust = 0),
  strip.placement = "outside",
  strip.background = element_blank(),
  strip.clip = "off"
)

umap_arrows <- list(
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(12, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(0, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  annotation_custom(grob = grid::segmentsGrob(x0 = unit(0, "mm"), x1 = unit(0, "mm"),
                                              y0 = unit(0, "mm"), y1 = unit(12, "mm"),
                                              arrow = arrow(length = unit(2.5, "mm"), ends = "last", type = "open"),
                                              gp = grid::gpar(col = "black", fill = "black", lwd = 1))),
  annotation_custom(grob = grid::textGrob(label = "UMAP 1", x = unit(0, "mm"), y = unit(0, "mm") - unit(2.5, "mm"),
                                          just = c(0, 1), gp = grid::gpar(fontsize = 10))),
  annotation_custom(grob = grid::textGrob(label = "UMAP 2", x = unit(0, "mm") - unit(2.5, "mm"), y = unit(0, "mm"),
                                          just = c(0, 0), rot = 90, gp = grid::gpar(fontsize = 10))),
  coord_cartesian(clip = "off")
)

theme_for_UMAPS <- theme(
  plot.background = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10, face = "bold"),
  plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.75, face = "bold"),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)


#####################
#     LOAD DATA     #
#####################

PCA_1to1 <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_1_Sycon_Spongilla.tsv",
                       header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

MCA_1to1 <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_1_Sycon_Spongilla_present_on_the_metazoan_stem.tsv",
                       header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

PCA_1toMany <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_many_Sycon_Spongilla.tsv",
                          header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

MCA_1toMany <- read.table("17_comparisons_1to1_orthologues/01_orthologue_tables/1_to_many_Sycon_Spongilla_present_on_the_metazoan_stem.tsv",
                          header = TRUE, sep = "\t", na.strings = c("", ".")) %>%
  select(-c(required_nodes_present, required_single_nodes_present, collected_specs_present)) %>%
  as_tibble()

# check overlaps
venn.diagram(list("PCA_1to1" = PCA_1to1$orthogroup_id,
                  "MCA_1to1" = MCA_1to1$orthogroup_id,
                  "PCA_1toMany" = PCA_1toMany$orthogroup_id,
                  "MCA_1toMany" = MCA_1toMany$orthogroup_id),
             filename = NULL, fill = c("red", "green", "blue", "yellow"))


# slac conversion table
slac_gene_names_conversion <- read.table("../spongilla_remapping/00_input/slac_genome/genes_to_transcript.tsv",
                                         sep = "\t", na.strings = "N/A") %>%
  left_join(read.table("../spongilla_remapping/00_input/slac_genome/transcripts_to_cdss.tsv",
                       sep = "\t", na.strings = "N/A"),
            by = join_by("V1" == "V2")) %>%
  as_tibble() %>%
  rename("transcript_id" = "V1", "gene_id" = "V2", "protein_id" = "V1.y") %>%
  arrange(gene_id)

# integrated samap
ScilSlac_samap <- schard::h5ad2seurat("05NEW_SAMap_porifera/ScilSlac_leiden3Clusters_samap.h5ad")
ScilSlac_samap[[]]$Scil_leiden_clusters <- as.factor(ScilSlac_samap[[]]$Scil_leiden_clusters)
ScilSlac_samap[[]]$Slac_leiden_clusters <- as.factor(ScilSlac_samap[[]]$Slac_leiden_clusters)
ScilSlac_samap[[]]$leiden_clusters      <- as.factor(ScilSlac_samap[[]]$leiden_clusters)


###############################################
#     COUNT HOW MANY GENES PER ORTHOGROUP     #
###############################################

# create a df with genes per species per og
orthogroup_full_set <- bind_rows(PCA_1to1, PCA_1toMany, MCA_1to1, MCA_1toMany) %>%
  # keep only one line per OG
  distinct(orthogroup_id, .keep_all = TRUE) %>%
  pivot_longer(c(LDemCA__Spongilla_lacustris, LCalcCA__Sycon_ciliatum),
               names_to = "species", values_to = "gene_ID") %>%
  # fix species anem
  mutate(species = case_when(str_detect(species, "Spongilla") ~ "Slac",
                             str_detect(species, "Sycon") ~ "Scil",
                             TRUE ~ species)) %>%
  # split multi-copy OG into multiple lines
  separate_rows(gene_ID, sep = ";") %>%
  # adjust slac gene names
  mutate(gene_ID = str_remove(gene_ID, "\\.[0-9]+$")) %>%
  left_join(slac_gene_names_conversion, by = join_by("gene_ID" == "protein_id")) %>%
  # add gene names for slac
  mutate(gene_id = if_else(is.na(gene_id), gene_ID, gene_id)) %>%
  select(-c(transcript_id, gene_ID)) %>%
  # remove duplicated genes (coming from slac isoforms)
  distinct()

# count genes per OG
orthogroup_full_set_with_geneCounts <- orthogroup_full_set %>%
  group_by(orthogroup_id, species) %>%
  mutate(n_genes = n()) %>%
  ungroup()


#######################################
#     FIND MARKERS PER SPECIES        #
#######################################

# skip this block if rerunning and marker genes can be input from tables

# subset sycon and SCTransform
scil_samap <- ScilSlac_samap %>% subset(species == "Scil") %>%
  SCTransform()

# find markers on sycon specific leiden clusters
scil_samap_allMarkers <- scil_samap %>%
  FindAllMarkers(group.by = "Scil_leiden_clusters") %>%
  mutate(gene = str_replace(gene, "Scil-", ""))

write.table(scil_samap_allMarkers,
            "17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# subset spongilla and SCTransform
slac_samap <- ScilSlac_samap %>% subset(species == "Slac") %>%
  SCTransform()

# find markers on spongilla specific leiden clusters
slac_samap_allMarkers <- slac_samap %>%
  FindAllMarkers(group.by = "Slac_leiden_clusters") %>%
  mutate(gene = str_replace(gene, "Slac-", "")) %>%
  # add gene_id conversion
  left_join(slac_gene_names_conversion, by = join_by("gene" == "gene_id"),
            relationship = "many-to-many") %>%
  select(-c(transcript_id, protein_id))

write.table(slac_samap_allMarkers,
            "17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)


#####################################################################
#     COMPUTE NUMBER OF CLUSTERS WHERE EACH OG IS EXPRESSED         #
#####################################################################

scil_samap_allMarkers <- read.table("17_comparisons_1to1_orthologues/scil_samap_allMarkers_leidenClusters.tsv",
                                    header = TRUE, sep = "\t") %>% as_tibble()
slac_samap_allMarkers <- read.table("17_comparisons_1to1_orthologues/slac_samap_allMarkers_leidenClusters.tsv",
                                    header = TRUE, sep = "\t") %>% as_tibble()

# create a single dataframe with all markers per cluster
allMarkers_per_cluster <- bind_rows(scil_samap_allMarkers %>%
                                      mutate(cluster = paste0("Scil_", cluster)),
                                    slac_samap_allMarkers %>%
                                      mutate(cluster = paste0("Slac_", cluster)) %>%
                                      distinct()) %>%
  filter(p_val_adj < 0.05, avg_log2FC > 1) %>%
  select(cluster, gene)
allMarkers_per_cluster

# add cluster-expression info at the gene level
orthogroup_full_set_with_clusterCount <- orthogroup_full_set %>%
  left_join(allMarkers_per_cluster, by = join_by("gene_id" == "gene"),
            relationship = "many-to-many") %>%
  group_by(orthogroup_id, species) %>%
  mutate(n_clusters = n_distinct(cluster, na.rm = TRUE)) %>%
  ungroup()

# compute og_type classification
og_type_table <- orthogroup_full_set_with_geneCounts %>%
  distinct(orthogroup_id, species, n_genes) %>%
  pivot_wider(names_from = species, values_from = n_genes) %>%
  mutate(og_type = case_when(Scil == 1 & Slac == 1 ~ "1 to 1",
                             Scil > 1  & Slac == 1 ~ "Single copy\nin Spongilla",
                             Scil == 1 & Slac > 1 ~ "Single copy\nin Sycon",
                             TRUE ~ "Many to many")) %>%
  select(orthogroup_id, og_type)

# get the total number of per-species leiden clusters
n_clusters_total <- tibble(species = c("Scil", "Slac"),
                           n_clusters_total = c(nlevels(droplevels(ScilSlac_samap[[]]$Scil_leiden_clusters)),
                                                nlevels(droplevels(ScilSlac_samap[[]]$Slac_leiden_clusters))))
n_clusters_total

# combine everything
orthogroup_full_set_with_geneCounts_clusterCounts <- orthogroup_full_set_with_geneCounts %>%
  left_join(orthogroup_full_set_with_clusterCount %>%
              select(orthogroup_id, species, cluster, n_clusters) %>%
              distinct(),
            by = c("orthogroup_id", "species"),
            relationship = "many-to-many") %>%
  left_join(og_type_table, by = "orthogroup_id") %>%
  left_join(n_clusters_total, by = "species") %>%
  mutate(n_cluster_prop = n_clusters / n_clusters_total) %>%
  group_by(orthogroup_id) %>%
  # remove all OGs where at least one gene is expressed in any cluster
  filter(!any(n_clusters == 0)) %>%
  ungroup() %>%
  mutate(species_label = case_when(species == "Scil" ~ "S. ciliatum",
                                   species == "Slac" ~ "S. lacustris"),
         og_type = factor(og_type, levels = c("1 to 1", "Single copy\nin Sycon",
                                              "Single copy\nin Spongilla", "Many to many")))

orthogroup_full_set_with_geneCounts_clusterCounts


#############################################################
#     GENE-LEVEL BREADTH METRICS (participation/redundancy) #
#############################################################

gene_breadth <- orthogroup_full_set_with_clusterCount %>%
  group_by(orthogroup_id, species, gene_id) %>%
  summarise(genewise_n_clusters = n_distinct(cluster, na.rm = TRUE),
            .groups = "drop")
gene_breadth

og_breadth_metrics <- gene_breadth %>%
  group_by(orthogroup_id, species) %>%
  summarise(n_genes = n(),
            n_genes_asMarkers = sum(genewise_n_clusters > 0),
            n_marker_events = sum(genewise_n_clusters),
            .groups = "drop") %>%
  left_join(orthogroup_full_set_with_geneCounts_clusterCounts %>%
              distinct(orthogroup_id, species, n_clusters, n_cluster_prop, og_type),
            by = c("orthogroup_id", "species")) %>%
  mutate(
    # fraction of an OG's genes that are a marker of *something*
    # (doesn't distinguish overlap vs. distinct clusters -- see redundancy)
    participation = n_genes_asMarkers / n_genes,
    # total marker-events / distinct clusters marked; ==1 no overlap between
    # paralogs, >1 some clusters marked by multiple paralogs, 0 no marker anywhere
    redundancy    = n_marker_events / pmax(n_clusters, 1)
  )

og_breadth_metrics


##############################
#     EXPRESSION BREADTH     #
##############################

# add og_type broad category
df_to_plot <- og_breadth_metrics %>%
  select(orthogroup_id, species, og_type, n_cluster_prop) %>%
  mutate(og_type_broad = case_when(og_type == "1 to 1" ~ "Single copy",
                                   TRUE ~ "Multiple copy"),
         species_label = if_else(species == "Scil", "S. ciliatum", "S. lacustris")) %>%
  filter(!(species == "Scil" & og_type == "Single copy\nin Sycon"),
         !(species == "Slac" & og_type == "Single copy\nin Spongilla")) %>%
  drop_na() 

df_to_plot %>%
  mutate(og_type_broad = factor(og_type_broad,
                                levels = c("Single copy", "Multiple copy"))) %>%
  ggplot(aes(og_type_broad, n_cluster_prop)) +
  geom_jitter(aes(col = species_label), width = 0.15, height = 0.01) +
  geom_boxplot(outliers = FALSE, fill = alpha("white", 0),
               linewidth = 0.7, width = 0.4) +
  ggpubr::stat_compare_means(aes(group = og_type_broad), label = "p.signif") +
  # scale_y_continuous(limits = c(0, 0.33)) +
  scale_color_manual(values = c("S. ciliatum" = "#ff9ebb",
                                "S. lacustris" = "#b9375e"),
                     guide = "none") +
  
  labs(y = "Proportion of clusters") +
  
  facet_wrap(~species_label, ncol = 2) +
  theme_bw(base_size = 12) +
  theme_for_plots +
  theme(axis.title.x = element_blank())




###############################################

correspondence_clusters_sycon <- scil_samap[[]] %>%
  select(Scil_leiden_clusters, Scil_seurat_clusters) %>%
  left_join(read.table("07_notableGenes_clusterAnnotation/cluster_identity.tsv", header = TRUE) %>%
              mutate(cluster_ID = as.factor(cluster_ID)),
            by = join_by("Scil_seurat_clusters" == "cluster_ID")) %>%
  count(Scil_leiden_clusters, Scil_seurat_clusters, value,
        name = "number_of_cells") %>%
  arrange(Scil_leiden_clusters, desc(number_of_cells)) %>% 
  rename("cell_type" = "value")
correspondence_clusters_sycon

correspondence_clusters_spongilla <- slac_samap[[]] %>%
  select(Slac_leiden_clusters, Slac_seurat_clusters_1, Slac_seurat_clusters_2) %>%
  left_join(
    data.frame(
      cell_type = c("Slac_0$" = "Archaeocytes S1", "Slac_1$" = "Archaeocytes S2",
                    "Slac_2$" = "Archaeocytes S3", "Slac_3$" = "Archaeocytes S4",
                    "Slac_4$" = "Archaeocytes S5", "Slac_5$" = "Archaeocytes S4",
                    "Slac_9$" = "Myopeptidocytes", "Slac_14$" = "Metabolocytes",
                    "Slac_15$" = "Myopeptidocytes", "Slac_17$" = "Pinacocytes",
                    "Slac_21$" = "Archaeocytes S6", "Slac_23$" = "Archaeocytes-like",
                    "Slac_24$" = "Myopeptidocytes", "Slac_25$" = "Choanocytes/-blasts",
                    "Slac_30$" = "Sclerocytes", "Slac_31$" = "Amoebocytes/Neuroid",
                    "Slac_35$" = "Basopinacocytes", "Slac_36$" = "Mesocytes",
                    "Slac_37$" = "Granulocytes-like", "Slac_41$" = "Mesocytes")
    ) %>%
      rownames_to_column() %>%
      mutate(rowname = str_remove(rowname, "\\$"),
             rowname = str_remove(rowname, "Slac_")),
    by = join_by(Slac_seurat_clusters_2 == rowname)
  ) %>%
  select(-Slac_seurat_clusters_1) %>% 
  count(Slac_leiden_clusters, Slac_seurat_clusters_2, cell_type,
        name = "number_of_cells") %>%
  arrange(Slac_leiden_clusters, desc(number_of_cells)) %>%
  replace_na(list(cell_type = "Unknown"))
correspondence_clusters_spongilla

write.table(correspondence_clusters_spongilla, "./correspondence_clusters_spongilla.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(correspondence_clusters_sycon, "./correspondence_clusters_sycon.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
