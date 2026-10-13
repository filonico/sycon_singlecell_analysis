#!/usr/bin/env Rscript

suppressMessages({library(tidyverse)
  library(topGO)})

args = commandArgs(trailingOnly = TRUE)

# args1 gene universe
# args2 genes of interest
# args3 output prefix


######################
#     READ INPUT     #
######################

universe_GOannotation <- read.table(file = args[1],
                                    header = FALSE,
                                    sep = "\t",
                                    na.strings = "-")
interest_GOannotation <- read.table(file = args[2],
                                    header = FALSE,
                                    sep = "\t",
                                    na.strings = "-")


#####################
#     FUNCTIONS     #
#####################

# function to perform GO enrichment and generate summary tables
performGOEnrichment <- function(GOannotation, ontology, algorithm) {
  
  # create topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = GOannotation)
  
  # perform GO enrichment
  result <- runTest(GOdata, algorithm = algorithm, statistic = "fisher")
  
  # get summary table
  results <- GenTable(GOdata,
                      classicFisher = result,
                      ranksOf = "classicFisher",
                      topNodes = 500)
  
  return(results)
}


############################
#     DATA PREPARATION     #
############################

# prepare the gene universe dataset as a named list
# this is with EggNog Mapper outputs
geneID2GO <- universe_GOannotation %>%
  dplyr::select(c(V1, V10)) %>%
  drop_na() %>%
  mutate(V10 = strsplit(V10, ",")) %>%
  tibble::deframe()

# # prepare the gene universe dataset as a named list
# # this is with OMA outputs
# geneID2GO <- universe_GOannotation %>%
#   dplyr::select(c(V2, V5)) %>%
#   group_by(V2) %>%
#   summarize(V5 = paste(V5, collapse = ","), .groups = 'drop') %>%
#   mutate(V5 = strsplit(V5, ",")) %>%
#   tibble::deframe()
  
# get list of genes
geneNames <- names(geneID2GO)

# get the vector of the genes of interest
# this is with EggNog Mapper outputs
gene_int_list <- as.vector(interest_GOannotation$V1)

# # get the vector of the genes of interest
# # this is with OMA outputs
# gene_int_list <- as.vector(interest_GOannotation$V2)

# get a factor list of interesting/not interesting genes 
geneList <- factor(as.integer(geneNames %in% gene_int_list))
names(geneList) <- geneNames


#########################
#     GO ENRICHMENT     #
#########################

# define GO ontologies
ontologies <- c("BP", "MF", "CC")

# define algorithms
algorithms <- c("classic", "elim")

# perform GO enrichment and generate summary tables for each ontology and each method
for (ontology in ontologies) {
  for (algorithm in algorithms) {
    
    result <- performGOEnrichment(geneID2GO,
                                  ontology,
                                  algorithm)
    file_name <- paste0(args[3],
                        "topGO_", ontology, "_", algorithm, ".txt")
    write.table(result,
                file = file_name,
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")
  }
}
