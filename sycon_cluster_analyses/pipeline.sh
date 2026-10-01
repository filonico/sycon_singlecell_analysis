#!/bin/bash

##############################
#     TIDY UP INPUT DATA     #
##############################

# tidy up raw data
# mind that data should be tar unzipped before
bash scripts/01_tidyup_data.sh


#####################
#     RUN BUSCO     #
#####################

mkdir 02_busco

# run busco on tidied proteomes
# REQUIRES: conda_envs/busco_env.yml
srun --partition=devel --job-name=busco --cpus-per-task=15 --time=02:00:00 --mem=4g --export=NONE --account=evassvis scripts/03_run_busco.sh

# summarise busco results into a tsv
grep "C:" 02_busco/*/short_summary*txt | sed -E 's/02_busco\///; s/_.+:\t/\t/; s/\t$//' > 02_busco/busco_result_summary.tsv

# compute species tree
# REQUIRES: conda_envs/phylo_env.yml
sbatch --job-name=physco --cpus-per-task=15 --time=04:00:00 --mem=4g --export=NONE --account=evassvis --mail-type=BEGIN,END,FAIL --mail-user=fn76@le.ac.uk scripts/04_compute_species_tree.sh

# gzip busco output directories
for i in 02_busco/*busco; do OUTNAME="$i".tar.gz; tar -cvf - "$i" | gzip -v -9 - > $OUTNAME && rm -rf $i; done


################################
#     RUN PAIRWISE DIAMOND     #
################################

mkdir 03_pairwise_diamond

# run all-vs-all diamond (exclude self comparison)
# REQUIRES: conda_envs/diamond_env.yml
srun --job-name=diamond --cpus-per-task=8 --time=04:00:00 --mem=0 --account=evassvis scripts/05_run_pairwise_diamond.sh


#############################
#     PREPARE H5AD FILES    #
############################

mkdir 04_preprocessed_scRNAseqs

# prepare .h5ad files
# REQUIRES: conda_envs/RSeurat_env.yml
Rscript scripts/06_prepare_h5ad_files.R


#####################
#     RUN SAMap     #
#####################

# run SAMap with the new pipeline that includes harmony correction of scil
python scripts/07_run_SAMap.py -s Aque,Scil,Slac -a both -i 04_preprocessed_scRNAseqs/ -d 03_pairwise_diamond/ -n 6 -o 05_SAMap_porifera/

# run SAMap for sponges vs sponges
# bash scripts/08_run_SAMap.sh Scil,Slac,Aque both 05b_SAMap_recodedSyconClusters | tee -a 05b_SAMap_recodedSyconClusters/Porifera_samap_both_leiden.log


################################
#     Get SAMap statistics     #
################################

# get mapping scores
for i in 05_SAMap_porifera/*pkl; do python scripts/09_get_SAMap_mappingTables.py -p $i -o 05_SAMap_porifera/01_mapping_scores -n 0; done

# get gene pairs for Placozoa and all Sycon clusters
for i in 05_SAMap_porifera/*pkl; do python scripts/10_get_SAMap_genePairs.py -p $i -o 05_SAMap_porifera/02_gene_pairs -t 0.2; done


###############################
#     SEQUENCE ANNOTATION     #
###############################

# annotate proteome with interproscan
bash /lustre/alice3/data/evassvis/fn76/SOFTWARES/InterProScan/interproscan-5.75-106.0/interproscan.sh -i 01_proteomes/Scil_ol.faa -goterms -b 09_gene_annotation/scil_proteome_interproscan


#######################################
#     GET GENES FOR GO ENRICHMENT     #
#######################################

# GO terms were annotated from the Sycon proteome with the OMA Web Server

# get the list of cluster markers and the gene universe
Rscript scripts/12_get_markers_forGOenrich.R


###################
#     hdWGCNA     #
###################

mkdir -p 12_hdWGCNA/{01_RNA_assay,02_SCT_assay}

Rscript scripts/13_hdWGCNA.R

# get GO annotation for each module
for i in 12_hdWGCNA/01_RNA_assay/*ls; do grep -wf $i 09_gene_annotation/GOterms_OMA.tsv > "${i%%.*}"_GOterms.tsv; done

# perform KO enrichment for each module
# for i in 12_hdWGCNA/01_RNA_assay/*ls; do Rscript scripts/19_perform_KEGGenrich.R 11_KEGG_enrichment/geneUniverse_KOterms.tsv $i "${i%%.*}"_KOenrich.tsv && echo done_$i; done


######################################
#     RECLUSTER THE CENTRAL BLOB     #
######################################

mkdir -p 13_recluster_blob/{01_onlyBlob_originalClusters,02_onlyBlob_newClusters,03_hdWGCNA}

# recluster the central blob and get cluster markers
Rscript scripts/14_recluster_blob.R

# get gene universe annotation
grep -wf 13_recluster_blob/sycon_onlyBlob_geneUniverse.ls 09_gene_annotation/GOterms_OMA.tsv > 13_recluster_blob/sycon_onlyBlob_geneUniverse_GOterms.tsv

# get GO annotation for each cluster marker and the gene universe
for i in 13_recluster_blob/*/*ls; do grep -wf $i 09_gene_annotation/GOterms_OMA.tsv > "${i%%.*}"_GOterms.tsv; done

# run hdWGCNA for reclusters and perform GO/KO enrich on module genes
mkdir -p 13_recluster_blob/03_hdWGCNA/01_RNA_assay/
Rscript scripts/25_hdWGCNA_blobOnly.R
for i in 13_recluster_blob/03_hdWGCNA/01_RNA_assay/*ls; do grep -wf $i 09_gene_annotation/GOterms_OMA.tsv > "${i%%.*}"_GOterms.tsv; done
# for i in 13_recluster_blob/03_hdWGCNA/01_RNA_assay/*ls; do Rscript scripts/19_perform_KEGGenrich.R 11_KEGG_enrichment/geneUniverse_KOterms.tsv $i "${i%%.*}"_KOenrich.tsv && echo done_$i; done

# run SAMap on re-clustered blob
python scripts/07_run_SAMap.py -s Scil,Aque -a pairwise -i 13_recluster_blob/04_SAMap/ -o 13_recluster_blob/04_SAMap/
python scripts/07_run_SAMap.py -s Scil,Slac -a pairwise -i 13_recluster_blob/04_SAMap/ -o 13_recluster_blob/04_SAMap/
python scripts/07_run_SAMap.py -s Scil,Aque,Slac -a pairwise -i 13_recluster_blob/04_SAMap/ -o 13_recluster_blob/04_SAMap/

# get SAMap statistics
for i in 13_recluster_blob/04_SAMap/*pkl; do python scripts/09_get_SAMap_mappingTables.py -p $i -o 13_recluster_blob/04_SAMap/01_mapping_scores -n 0; done
for i in 13_recluster_blob/04_SAMap/*pkl; do python scripts/10_get_SAMap_genePairs.py -p $i -o 13_recluster_blob/04_SAMap/02_gene_pairs -t 0.2; done


#################################
#     RUN SAMap on cnidaria     #
#################################

# run SAMap for cnidaria
bash scripts/08_run_SAMap.sh Nvec,Spis,Hvul,Xesp pairwise 15_SAMap_cnidaria/ | tee -a 15_SAMap_cnidaria/Cnidaria_samap_pairwise_leiden.log
bash scripts/08_run_SAMap.sh Aque,Slac,Scil,Nvec,Spis,Hvul,Xesp pairwise 15_SAMap_cnidaria/ | tee -a 15_SAMap_cnidaria/CnidariaPorifera_samap_pairwise_leiden.log

# get mapping scores
for i in 15_SAMap_cnidaria/*pkl; do python scripts/09_get_SAMap_mappingTables.py -p $i -o 15_SAMap_cnidaria/01_mapping_scores -n 100; done
for i in 15_SAMap_cnidaria/*pkl; do python scripts/09b_get_SAMap_mappingTables_leidenClusters.py -p $i -o 15_SAMap_cnidaria/01_mapping_scores -n 100; done

# get gene pairs
for i in 15_SAMap_cnidaria/*pkl; do python scripts/10_get_SAMap_genePairs.py -p $i -o 15_SAMap_cnidaria/02_gene_pairs -t 0.2; done