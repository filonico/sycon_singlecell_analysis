from samap.analysis import (GenePairFinder)

from argparse import RawTextHelpFormatter
import pandas as pd
import sys, os, subprocess, pickle, argparse


#####################
#     ARGUMENTS     #
#####################

parser = argparse.ArgumentParser(description = "Extract the gene pairs from SAMap results. NOTE that this script will work only with species pairs at the moment.")
    
parser.add_argument("-p", "--input_pickle", required = True,
                    help = "Path to the pickle file to process.")

parser.add_argument("-o", "--output_dir", required = True,
                    help = "Path to the output directory.")

parser.add_argument("-n", "--n_top", default = "0",
                    help = "n_top parameter to use in the get_mapping_scores function. [default = 0]")

parser.add_argument("-c", "--cell_clusters", default = "all",
                    help = "Name of the cell cluster of which to get gene pairs. Pprepend the species ID to the name (e.g., \"Scil_pinacocytes\"). [default = all]")

parser.add_argument("-t", "--threshold", default = 0.4, type = float, 
                    help = "Mapping score threshold. [default = 0.4]")

parser.add_argument("-l", "--use_leiden_clusters", action = "store_true", 
                    help = "Include this flag to use Leiden clusters to compute gene pairs. [default = False]")


# check if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

args = parser.parse_args()


#####################
#     FUNCTIONS     #
#####################

# def add_speciesID_to_cluster_names(species_id):
#     obs = samap_obj.samap.adata.obs.copy()
    
#     # define column names
#     cluster_col = f"{species_id}_N1.Int.Clusters"   # for sycon
#     cell_type_col = f"{species_id}_cell_type"       # for any other species
    
#     # define the criterion to match only corresponding lines
#     mask = obs['species'] == species_id
    
#     # add prefixes
#     if cluster_col in obs.columns:
#         obs.loc[mask, cluster_col] = (
#             species_id + '_' + obs.loc[mask, cluster_col].astype(str)
#         )
    
#     if cell_type_col in obs.columns:
#         obs.loc[mask, cell_type_col] = (
#             species_id + '_' + obs.loc[mask, cell_type_col].astype(str)
#         )
    
#     # update the original dataframe
#     samap_obj.samap.adata.obs = obs


#######################
#     DEFINE I/Os     #
#######################

# define input file and output directory
output_dir = args.output_dir
input_pkl = args.input_pickle
cell_cluster = args.cell_clusters
threshold = args.threshold

# define other parameters which will be useful later on
ID = os.path.basename(input_pkl).split("_")[0]              # ID of the SAMap run
species = [ID[:4], ID[-4:]]                                 # included species in the analysis
if args.use_leiden_clusters:
    output_suffix = os.path.basename(input_pkl).split("_")[1] + "_leiden_" + \
        str(threshold).replace(".", "") + "threhsold_" + \
        args.n_top + "topCells_" + cell_cluster.replace("_", "")
else:
    output_suffix = os.path.basename(input_pkl).split("_")[1] + "_costum_" + \
        str(threshold).replace(".", "") + "threhsold_" + \
        args.n_top + "topCells_" + cell_cluster.replace("_", "")

# load the table with species metadata
species_metadata = pd.read_table("00_input/species_metadata.tsv")

# create output dir if does not exist
if not os.path.isdir(output_dir):
    subprocess.run(f"mkdir -p {output_dir}",
                   shell = True)
    

##########################
#     GET GENE PAIRS     #
##########################

# load the samap obj stored in the pickle file
with open(input_pkl, "rb") as file:
    samap_obj = pickle.load(file)

# add prefixes for the two species
# for sp in species:
#     add_speciesID_to_cluster_names(sp)

if args.use_leiden_clusters:
    # create a dictionary for the cell cluster column name
    leidCluster_dict = dict.fromkeys(species_metadata["speciesID"], "leiden_clusters")

    # find markers on Leiden clusters
    gpf = GenePairFinder(samap_obj, keys = leidCluster_dict)
else:
    # create a dictionary for the cell cluster column name
    cellCluster_dict = species_metadata.set_index("speciesID")["cellCluster_annotation_name"].to_dict()

    # find markers on pre-computed clusters
    gpf = GenePairFinder(samap_obj, keys = cellCluster_dict)

# find gene pairs
if cell_cluster == "all":
    gene_pairs = gpf.find_all(align_thr = threshold, n_top = int(args.n_top))
else:
    gene_pairs = gpf.find_all(n = cell_cluster, align_thr = threshold, n_top = int(args.n_top))

# set the gene pair matrix output file name
matrix_gene_pairs = os.path.join(output_dir, ID + "_" +
                                 output_suffix + "Cluster_samapGenePairs.tsv")

# save the gene pair matrix to a file
gene_pairs.to_csv(matrix_gene_pairs, index = False, sep = '\t')