from samap.analysis import (get_mapping_scores)

from argparse import RawTextHelpFormatter
import pandas as pd
import sys, os, subprocess, pickle, argparse


#####################
#     ARGUMENTS     #
#####################

parser = argparse.ArgumentParser(description = "Extract the SAMap mapping table from a pickle file.")
    
parser.add_argument("-p", "--input_pickle", required = True,
                    help = "Path to the pickle file to process.")

parser.add_argument("-o", "--output_dir", required = True,
                    help = "Path to the output directory.")

parser.add_argument("-n", "--n_top", default = "0",
                    help = "n_top parameter to use in the get_mapping_scores function. [default = 0]")

parser.add_argument("-s", "--species_metadata", required = True,
                    help = "Tsv file with information about single-cell annotation to use; needs to have two columns named: \"speciesID\" and \"cellCluster_annotation_name\"")

parser.add_argument("-l", "--use_leiden_clusters", action = "store_true",
                    help = "Use Leiden clusters as computed by SAMap to get mapping scores. Do not invoke this flag if want to use costum cell clusters.")


# check if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

args = parser.parse_args()


#######################
#     DEFINE I/Os     #
#######################

# define input file and output directory
output_dir = args.output_dir
input_pkl = args.input_pickle
species_metadata = pd.read_table(args.species_metadata)

# define other parameters which will be useful later on
ID = os.path.basename(input_pkl).split("_")[0]
if args.use_leiden_clusters:
    output_suffix = os.path.basename(input_pkl).split("_")[1] + "_leiden_" + args.n_top + "topCells"
else:
    output_suffix = os.path.basename(input_pkl).split("_")[1] + "_costum_" + args.n_top + "topCells"

# create output dir if does not exist
if not os.path.isdir(output_dir):
    subprocess.run(f"mkdir -p {output_dir}",
                   shell = True)
    

##############################
#     GET MAPPING SCORES     #
##############################

# load the samap obj stored in the pickle file
with open(input_pkl, "rb") as file:
    samap_obj = pickle.load(file)

# create a dictionary for the cell cluster column name
if args.use_leiden_clusters:
    cellCluster_dict = {k: "leiden_clusters" for k in species_metadata["speciesID"]}
else:
    cellCluster_dict = species_metadata.set_index("speciesID")["cellCluster_annotation_name"].to_dict()

# set the mapping score matrix name
matrix_filename = os.path.join(output_dir, ID + "_" +
                               output_suffix + "_samapMappingTable.tsv")
scoringAln_filename = os.path.join(output_dir, ID + "_" +
                                   output_suffix + "_" +
                                   args.n_top + "topCells_samapScoringAln.tsv")
# get the mapping scores
D,MappingTable = get_mapping_scores(samap_obj, cellCluster_dict, n_top = int(args.n_top))

# save to files
MappingTable.to_csv(matrix_filename, index = True, sep = '\t')
D.to_csv(scoringAln_filename, index = False, sep = '\t')