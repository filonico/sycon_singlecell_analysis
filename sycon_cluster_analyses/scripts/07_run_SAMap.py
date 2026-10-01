#!/usr/bin/env python

import subprocess

from samap import sam, SAMAP
from samap.utils import save_samap
import pandas as pd
import sys, os, pickle, argparse
from argparse import RawTextHelpFormatter


#####################
#     FUNCTIONS     #
#####################

# check if input species IDs are valid or not (based on species_metadata file)
def validate_species(species_list):
    allowed_species = ["Mlei", "Ppil", "Aque", "Scil", "Slac", "Hvul", "Nvec", "Spis", "Xesp"]
    species = sorted([s.strip() for s in species_list.split(",")])
    if not set(species).issubset(allowed_species):
        invalid_species = set(species).difference(allowed_species)
        raise argparse.ArgumentTypeError(f"Invalid species ID(s): {', '.join(invalid_species)}")
    if len(species) <= 1:
        raise argparse.ArgumentTypeError("You must select at least two species.")
    return species

# check if number of iterations is >0
def positive_int(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{value}' is not a valid integer.")
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"Number of SAMap iterations must be greater than 0, got {ivalue}.")
    return ivalue

# create samap object with samap-computed cell clusters
def create_samap_with_leiden_clustering(input):
    # assign the default resolution value to all species
    resolution_dict = {key: 3 for key in species_to_process}
    # create samap object
    samap_object = SAMAP(input,
                         f_maps = input_diamond_dir,
                         resolutions = resolution_dict,
                         save_processed = False)
    return(samap_object)


#####################
#     ARGUMENTS     #
#####################

parser = argparse.ArgumentParser(description = "Run SAMap on pre-filtered single-cell h5ad files.")
    
parser.add_argument("-s", "--species",
                    type = validate_species, default = "Aque,Scil,Slac,Hvul,Nvec,Spis,Xesp",
                    help = "Comma-separated list of species IDs (e.g., species1,species2,species3). [default: \"Aque,Scil,Slac,Hvul,Nvec,Spis,Xesp\"]")

parser.add_argument("-a", "--analysis",
                    choices = ["pairwise", "stitched", "both"], default = "both",
                    help = "Type of analysis. Choose among \"pairwise\", \"stitched\", and \"both\". [default: both]")

parser.add_argument("-i", "--input_h5ad_dir",
                    default = "./04_preprocessed_scRNAseqs/",
                    help = "Input directory where h5ad files are stored. [default: ./04_preprocessed_scRNAseqs/]")

parser.add_argument("-d", "--input_pairwise_diamond_dir",
                    default = "./03_pairwise_diamond/",
                    help = "Input directory where pairwise diamond outputs are stored. [default: ./03_pairwise_diamond/]")

parser.add_argument("-n", "--number_of_SAMap_iterarions",
                    type = positive_int, default = 6,
                    help = "Number of SAMap iterations to perform. [default: 6]")

parser.add_argument("-l", "--cluster_stitched_manifold", action = "store_true",
                    help = "Perform Leiden clustering with default parameters on the stitched manifold.")

parser.add_argument("-o", "--output_dir", required = True,
                    help = "Output directory.")


# check if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

args = parser.parse_args()

type_of_analysis_dict = {
    "pairwise": "SAMap true pairwise",
    "stitched": "SAMap stitched species",
    "both": "SAMap pairwise + stitched"
}


##########################
#     ANALYSIS RECAP     #
##########################

input_h5ad_dir = args.input_h5ad_dir
input_diamond_dir = args.input_pairwise_diamond_dir
output_dir = args.output_dir
species_to_process = args.species
number_of_SAMap_iterarions = args.number_of_SAMap_iterarions
output_suffix = 'leiden3Clusters'

# create output dir if does not exist
if not os.path.isdir(output_dir):
    subprocess.run(f"mkdir -p {output_dir}",
                   shell = True)

print("\n#####   ANALYSIS RECAP   #####")
print(f"Species to process: {', '.join(args.species)}")
print(f"Input h5ad directory: {input_h5ad_dir}")
print(f"Input diamond directory: {input_diamond_dir}")
print(f"Type of analysis: {type_of_analysis_dict[args.analysis]}")
print(f"Number of SAMap iterarions: {number_of_SAMap_iterarions}")
print(f"Output directory: {output_dir}")
print("##############################\n")

sys.stdout.flush()


########################
#     IMPORT FILES     #
########################

# define the dictionary where spIDs and path to h5ad files will be stored
data = {}

# populate the dictionary
for filepath, dirnames, filenames in os.walk(input_h5ad_dir):

    for filename in filenames:

        # load only files of chosen species
        if any(species in filename for species in species_to_process) and filename.endswith(".h5ad"):
            print("Loading raw data...")
            # extract the species ID
            spID = filename[0:4]
            # get the path to the file
            h5ad_file = os.path.join(filepath, filename)
            # perform harmony batch correction on Ppil
            
            if spID == "Scil":
                print("Preprocessing Scil with SAM (batch correction with harmonypy)...")
                sam_ppil = sam.SAM()
                sam_ppil.load_data(h5ad_file)
                sam_ppil.run(batch_key = "orig.ident")
                data[spID] = sam_ppil
            else:
                data[spID] = h5ad_file
            
            #data[spID] = h5ad_file
            

print("\n")

print(data)

# sys.exit("\nNothing more to do for now :-) we are tetsing the code. Exiting...\n")

sys.stdout.flush()

###################################
#     RUN THE SAMAP ALGORITHM     #
###################################

# create a dictionary with True booleans to be used in the neigh_from_keys SAMap algorithm
speciesTrue_dict = {key: True for key in species_to_process}

# STITCHED SPECIES ANALYSIS
if args.analysis == "stitched" or args.analysis == "both":

    # SAMAP CELL CLUSTERING
    # use the SAMap clustering method, assuming a resolution parameter of three for all species
    sam_mapping = create_samap_with_leiden_clustering(data)
    # run the SAMAP algorithm
    sam_mapping.run(pairwise = True,
                    neigh_from_keys = speciesTrue_dict,
                    n_iterations = number_of_SAMap_iterarions)

    sys.stdout.flush()

    if args.cluster_stitched_manifold:
        sam_mapping.samap.leiden_clustering(res = 4)

    # save the computed SAMAP object
    save_samap(sam_mapping, os.path.join(output_dir, "".join(species_to_process) + "_" + output_suffix + "_samap.pkl"))
    sam_mapping.samap.save_anndata(os.path.join(output_dir, "".join(species_to_process) + "_" + output_suffix + "_samap.h5ad"))

# TRUE PAIRWISE ANALYSIS
if args.analysis == "pairwise" or args.analysis == "both":
    analysis_pairs = [[x, y] for i, x in enumerate(species_to_process) for y in species_to_process[i + 1:]]

    for pair in analysis_pairs:
        # subset the data on the considered pair
        data_pairs = {species: data[species] for species in pair}

        # SAMAP CELL CLUSTERING
        # use the SAMap clustering method, assuming a resolution parameter of three for all species
        sam_mapping = create_samap_with_leiden_clustering(data_pairs)
        # run the SAMAP algorithm
        sam_mapping.run(pairwise = True,
                        neigh_from_keys = speciesTrue_dict,
                        n_iterations = number_of_SAMap_iterarions)

        if args.cluster_stitched_manifold:
            sam_mapping.samap.leiden_clustering(res = 4)

        # save the computed SAMAP object
        save_samap(sam_mapping, os.path.join(output_dir, "".join(pair) + "_" + output_suffix + "_samap.pkl"))
        sam_mapping.samap.save_anndata(os.path.join(output_dir, "".join(pair) + "_" + output_suffix + "_samap.h5ad"))

        sys.stdout.flush()
