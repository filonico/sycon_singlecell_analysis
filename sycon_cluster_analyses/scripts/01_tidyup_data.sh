#!/bin/bash

#####################################
#     TIDY UP S. LACUSTRIS DATA     #
#####################################

# proteome taken from Matt
# remove annotation from gene IDs in UMI matrix
zcat 00_input/scRNAseqs_rawTables/S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix.txt.gz | sed -E 's/ [^\t]+//' | gzip -v -9 - > 00_input/scRNAseqs_rawTables/S_lacustris/GSE134912_Slac_spongilla_10x_count_matrix_edited.txt.gz

# edit gene IDs in the proteome fasta
zcat 00_input/genomic_data/S_lacustris/Spongilla.fasta.gz | sed -E 's/-/_/' > 01_proteomes/Slac.faa


#############################################
#     TIDY UP A. QUEENSLANDICA PROTEOME     #
#############################################

# proteome taken from Matt
# edit gene IDs in the proteome fasta
zcat 00_input/genomic_data/A_queenslandica/Aque.fasta.gz | sed -E 's/-/_/' > 01_proteomes/Aque.faa


##########################################
#     ONELINE AND REMOVE UNDERSCORES     #
##########################################

# underscores are replaced with "-" because Seurat does not like them
for i in 01_proteomes/*faa; do FILENAME="${i%.*}"; FILEEXTENSION="${i##*.}"; awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | tail -n +2 | sed -E 's/_/-/g' > "$FILENAME"_ol."$FILEEXTENSION"; done
