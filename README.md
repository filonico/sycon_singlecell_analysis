# The single-cell RNA-seq atlas of the calcarean sponge *Sycon ciliatum*

In this repository, you can find all the **code and scripts used to analysed the single-cell atlas of *Sycon ciliatum***, together with the re-analysis of the *Spongilla lacustris* single-cell atlas published by Musser et al. (2021; [10.1126/science.abj2949](https://doi.org/10.1126/science.abj2949)).

## Repository content
The repository is structured into three main folders:

- [`sycon_cluster_analyses/`](.sycon_cluster_analyses/) contains all the scripts and data for the annotation and molecular characterisation of the *S. ciliatum* single-cell dataset;
> [!CAUTION]
>
> POSSIBLY TO ADD
> - [`sycon_genome_annotation/`](./sycon_genome_annotation/) contains all the scripts and data for the annotation of the *S. ciliatum* genome assembly, which has been produced by the "Darwin Tree of Life" initiative;
> - [`sycon_scRNAseq_preprocessing/`](./sycon_scRNAseq_preprocessing/) contains all the scripts and data for raw-read mapping and single-cell data preprocessing of *S. ciliatum*;

- [`spongilla_remapping/`](./spongilla_remapping/) contains all the scripts and data for the re-analysis of the *S. lacustris* single-cell dataset, which consists of re-mapping raw reads to the newly realesed reference genome.

Each folder contains its own `README.md` file.