# Genome decontamination pipeline

<img src="bacteria_logo.png" alt="Logo" width="200" height="160"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-brightgreen)

# Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Example usage](#example-usage)
- [Detailed options](#detailed-options)
- [Output explanation](#output-explanation)
- [Citation](#citation)
- [Troubleshooting](#troubleshooting)
- [Changelog](#changelog)

# Introduction

This repository provides a reproducible genome decontamination pipeline developed for large-scale analysis of sponge genome assemblies and designed to be broadly applicable to other non-model organisms

The pipeline combines **three complementary lines of evidence** to identify contaminant scaffolds:

- **Compositional score** — based on sequence composition. A scaffold is flagged when either k-mer profile clustering or GC-content outlier detection indicates contamination
- **Protein taxonomy score** — based on protein similarity searches. A scaffold is flagged when either the DIAMOND consensus method or MEGAN taxonomic assignment indicates a non-target origin
- **Nucleotide taxonomy score** — based on scaffold k-mer classification. A scaffold is flagged when k-mers classified by Kraken2 suggest a non-target taxon

A scaffold is classified as a contaminant only when supported by **at least two independent signals**, providing a consensus-based decision framework that balances precision and recall.

This pipeline is intended for researchers working with draft genome assemblies, particularly for non-model organisms or systems where contamination from symbionts or environmental DNA is difficult to avoid. It provides an automated framework for generating more reliable genome assemblies. Notably, the pipeline can only detect contamination originating from different taxonomic kingdoms (e.g., Metazoa vs. Bacteria, Metazoa vs. Fungi...).

The pipeline produces a comprehensive set of output files to facilitate downstream analysis and reproducibility. These include:

- **Detailed contamination statistics** summarizing contamination evidence per scaffold
- A **taxonomy report of detected contaminants**, including taxonomic assignments for contaminant scaffolds
- A **FASTA file containing the decontaminated genome assembly**
- A **FASTA file of contaminant scaffolds**

More information on the pipeline methodology and benchmarking can be found im the preprint (coming soon)

# Installation

The decontamination pipeline consists of a Bash and an R script located in the `scripts` directory of this repository. After cloning the repository. both scripts must be included in the `PATH` environment variable.
The following dependencies are required and should be included in `PATH`:

| **Dependency** | **Tested Version** | **Homepage**                           | **Conda Installation**              |
|----------------|--------------------|----------------------------------------|-------------------------------------|
| R              | 4.5.2              | https://www.r-project.org              | `conda install r-base`              |
| DIAMOND        | 2.1.13             | https://github.com/bbuchfink/diamond   | `conda install -c bioconda diamond` |
| MEGAN          | 7.1.1              | https://github.com/husonlab/megan-ce   | `conda install -c bioconda megan  ` |
| Kraken2        | 2.1.6              | https://github.com/DerrickWood/kraken2 | `conda install -c bioconda kraken2` |

The following R packages are also required:

| **R package**            | **Tested Version** | **Homepage**                                                       | **R installation**                      |
|--------------------------|--------------------|--------------------------------------------------------------------|-----------------------------------------|
| data.table               | 1.17.8             | https://cran.r-project.org/package=data.table                      | `install.packages("data.table")`        |
| taxonomizr               | 0.11.1             | https://cran.r-project.org/web/packages/taxonomizr                 | `install.packages("taxonomizr")`        |
| Biostrings               | 2.76.0             | https://bioconductor.org/packages/release/bioc/html/Biostrings/    | `BiocManager::install("Biostrings")`    |
| compositions             | 2.0.9              | https://cran.r-project.org/web/packages/compositions               | `install.packages("compositions")`      |
| mclust                   | 6.1.2              | https://cran.r-project.org/web/packages/mclust                     | `install.packages("mclust")`            |
| umap                     | 0.2.10             | https://cran.r-project.org/web/packages/umap                       | `install.packages("umap")`              |
| ggplot2                  | 4.0.2              | https://cran.r-project.org/web/packages/ggplot2                    | `install.packages("ggplot2")`           |
| GenomicRanges            | 1.62.1             | https://bioconductor.org/packages/release/bioc/html/GenomicRanges/ | `BiocManager::install("GenomicRanges")` |

# Example usage

The genome decontamination pipeline requires two mandatory arguments:

1. Genome assembly in FASTA format
2. Expected (target) kingdom-level taxon of the analysed genome assembly

The following kingdom-level taxons are supported:

Viruses, Bacteria, Archaea, Metazoa, Fungi, Viridiplantae, Haptophyta, Rhodophyta, Stramenopiles, Alveolata, Excavata, Rhizaria, Choanoflagellata, and Amoebozoa


Example usage:

```bash
decontamination_pipeline [OPTIONS] GENOME TARGET_TAXON
```

# Detailed options

The decontamination pipeline offers a comprehensive list of options which allow users to control the analysis parameters

## Supplementary files

If available, users can supply several preprocessed files to reduce the pipeline runtime

These files can be supplied with the following flags:

`-n`: Pre-generated DIAMOND NR database

`-d`: Pre-generated unmeganized DIAMOND results

`-m`: Pre-generated meganized DIAMOND results

`-k`: Pre-generated Kraken2 core NT database index

`-K`: Pre-generated Kraken2 results

`-T`: Pre-generated Taxonomizr database

If database files are not provided, the pipeline will download the current database versions from NCBI (requires an active internet connection)

If DIAMOND results are provided by the user, the NR database does not need to be supplied. Similarly, if Kraken2 results are provided, the Kraken2 database index is not required

For best performance, the provided DIAMOND output files should be obtained by using the same parameters employed in this pipeline. These include:

| Tool     | Parameters / Database |
|----------|-----------------------------------------------------------------|
| DIAMOND  | `-e 1e-5 --range-culling -F 15 --top 10`; database: **NCBI NR** |
| MEGAN    | `-lg -alg longReads -me 1e-5`; mapping file: **MEGAN7**         |
| Kraken2  | Default parameters; database: **NCBI core NT**                  |

## DIAMOND run options

The pipeline provides the following option to control DIAMOND speed and CPU/RAM usage:

`-f`: Fast DIAMOND parameter preset, default: off

By default, DIAMOND runs with the `--sensitive` parameter preset. Enabling the fast preset speeds up the run at the cost of reduced sensitivity

`-g`: Number of splits performed on the genome, default: 1

`-s`: Number of splits performed on the DIAMOND NR database, default: 2

Users can increase the number of splits in the query (genome) and target (NR database) to reduce RAM usage, at the cost of longer runtimes

## Contamination analysis options

`-S`: Random seed for reproducibility, default: 100

Random seed is defined to ensure reproducible k-mer-based clustering

`-L`: Minimum scaffold length for taxonomy determination (in bp), default: 1000

Length threshold is required as very short scaffolds cannot be reliably classified using k-mer profiles or taxonomic signals. Scaffolds shorter than this threshold are therefore excluded from the analysis and the final decontaminated assembly

`-G`: Maximum number of Gaussian mixture clusters used for sequence k-mer clustering, default: 10

This parameter defines the maximum number of scaffold clusters with distinct k-mer profiles. The value may be increased for assemblies where a high level of contamination is expected

`-C`: Significance level of the chi-square distribution used to filter k-mer clusters by Mahalanobis distance, default: 0.005

Scaffold clusters are classified as contaminants based on the Mahalanobis distance between their centroid and the reference centroid representing the target genome. This parameter defines the significance threshold used to determine whether a cluster is considered an outlier

`-Z`: Z-score threshold for GC-content outlier classification, default: 2.5

Scaffolds are classified as GC-content outliers if their GC content deviates significantly from the genome mean. This parameter defines the Z-score threshold used to determine whether a scaffold is considered an outlier

`-W`: Multiplier for the DIAMOND hit set weight comparison; a scaffold is flagged as a contaminant if non-target weight > W × target weight (must be satisfied together with `-c`), default: 1.2

`-c`: Minimum non-target DIAMOND hit set coverage fraction of a scaffold required for contaminant classification (must be satisfied together with `-W`), default: 0.2

Overlapping DIAMOND hits on each scaffold are merged into hit sets. The weight of each hit set is proportional to the mean bitscore of the DIAMOND hits and the length of the hit set
The DIAMOND consensus method assigns a scaffold a non-target taxonomic origin if the total weight of a non-target taxon is significantly higher than the weight of the target taxon, and if the coverage of the non-target hit sets exceeds the required threshold

`-R`: Multiplier for the Kraken2 k-mer fraction comparison; a scaffold is flagged as a contaminant if non-target fraction > R × target fraction (must be satisfied together with `-r`), default: 1.2

`-r`: Minimum non-target Kraken2 k-mer fraction required for contaminant classification (must be satisfied together with `-R`), default: 0.05

The nucleotide-based methodology assigns a scaffold a non-target taxonomic origin if the total k-mer fraction of a non-target taxon is significantly higher than that of the target taxon, and if the non-target k-mer fraction exceeds the required threshold.

`-A`: Minimum Kraken2 k-mer fraction required to assign a taxon to a scaffold, default: 0.1

`-a`: Minimum Kraken2 relative k-mer fraction required to assign a taxon to a scaffold (normalized by the total number of classified k-mers per scaffold), default: 0.5

`-B`: Minimum DIAMOND hit set coverage fraction required to assign a taxon to a scaffold, default: 0.25

`-b`: Minimum DIAMOND hit set relative coverage fraction required to assign a taxon to a scaffold (normalized by the total length of DIAMOND hit sets per scaffold), default: 0.5

Contaminant scaffolds are assigned taxonomy with the highest possible resolution (maximum rank: genus). If the nucleotide-based method suggests contamination, the scaffold is assigned the most resolved taxon whose k-mer fraction satisfies the thresholds defined above. If nucleotide evidence is missing but the protein-based method suggests contamination, the scaffold is assigned the most resolved fraction for which DIAMOND hit set coverage meets the required thresholds. If none of the criteria are met at any taxonomic rank, the scaffold taxonomy is set to unknown.

## General options

`-t`: Number of CPU threads, default: 10

Several steps of the decontamination pipeline are parallelized. Recommended number of threads: 10-20.

`-D`: Output directory name, default: GENOME_decontamination_dir

`-O`: Overwrite the output directory, default: off

`-h`: Show usage information

# Output explanation

The pipeline produces several output files outlined below:

`scaffold_contamination_info.tsv`: Contains the main contamination statistics for each scaffold:

 | **Column**                                | **Description**                                                                                                 |
 |-------------------------------------------|-----------------------------------------------------------------------------------------------------------------|
 | `scaffold_name`                           | Scaffold name                                                                                                   |
 | `scaffold_len`                            | Scaffold length (in bp)                                                                                         |
 | `cluster`                                 | K-mer-based cluster assigned by mclust                                                                          |
 | `outlier_kmer`                            | Binary indicator of k-mer profile outliers                                                                      |
 | `gc_content`                              | Scaffold GC content                                                                                             |
 | `gc_z`                                    | Z-score representing the deviation of GC content from the genome mean                                           |
 | `outlier_gc`                              | Binary indicator of GC content outliers                                                                         |
 | `dia_hit_set_weight_ratio_target`         | Total DIAMOND hit set weigth ratio of the target taxon                                                          |
 | `dia_hit_set_cov_fraction_target`         | Fraction of scaffold length contained in DIAMOND hit sets of the target taxon                                   |
 | `rel_dia_hit_set_cov_fraction_target`     | Relative length fraction of DIAMOND hit sets of the target taxon (compared to all hit sets)                     |
 | `dia_hit_set_weight_ratio_non_target`     | Total DIAMOND hit set weigth of highest-weigth non-target taxon                                                 |
 | `dia_hit_set_cov_fraction_non_target`     | Fraction of scaffold length contained in DIAMOND hit sets of the highest-weigth non-target taxon                |
 | `rel_dia_hit_set_cov_fraction_non_target` | Relative length fraction of DIAMOND hit sets of the highest-weigth non-target taxon (compared to all hit sets)  |
 | `outlier_diamond`                         | Binary indicator of DIAMOND consensus method outliers                                                           |
 | `taxon_lca`                               | Taxon assigned by MEGAN LCA                                                                                     |
 | `outlier_lca`                             | Binary indicator of MEGAN LCA outliers                                                                          |
 | `kraken_kmer_fraction_target`             | K-mer fraction of the target taxon (compared to all k-mers)                                                     |
 | `kraken_classif_kmer_fraction_target`     | K-mer fraction of the target taxon (compared to taxonomically classified k-mers)                                |
 | `kraken_kmer_fraction_non-target`         | K-mer fraction of the highest-fraction non-target taxon (compared to all k-mers)                                |
 | `kraken_classif_kmer_fractiob_non-target` | K-mer fraction of the highest-fraction non-target taxon (compared to taxonomically classified k-mers)           |
 | `outlier_kraken`                          | Binary indicator of Kraken2 outliers                                                                            |
 | `compositional_score`                     | Compositional score                                                                                             |
 | `prot_taxonomy_score`                     | Protein taxonomy score                                                                                          |
 | `nucl_taxonomy_score`                     | Nucleotide taxonomy score                                                                                       |
 | `contamination_score`                     | Contamination score                                                                                             |
 | `contamination`                           | Binary indicator of contamination                                                                               |


 `contaminant_scaffolds_taxonomy.tsv`: Contains the taxonomic classification of contaminant scaffolds:

 | **Column**                                 | **Description**                                                                                 |
 |--------------------------------------------|-------------------------------------------------------------------------------------------------|
 | `scaffold_name`                            | Scaffold name                                                                                   |
 | `scaffold_len`                             | Scaffold length (in bp)                                                                         |
 | `dia_hit_set_cov_fraction_taxon_rank`*     | Fraction of scaffold length contained in DIAMOND hit sets of the taxon rank                     |
 | `rel_dia_hit_set_cov_fraction_taxon_rank`* | Relative length fraction of DIAMOND hit sets of the taxon rank (compared to all other hit sets) |
 | `taxon_rank_prot`*                         | Protein-based taxon assignment                                                                  |
 | `kraken_kmer_fraction_tax_rank`*           | K-mer fraction of the taxon rank (compared to all k-mers)                                       |
 | `kraken_classif_kmer_fraction_tax_rank`*   | K-mer fraction of the taxon rank (compared to taxonomically classified k-mers)                  |
 | `taxon_rank_nucl`*                         | Nucleotide-based taxon assignment                                                               |
 | `assigned_tax`                             | Final taxonomy call                                                                             |
 * Appleis to multiple columns for the following taxonomic ranks: genus, family, order, class, phylum, domain

 `decontaminated_genome.fa`: FASTA file of the decontaminated genome assembly

 `contaminant_sequences.fa`: FASTA file of contaminant sequences

 `mc_umap_dt.tsv`: UMAP projection coordinates for principal components derived from scaffold k-mer profiles

 `mc_umap_plot.tiff`: UMAP projection plot for principal components derived from scaffold k-mer profiles (colored by mclust clusters)

 # Citation

 The genome decontamination pipeline is distributed under the MIT license.

 Copyright © 2026 Kristian Bodulić

 if you use the genome contamination pipeline, please cite the following preprint: (coming soon)

 # Troubleshooting

 Please report all potential bugs in the Issues tracker.

 # Changelog

 Version 1.0.0: Bug fixes regarding supplementary file simbolic linking, March 20, 2026.
 Version 1.0.0: Initial commit, March 16, 2026.
