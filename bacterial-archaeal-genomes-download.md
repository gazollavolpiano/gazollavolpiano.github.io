# How to Easily Download Public Bacterial and Archaeal Genomes

## 🚧 Problem
Downloading genomes can be very frustrating when navigating back and forth through the [NCBI Assembly Database](https://www.ncbi.nlm.nih.gov/assembly/) website. 

## 💡 What You Will Learn
- **Selecting Genomes**: Filter based on taxonomy as defined by the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/), quality, isolation source, or other relevant metadata.
- **Downloading Files**: Use a download tool to obtain the selected genomes.

## 💾 Storage Requirements

Rough estimates of storage requirements for  (`.fna`) genomes:

- **100 Genomes**: ≈250 MB
- **1,000 Genomes**: ≈2.5 GB
- **Whole GTDB r220 (~600K Genomes)**: ≈1.5 TB

## 🧬 Why GTDB Labels?

The **Genome Taxonomy Database (GTDB)** offers a standardized and up-to-date taxonomy. Here's why using it is beneficial:

- **Consistency**: GTDB provides a uniform naming convention, minimizing confusion from outdated/conflicting classifications.
- **Accuracy**: Regularly updated with the latest sequences, including draft genomes from metagenomes and isolates.
- **Interoperability**: Facilitates data comparison and integration across various studies.
- **Quality Control**: All genomes undergo quality checks using CheckM.

For more details, visit [GTDB Methods](https://gtdb.ecogenomic.org/methods).

## Selecting Genomes

To start our process, we are first going to download a table with the GTDB metadata for the desired GTDB release. In this guide, I will be working with **r220**.

```bash
# Download the GTDB metadata for bacteria (r220) and decompress it:
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz
gunzip bac120_metadata_r220.tsv.gz

# If you need Archaea genomes:
# wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/ar53_metadata_r220.tsv.gz
# gunzip ar53_metadata_r220.tsv.gz
```
### Filtering the Metadata

Now, you can use the data in `bac120_metadata_r220.tsv` (or `ar53_metadata_r220.tsv` if analyzing Archaea) to filter according to desired criteria. Some of the most common criteria include:

- **`gtdb_taxonomy`**: Use this column to select genomes from specific GTDB labels.

- **`gtdb_type_designation_ncbi_taxa`**: This can be used to select genomes from a type strain (strains that serves as the reference for that species' characteristics).

- **`ncbi_assembly_level`**: Use this to select the assembly level (e.g., Complete Genome, Contig, etc.).

- **`checkm2_completeness`** and **`checkm2_contamination`**: These columns indicate the completeness and contamination of the genome assembly. GTDB already filters for quality (completeness estimate >50%, contamination estimate <10%; see [GTDB Methods](https://gtdb.ecogenomic.org/methods) for more details), but additional filtering can help reduce the size of your database.

- **`contig_count`**: The number of contigs in the assembly. GTDB already filters for assemblies with <1,000 contigs, but this can be useful for further refining your dataset.

- **`ncbi_country`** and **`ncbi_isolation_source`**: These columns provide context about the origin of the genome, which can be important for ecological or epidemiological studies.

📝 There are several other columns that may be relevant, make sure to check them all!

#### Example: Downloading *Staphylococcus aureus* genomes that are complete, from Australia, and isolated from blood.

Use R:

```R
# Load necessary library
library(tidyverse)

# Read GTDB metadata
bac120 <- read_tsv("bac120_metadata_r220.tsv")
dim(bac120) # 584,382 genomes with 113 columns of information

# Inspect the first two rows of the gtdb_taxonomy column
head(bac120$gtdb_taxonomy, 2)
# Example output:
# [1] "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Bordetella;s__Bordetella pseudohinzii"
# [2] "d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus epidermidis"        

# Search for Staphylococcus aureus in the gtdb_taxonomy column
staph_aureus_taxonomy <- grep("Staphylococcus aureus", bac120$gtdb_taxonomy, value = TRUE)[1]
staph_aureus_taxonomy
# [1] "d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus"

# Filter for Staphylococcus aureus genomes
staph_aureus <- bac120 %>%
  filter(gtdb_taxonomy == staph_aureus_taxonomy)
dim(staph_aureus) # 16,021 Staphylococcus aureus genomes in metadata

# Further filter by completeness (Complete Genome)
staph_aureus <- staph_aureus %>%
  filter(ncbi_assembly_level == "Complete Genome")
dim(staph_aureus) # 1,106 complete Staphylococcus aureus genomes

# Filter by country: Australia
staph_aureus <- staph_aureus %>%
  filter(str_detect(ncbi_country, "Australia"))
dim(staph_aureus) # 28 complete genomes from Australia

# Inspect isolation_source and filter for blood-related sources
table(staph_aureus$ncbi_isolation_source)
# Possible output categories like blood, Blood Culture, bloodstream, CA-MRSA blood site, etc.

staph_aureus <- staph_aureus %>%
  filter(str_detect(ncbi_isolation_source, regex("blood", ignore_case = TRUE)))
dim(staph_aureus) # 12 complete genomes from Australia isolated from blood

# Write the filtered metadata to a TSV file
write_tsv(staph_aureus, "staphaureus_12_complete_genomes_from_australia_isolated_from_blood.tsv")

# Save the assembly_accession column to a separate file for downloading genomes
staph_aureus %>%
  pull(ncbi_genbank_assembly_accession) %>%
  write_lines("assembly_accessions_12_genomes.txt")
```

## Downloading Genomes with NCBI Datasets 

Use the **NCBI Datasets** tool to download 12 genomes.

```bash
# Install and activate ncbi datasets tool with micromamba
micromamba create -n ncbi_datasets -c conda-forge ncbi-datasets-cli 
micromamba activate ncbi_datasets 
datasets --version # datasets version: 16.31.0

# Create directory to store the genomes
mkdir genomes_data

# Inspect file with genome accessions
head assembly_accessions_12_genomes.txt -n 3
#GCA_900607275.1
#GCA_013414865.1
#GCA_900607295.1

# Download the genomes
cat assembly_accessions_12_genomes.txt | while read accession; do  
  datasets download genome accession "$accession" --filename "genomes_data/${accession}.zip" 
  unzip -o "genomes_data/${accession}.zip" -d genomes_data 
  mv genomes_data/ncbi_dataset/data/*/*.fna genomes_data/ 
  rm -rf "genomes_data/${accession}.zip" genomes_data/*md genomes_data/*txt genomes_data/ncbi_dataset
done 

# Count the number of genomes
ls genomes_data/*.fna | wc -l # 12, as expected
```

Now you can use the `.fna` files for your downstream analysis. The metadata file (`staphaureus_12_complete_genomes_from_australia_isolated_from_blood.tsv`) contains detailed information about each genome, which can be used for further analysis. 🎉







