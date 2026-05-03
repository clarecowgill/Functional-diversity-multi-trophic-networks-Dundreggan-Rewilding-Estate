# Dundreggan Multi-Trophic Biodiversity Dataset and Analysis

## Overview

This repository contains all data and code required to reproduce analyses from:

**“Integrated monitoring reveals recovery of functional diversity and multi-trophic networks across a woodland restoration chronosequence”**

The study integrates **environmental DNA (fungi, invertebrates)** and **vertebrate monitoring (acoustic + tree sampling)** to quantify biodiversity patterns and reconstruct **multi-trophic interaction networks (metawebs)**.

A key contribution of this repository is a **bespoke functional trait database for soil invertebrates**, developed to enable ecological inference from 18S metabarcoding data where trait information is otherwise sparse.

---

## 1. Data Files

### Fungi

- `OTU_table_ITS2.tsv`  
  ITS2 OTU table (97% clustered sequences)

- `soil_meta.csv`  
  Sample metadata

- `wca.csv`  
  Woodland condition metrics

---

### Invertebrates (18S)

- `18S_ASV_table_with_taxonomy.tsv`  
  Final ASV table with taxonomic assignments and confidence scores

- `dundreggan_invert_taxonomy_refine.csv`  
  Curated taxonomy table following manual refinement

- **`invertebrate_trait_database.xlsx`**  
  Bespoke functional trait database for soil invertebrates.

  This dataset assigns ecological traits (feeding guild, microhabitat/ soil strata, moisture affinity and body length) to invertebrate taxa detected via 18S metabarcoding, including family-level and higher assignments where species-level traits could not be found.

  This database was developed due to a lack of existing trait resources for soil invertebrates and is intended as a reusable resource for the wider community.

- `soil_meta.csv`  
  Sample metadata

- `wca.csv`  
  Woodland condition metrics

---

### Vertebrates

- `dundreggan_matrix.csv`  
  Tree roller eDNA OTU read table

  - `all_acoustics.csv`  
  Acoustic detections (BirdNET outputs)

- `blank_association.csv`  
  Filtering/association reference

- `dun_meta.csv`  
  Metadata for vertebrate sampling

- `vert_matrix.csv`  
  Species detection matrix

- `vert_traitdata.csv`  
  Vertebrate trait data (diet, foraging strata, body size)

- `verts_nvdi_stats.csv`  
  NDVI-derived habitat metrics

- `wca.csv`  
  Woodland condition metrics

---

### Metaweb Analysis

Pre-processed datasets used to construct interaction networks:

- `acoustic_pa_long.csv` – vertebrate presence–absence (long format)  
- `invert_pa_long.csv` – invertebrate presence–absence  
- `fungi_genus_long.csv` – fungal presence–absence  
- `tree_pa_long.csv` – tree species presence  

Trait datasets:

- `invert_traits_networks.csv`  
- `fungi_traits_networks.csv`  
- `verts_traitdata.csv`  

---

## 2. R Scripts

### `fungi analysis.R`
Processes ITS2 fungal OTU data and calculates all diversity metrics and functional space analysis.

### `invertebrate taxonomy refinement.R`
Applies manual curation and standardisation of 18S taxonomic assignments, intended to be run prior to invertebrate analysis.

### `invertebrate analysis.R`
Processes 18S invertebrate ASV data, integrates taxonomy and traits, and derives community-level metrics.

### `vertebrate analysis.R`
Processes acoustic and tree-derived vertebrate detections, including analysis of functional diveristy.

### `MetaWebs.R`
Constructs multi-trophic interaction networks (metawebs) using:
- presence–absence data  
- trait-based interaction rules  
- inferred trophic links
Note that networks were created for each soil sample (for co-occurring soil fauna), but importantly vertebrate taxa were incorporated at the site level (i.e. the same across all networks from the same site)

Visualises example networks and calculates network metrics including:
- species richness (S)  
- number of links (L)  
- connectance  
- modularity  
- robustness (R50)

---

## 3. Taxonomic Assignment

### ITS2 (Fungi)

- `ITS2_taxonomic_assignment.sh`  
  QIIME2-based pipeline using the **UNITE database (v2025-02-19, 97% clustering)**.

---

### 18S (Invertebrates)

- `18S_taxonomic_assignment.sh`  
  QIIME2 pipeline using **SILVA 138.1**, including:
  - primer trimming  
  - classifier training  
  - taxonomic assignment  

---

### 12S (Vertebrates)

- `12s_verts.fasta`  
  Reference database of sequences for UK vertebrate taxa

- `12s_verts_tax_map.txt`  
  Mapping file linking sequences to taxonomic IDs (GenBank)

These files were used for taxonomic identification of UK vertebrates from 12S data.

---

## Bioinformatic Processing Notes

- **12S (vertebrates):** Processed using Tapirs pipeline, including quality filtering, merging, clustering (VSEARCH), and BLAST-based taxonomic assignment against a curated UK vertebrate database.

- **18S (invertebrates) and ITS2 (fungi):** Processed using **DADA2**, including:
  - quality filtering  
  - error correction  
  - chimera removal  

  - 18S data were retained as **ASVs**  
  - ITS2 data were clustered to **97% OTUs** prior to QIIME2 classification  

- Taxonomic assignment was performed in **QIIME2 (v2025.10)** using:
  - **SILVA (v138.1)** for 18S  
  - **UNITE (v2025-02-19)** for ITS2  

---

## Reproducibility

To reproduce analyses:

1. Run taxonomic assignment scripts in `Taxonomic assignment/`  
2. Use processed data in `Data files/`  
3. Run R scripts in the following order:
   - `invertebrate taxonomy refinement.R`  
   - `fungi analysis.R`  
   - `invertebrate analysis.R`  
   - `vertebrate analysis.R`  
   - `MetaWebs.R`  

---
