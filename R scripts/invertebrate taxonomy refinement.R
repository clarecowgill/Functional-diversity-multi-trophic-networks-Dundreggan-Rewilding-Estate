library(dplyr)
library(stringr)
library(taxize)
library(purrr)

## ============================================================
## 0. Load taxonomy
## ============================================================
taxa <- read.csv("dundreggan_invert_database.csv")

## ============================================================
## 1. Extract species_clean from SILVA taxonomy
## ============================================================
taxa_clean <- taxa %>%
  mutate(
    species_clean = taxonomy %>%
      str_extract("s__[^;]+$") %>%   
      str_replace("^s__", "") %>%    
      str_replace_all("_", " ") %>%  
      str_squish()
  )

## ============================================================
## 2. SAFE GBIF LOOKUP
## ============================================================

get_gbif_lineage <- function(name) {
  if (is.na(name) || name == "" || nchar(name) < 2) {
    return(tibble(
      species_clean = name,
      kingdom = NA, phylum = NA, class = NA, order = NA,
      family = NA, genus = NA, species = NA
    ))
  }
  
  out <- try(classification(name, db="gbif", rows=1)[[1]], silent=TRUE)
  
  if (inherits(out, "try-error") || is.null(out) || !is.data.frame(out)) {
    return(tibble(
      species_clean = name,
      kingdom = NA, phylum = NA, class = NA, order = NA,
      family = NA, genus = NA, species = NA
    ))
  }
  
  get_rank <- function(rank_name) {
    val <- out$name[out$rank == rank_name]
    if (length(val) > 0) val[[1]] else NA
  }
  
  tibble(
    species_clean = name,
    kingdom = get_rank("kingdom"),
    phylum  = get_rank("phylum"),
    class   = get_rank("class"),
    order   = get_rank("order"),
    family  = get_rank("family"),
    genus   = get_rank("genus"),
    species = name
  )
}

## ============================================================
## 3. Run GBIF on unique species names
## ============================================================
species_vec <- taxa_clean$species_clean %>% unique() %>% sort()

lineage_df <- map_df(species_vec, get_gbif_lineage)

## ============================================================
## 4. JOIN WITHOUT RENAMING
## ============================================================
taxa_joined <- taxa_clean %>%
  left_join(lineage_df, by = "species_clean")

## ============================================================
## 5. BUILD FINAL TAXONOMY
## ============================================================
taxa_final <- taxa_joined %>%
  mutate(
    kingdom_final = kingdom,   # SILVA doesn't have kingdom → GBIF only
    phylum_final  = coalesce(phylum.x, phylum.y),
    class_final   = coalesce(class.x, class.y),
    order_final   = coalesce(order.x, order.y),
    family_final  = coalesce(family.x, family.y),
    genus_final   = coalesce(genus.x, genus.y),
    species_final = species_clean
  ) %>%
  select(
    ASV_id,
    kingdom_final, phylum_final, class_final, order_final,
    family_final, genus_final, species_final,
    everything()
  )

write.csv(taxa_final, "dundreggan_invert_taxonomy_refine.csv", row.names = FALSE)
