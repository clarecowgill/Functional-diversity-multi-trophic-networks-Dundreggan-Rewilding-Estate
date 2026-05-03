library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(igraph)
library(vegan)
library(RColorBrewer)
library(patchwork)
library(ggplot2)

# presence-absence long format for each taxa group
verts_pa <- read.csv('acoustic_pa_long.csv')
inverts_pa <- read.csv('invert_pa_long.csv')
fungi_pa <- read.csv('fungi_genus_long.csv')
tree_pa <- read.csv('tree_pa_long.csv')

# traits
verts_traits <- read.csv("verts_traitdata.csv") %>%
  dplyr::mutate(
    scientificName = gsub(" ", "_", scientificName)
  )
inverts_traits <- read.csv('invert_traits_networks.csv')
fungi_traits <- read.csv('fungi_traits_networks.csv')

fungi_traits <- fungi_traits %>%
  dplyr::filter(
    !str_detect(genus, "Incertae_sedis"))
valid_fungal_genera <- fungi_traits %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::pull(genus) %>%
  unique()
fungi_pa <- fungi_pa %>%
  dplyr::filter(genus %in% valid_fungal_genera)

#_____________________________________
# Data prep ----
#_____________________________________
# Regroup fungal traits into the broader lifestyles used in trait diversity analysis
fungi_traits <- fungi_traits %>%
  dplyr::distinct(genus, primary_lifestyle, Decay_substrate_template, 
          Plant_pathogenic_capacity_template, 
          Animal_biotrophic_capacity_template, Fruitbody_type_template) %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::mutate(
    broad_lifestyle = case_when(
      primary_lifestyle == "ectomycorrhizal" ~ "ecm",
      primary_lifestyle == "arbuscular_mycorrhizal" ~ "am",
      
      primary_lifestyle %in% c("root_endophyte", "foliar_endophyte") ~ "endophyte",
      
      primary_lifestyle == "lichenized" ~ "lichenized",
      primary_lifestyle %in% c("algal_parasite", "moss_symbiont") ~ "phototroph_symb",
      
      primary_lifestyle == "soil_saprotroph"   ~ "sapro_soil",
      primary_lifestyle == "litter_saprotroph" ~ "sapro_litter",
      primary_lifestyle == "wood_saprotroph"   ~ "sapro_wood",
      
      primary_lifestyle %in% c(
        "dung_saprotroph",
        "pollen_saprotroph",
        "nectar/tap_saprotroph",
        "unspecified_saprotroph"
      ) ~ "sapro_unspec",
      
      primary_lifestyle == "plant_pathogen" ~ "plant_pathogen",
      primary_lifestyle %in% c(
        "animal_parasite",
        "animal_endosymbiont",
        "mycoparasite",
        "lichen_parasite"
      ) ~ "other_parasites",
      
      TRUE ~ "other"
    )
  ) %>%
  dplyr::mutate(
    broad_lifestyle = factor(broad_lifestyle),
  ) %>%
  dplyr::select(genus, broad_lifestyle, Decay_substrate_template,
                Plant_pathogenic_capacity_template,
                Animal_biotrophic_capacity_template, Fruitbody_type_template)

fungi_traits <- fungi_traits %>%
  dplyr::filter(!is.na(broad_lifestyle))
#_____________________________________
# Create possible edges ----
#_____________________________________
# Fungi links ----
fungi_substrate_basal <- fungi_traits %>%
  dplyr::select(genus, Decay_substrate_template) %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::mutate(
    Decay_substrate_template = str_replace_all(
      Decay_substrate_template,
      "algal_material",
      ""
    )
  ) %>%
  tidyr::separate_rows(Decay_substrate_template, sep = ",") %>%
  dplyr::mutate(substrate = str_trim(Decay_substrate_template)) %>%
  dplyr::filter(substrate != "") %>%
  dplyr::mutate(
    to = case_when(
      substrate == "roots" ~ "PLANTS",
      substrate %in% c(
        "soil",
        "leaf/fruit/seed",
        "wood",
        "dung",
        "animal_material",
        "fungal_material"
      ) ~ "SOIL_ORGANIC_MATTER",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(to)) %>%
  dplyr::select(genus, to) %>%
  distinct()
fungi_no_substrate <- fungi_traits %>%
  dplyr::filter(
    is.na(Decay_substrate_template) | Decay_substrate_template == ""
  )
fungi_lifestyle_basal <- fungi_no_substrate %>%
  dplyr::mutate(
    to = case_when(
      broad_lifestyle %in% c("ecm", "am", "endophyte", "plant_pathogen") ~ "PLANTS",
      broad_lifestyle %in% c(
        "sapro_soil",
        "sapro_litter",
        "sapro_wood",
        "sapro_unspec"
      ) ~ "SOIL_ORGANIC_MATTER",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(to)) %>%
  dplyr::select(genus, to) %>%
  distinct()
edges_fungi_basal <- bind_rows(
  fungi_substrate_basal,
  fungi_lifestyle_basal
) %>%
  distinct() %>%
  dplyr::transmute(
    from = genus,
    to   = to,
    interaction = "basal_feeding",
    rule = "fungi_substrate_or_lifestyle"
  )

# Invertebrate consumers ----
## Invert Fungivory ----
# align invertebrate microhabitats with decay substrates for potential fungivory interactions
# 1) Microhabitat -> allowed fungal substrates
invert_microhabitat_substrate <- tribble(
  ~microhabitat, ~substrate,
  "ep", "leaf/fruit/seed",
  "ep", "wood",
  "ep", "dung",
  "ep", "animal_material",
  "ep", "fungal_material",
  "ep", "soil",
  
  "he", "leaf/fruit/seed",
  "he", "soil",
  "he", "roots",
  "he", "wood",
  "he", "dung",
  "he", "animal_material",
  "he", "fungal_material",
  
  "eu", "soil",
  "eu", "roots",
  "eu", "fungal_material"
)

# 2) Parse fungal substrates (drop algal_material, blanks)
fungal_substrate_long <- fungi_traits %>%
  dplyr::select(genus, Decay_substrate_template) %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::mutate(Decay_substrate_template = str_replace_all(Decay_substrate_template,
                                                           "algal_material", "")) %>%
  separate_rows(Decay_substrate_template, sep = ",") %>%
  dplyr::mutate(substrate = str_trim(Decay_substrate_template)) %>%
  dplyr::filter(substrate != "") %>%
  dplyr::distinct(genus, substrate)

# regroup substrates (pollen, sugary to plant material and drop protist material)
fungal_substrate_long <- fungal_substrate_long %>%
  mutate(
    substrate = case_when(
      substrate %in% c("sugar-rich_substrates", "pollen") ~ "leaf/fruit/seed",
      substrate == "protist_material" ~ NA_character_,
      TRUE ~ substrate
    )
  ) %>%
  filter(!is.na(substrate)) %>%
  distinct(genus, substrate)

# 3) Identify fungivores (strict: feeding_guild == "f")
fungivores <- inverts_traits %>%
  dplyr::rename(taxon_id = X) %>%
  dplyr::filter(
    feeding_guild == "f" |
      (
        feeding_guild == "o" &
          !is.na(secondary_guilds) &
          str_detect(
            paste0(",", secondary_guilds, ","),
            ",f,"
          )
      )
  ) %>%
  dplyr::filter(!is.na(microhabitat)) %>%
  dplyr::select(taxon_id, feeding_guild, secondary_guilds, microhabitat) %>%
  distinct()

edges_fungivory <- fungivores %>%
  inner_join(
    invert_microhabitat_substrate,
    by = "microhabitat",
    relationship = "many-to-many"
  ) %>%
  inner_join(
    fungal_substrate_long,
    by = "substrate",
    relationship = "many-to-many"
  ) %>%
  transmute(
    from = taxon_id,
    to   = genus,
    interaction = "fungivory",
    rule = paste0(
      "guild=", feeding_guild,
      ifelse(!is.na(secondary_guilds),
             paste0("(", secondary_guilds, ")"),
             ""),
      "; microhabitat=", microhabitat,
      "; substrate=", substrate
    )
  ) %>%
  distinct()

## Inverts --> basal nodes ----
basal_nodes <- tibble(
  node_id = c(
    "SOIL_ORGANIC_MATTER",
    "BACTERIA",
    "PLANTS"),
  node_type = "basal"
)

edges_invert_basal_core <- inverts_traits %>%
  dplyr::rename(taxon_id = X) %>% 
  dplyr::filter(feeding_guild %in% c("d", "b", "h")) %>%
  dplyr::mutate(
    to = case_when(
      feeding_guild == "d" ~ "SOIL_ORGANIC_MATTER",
      feeding_guild == "b" ~ "BACTERIA",
      feeding_guild == "h" ~ "PLANTS"
    )
  ) %>%
  dplyr::transmute(
    from = taxon_id,
    to   = to,
    interaction = "basal_feeding",
    rule = paste0("feeding_guild=", feeding_guild)
  ) %>%
  distinct()
edges_invert_basal_omni <- inverts_traits %>%
  dplyr::rename(taxon_id = X) %>%
  dplyr::filter(
    feeding_guild == "o",
    !is.na(secondary_guilds)
  ) %>%
  tidyr::separate_rows(secondary_guilds, sep = ",") %>%
  dplyr::mutate(
    secondary_guilds = str_trim(secondary_guilds),
    to = case_when(
      secondary_guilds == "d" ~ "SOIL_ORGANIC_MATTER",
      secondary_guilds == "b" ~ "BACTERIA",
      secondary_guilds == "h" ~ "PLANTS",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(to)) %>%
  dplyr::transmute(
    from = taxon_id,
    to   = to,
    interaction = "basal_feeding",
    rule = paste0("omnivore_secondary=", secondary_guilds)
  ) %>%
  distinct()
edges_invert_basal <- bind_rows(
  edges_invert_basal_core,
  edges_invert_basal_omni
) %>%
  distinct()

## Invertebrate predators ----
size_hierarchy <- tribble(
  ~pred_size, ~prey_size,
  "micro", "micro",
  "meso",  "micro",
  "meso",  "meso",
  "macro", "micro",
  "macro", "meso",
  "macro", "macro"
)
# 1) define microhabitat accessibility
pred_prey_micro <- tibble::tribble(
  ~pred_micro, ~prey_micro,
  "eu", "eu",
  "eu", "he",
  "he", "he",
  "he", "eu",
  "he", "ep",
  "ep", "ep",
  "ep", "he"
)

predators <- inverts_traits %>%
  rename(taxon_id = X) %>%
  filter(feeding_guild == "c" | (feeding_guild == "o" & str_detect(paste0(",",secondary_guilds,","), ",c,"))) %>%
  filter(!is.na(size_class), !is.na(microhabitat)) %>%
  transmute(predator_id = taxon_id,
            predator_size = size_class,
            predator_micro = microhabitat)
prey <- inverts_traits %>%
  rename(taxon_id = X) %>%
  filter(!is.na(size_class), !is.na(microhabitat)) %>%
  transmute(prey_id = taxon_id,
            prey_size = size_class,
            prey_micro = microhabitat)
edges_predation <- predators %>%
  inner_join(size_hierarchy, by = c("predator_size" = "pred_size")) %>%
  inner_join(prey, by = c("prey_size" = "prey_size")) %>%
  inner_join(pred_prey_micro, by = c("predator_micro" = "pred_micro",
                                    "prey_micro" = "prey_micro")) %>%
  filter(predator_id != prey_id) %>%
  transmute(from = predator_id,
            to   = prey_id,
            interaction = "invertivory",
            rule = paste0("pred_size=", predator_size,
                          "; prey_size=", prey_size,
                          "; pred_micro=", predator_micro,
                          "; prey_micro=", prey_micro)) %>%
  distinct()

# Vertebrate consumers ----
verts_diet_long <- verts_traits %>%
  dplyr::select(scientificName, diet_component) %>%
  dplyr::filter(!is.na(diet_component)) %>%
  tidyr::separate_rows(diet_component, sep = ",") %>%
  dplyr::mutate(
    diet_component = stringr::str_trim(diet_component)
  ) %>%
  distinct()
ground_verts <- verts_traits %>%
  dplyr::select(scientificName, ForStrat.ground) %>%
  dplyr::filter(!is.na(ForStrat.ground)) %>%
  dplyr::filter(ForStrat.ground > 0) %>%
  dplyr::select(scientificName) %>%
  distinct()

## Vertebrates --> basal nodes
edges_vert_plants <- verts_diet_long %>%
  dplyr::filter(diet_component %in% c("plant", "seeds")) %>%
  dplyr::transmute(
    from = scientificName,
    to   = "PLANTS",
    interaction = "herbivory",
    rule = paste0("diet=", diet_component)
  ) %>%
  distinct()
# edges_vert_non_soil_inverts <- verts_diet_long %>%
#   dplyr::filter(diet_component %in% c(
#     "foliage_inverts",
#     "aerial_inverts",
#     "aquatic_inverts"
#   )) %>%
#   dplyr::transmute(
#     from = scientificName,
#     to   = "NON_SOIL_INVERTEBRATES",
#     interaction = "invertivory",
#     rule = paste0("diet=", diet_component)
#   ) %>%
#   distinct()

## Vertebrate invertivores ----
invert_prey_pool <- inverts_traits %>%
  dplyr::rename(taxon_id = X) %>%
  dplyr::filter(size_class %in% c("meso", "macro")) %>%
  dplyr::select(taxon_id) %>%
  distinct()
edges_vert_soil_inverts <- verts_diet_long %>%
  dplyr::filter(diet_component == "soil_inverts") %>%
  dplyr::inner_join(ground_verts, by = "scientificName") %>%
  dplyr::select(scientificName) %>%
  distinct() %>%
  dplyr::mutate(key = 1) %>%
  dplyr::inner_join(
    invert_prey_pool %>% dplyr::mutate(key = 1),
    by = "key",
    relationship = "many-to-many"
  ) %>%
  dplyr::transmute(
    from = scientificName,
    to   = taxon_id,
    interaction = "invertivory",
    rule = "diet=soil_inverts; ground_use>0; prey_size=meso+macro"
  ) %>%
  distinct()
# edges_vert_nonsoil_inverts <- verts_diet_long %>%
#   dplyr::filter(
#     diet_component %in% c(
#       "foliage_inverts",
#       "aquatic_inverts",
#       "aerial_inverts"
#     )
#   ) %>%
#   dplyr::transmute(
#     from = scientificName,
#     to   = "NON_SOIL_INVERTEBRATES",
#     interaction = "invertivory",
#     rule = paste0("diet=", diet_component)
#   ) %>%
#   distinct()

## Vertebrate fungivores ----
vertebrate_fungivory_fruitbodies <- c(
  "agaricoid",
  "clavarioid",
  "phalloid",
  "tremelloid",
  "cyphelloid",
  "polyporoid",
  "corticioid",
  "gasteroid",
  "gasteroid-hypogeous",
  "apothecium_(hymenium_on_surface)"
)
fungi_prey_pool_vert <- fungi_traits %>%
  dplyr::filter(
    Fruitbody_type_template %in% vertebrate_fungivory_fruitbodies
  ) %>%
  dplyr::select(genus) %>%
  distinct()
edges_vert_fungi <- verts_diet_long %>%
  dplyr::filter(diet_component == "fungi") %>%
  dplyr::inner_join(ground_verts, by = "scientificName") %>%
  dplyr::select(scientificName) %>%
  distinct() %>%
  dplyr::mutate(key = 1) %>%
  dplyr::inner_join(
    fungi_prey_pool_vert %>% mutate(key = 1),
    by = "key",
    relationship = "many-to-many"
  ) %>%
  dplyr::transmute(
    from = scientificName,
    to   = genus,
    interaction = "fungivory",
    rule = "diet=fungi; ground_use>0; filter=fruitbody_accessible"
  ) %>%
  distinct()

## Vertebrate vertivores ----
strata_cols <- c(
  "ForStrat.watbelowsurf",
  "ForStrat.wataroundsurf",
  "ForStrat.ground",
  "ForStrat.understory",
  "ForStrat.midhigh",
  "ForStrat.canopy",
  "ForStrat.aerial"
)
vert_strata_long <- verts_traits %>%
  dplyr::select(scientificName, all_of(strata_cols)) %>%
  tidyr::pivot_longer(
    cols = all_of(strata_cols),
    names_to = "stratum",
    values_to = "use"
  ) %>%
  dplyr::filter(!is.na(use), use > 0) %>%
  dplyr::select(scientificName, stratum) %>%
  distinct()
vert_predators <- verts_diet_long %>%
  dplyr::filter(diet_component == "vertebrates") %>%
  dplyr::distinct(scientificName) %>%
  dplyr::rename(predator = scientificName)
edges_vert_verts <- vert_predators %>%
  dplyr::inner_join(
    vert_strata_long,
    by = c("predator" = "scientificName")
  ) %>%
  dplyr::inner_join(
    vert_strata_long,
    by = "stratum",
    relationship = "many-to-many"
  ) %>%
  dplyr::rename(prey = scientificName) %>%
  dplyr::filter(predator != prey) %>%   # no self-links
  dplyr::transmute(
    from = predator,
    to   = prey,
    interaction = "vertivory",
    rule = paste0("diet=vertebrates; shared_stratum=", stratum)
  ) %>%
  distinct()

# All possible links
edges_metaweb <- bind_rows(
  # edges_fungi_basal,
  # edges_fungivory,
  # edges_invert_basal,
  # edges_predation,
  edges_vert_plants,
  # edges_vert_soil_inverts,
  # edges_vert_fungi,
  edges_vert_verts
)
edges_metaweb %>%
  count(interaction)

## Network Metrics ----
# Build igraph for a given site (two versions)
make_site_graphs <- function(site_id, edges_site, nodes_site, basal_nodes) {
  
  nodes_s <- nodes_site %>%
    filter(site == site_id) %>%
    distinct(node)
  
  edges_s <- edges_site %>%
    filter(site == site_id) %>%
    semi_join(nodes_s, by = c("from" = "node")) %>%
    semi_join(nodes_s, by = c("to"   = "node")) %>%
    distinct(from, to, interaction, rule)
  
  # keep only nodes that are actually in at least one edge
  nodes_in_edges <- tibble(node = unique(c(edges_s$from, edges_s$to)))
  
  nodes_s <- nodes_s %>%
    inner_join(nodes_in_edges, by = "node") %>%
    mutate(
      is_basal = node %in% basal_nodes$node_id,
      type = case_when(
        node %in% verts_pa$species    ~ "verts",
        node %in% inverts_pa$taxon_id ~ "inverts",
        node %in% fungi_pa$genus      ~ "fungi",
        node %in% basal_nodes$node_id ~ "basal",
        TRUE ~ NA_character_
      )
    )
  
  # consumer -> resource (your current orientation)
  g_cr <- graph_from_data_frame(
    edges_s %>% dplyr::select(from, to),
    directed = TRUE,
    vertices = nodes_s %>% transmute(name = node, is_basal)
  )
  
  # resource -> consumer (flip) for trophic metrics / coherence / robustness
  edges_rc <- edges_s %>% transmute(from = to, to = from)
  g_rc <- graph_from_data_frame(
    edges_rc,
    directed = TRUE,
    vertices = nodes_s %>% transmute(name = node, is_basal, type)
  )
  
  list(g_cr = g_cr, g_rc = g_rc, edges_s = edges_s, nodes_s = nodes_s)
}

# Robustness R50 by random primary removals (excluding basal),
# cascading extinctions remove any non-basal consumer with in-degree 0.
robustness_R50 <- function(g_rc, basal_nodes_vec, n_sims = 200) {
  
  all_nodes <- V(g_rc)$name
  basal <- intersect(all_nodes, basal_nodes_vec)
  candidates <- setdiff(all_nodes, basal)
  if (length(candidates) < 2) return(NA_real_)
  prop_primary_to_50 <- numeric(n_sims)
  for (sim in seq_len(n_sims)) {
    g <- g_rc
    removed_primary <- character(0)
    # random removal order
    order <- sample(candidates)
    # track extinctions as nodes removed (primary + secondary)
    extinct <- character(0)
    for (v in order) {
      if (!(v %in% V(g)$name)) next
      
      # primary removal
      g <- delete_vertices(g, v)
      removed_primary <- c(removed_primary, v)
      
      repeat {
        non_basal_now <- setdiff(V(g)$name, basal)
        if (length(non_basal_now) == 0) break
        
        # secondary extinctions: consumers without any incoming resources
        indeg <- degree(g, mode = "in")
        to_drop <- intersect(names(indeg)[indeg == 0], non_basal_now)
        
        if (length(to_drop) == 0) break
        g <- delete_vertices(g, to_drop)
      }
      
      # evaluate loss threshold
      extinct_now <- setdiff(all_nodes, V(g)$name)
      if (length(extinct_now) / length(all_nodes) >= 0.5) {
        prop_primary_to_50[sim] <- length(removed_primary) / length(candidates)
        break
      }
    }
    
    if (prop_primary_to_50[sim] == 0) prop_primary_to_50[sim] <- 1
  }
  
  mean(prop_primary_to_50, na.rm = TRUE)
}

# Modularity: Blackman uses igraph multilevel community (Louvain) (undirected). :contentReference[oaicite:6]{index=6}
modularity_louvain <- function(g_rc) {
  if (vcount(g_rc) < 3 || ecount(g_rc) < 2) return(NA_real_)
  gu <- as_undirected(g_rc, mode = "collapse")

  modularity(cluster_louvain(gu))
}

calc_site_metrics <- function(site_id, edges_site, nodes_site, basal_nodes, n_sims = 200) {
  gs <- make_site_graphs(site_id, edges_site, nodes_site, basal_nodes)
  g_rc <- gs$g_rc
  if (vcount(g_rc) < 2 || ecount(g_rc) < 1) {
    return(tibble(
      site = site_id,
      S = vcount(g_rc),
      L = ecount(g_rc),
      link_density = NA_real_,
      connectance = NA_real_,
      modularity = NA_real_,
      robustness_R50 = NA_real_
    ))
  }
  S <- vcount(g_rc)
  L <- ecount(g_rc)
  if (S == 0) {
    return(tibble(site = site_id, S = 0, L = 0))
  }
  link_density <- L / S
  connectance  <- if (S > 1) L / (S^2) else NA_real_
  s <- trophic_levels_prey_avg(g_rc)
  mod  <- modularity_louvain(g_rc)
  R50 <- robustness_R50(g_rc, basal_nodes$node_id, n_sims = n_sims)
  tibble(
    site = site_id,
    S = S,
    L = L,
    link_density = link_density,
    connectance  = connectance,
    modularity = mod,
    robustness_R50 = R50
  )
}
# 
# # run across all sites
# all_sites <- sort(unique(edges_site$site))

# metrics_df <- bind_rows(lapply(all_sites, calc_site_metrics,
#                                edges_site = edges_site,
#                                nodes_site = nodes_site,
#                                basal_nodes = basal_nodes,
#                                n_sims = 200))
# metrics_df

# Sample-level networks ----
# Each soil sample has its own network, but all are overlayed with the vertebrate interactions for the site

# soil sample nodes
# samples retained after QC in each dataset
fungi_samples <- fungi_pa %>%
  filter(presence > 0) %>%
  distinct(sample) %>%
  pull(sample)

invert_samples <- inverts_pa %>%
  filter(presence > 0) %>%
  distinct(sample) %>%
  pull(sample)

# intersection: samples with BOTH fungi and inverts
soil_samples_valid <- intersect(fungi_samples, invert_samples)
fungi_pa_net <- fungi_pa %>%
  filter(sample %in% soil_samples_valid)

inverts_pa_net <- inverts_pa %>%
  filter(sample %in% soil_samples_valid)

nodes_soil <- bind_rows(
  inverts_pa_net %>% filter(presence > 0) %>% transmute(sample, site, node = taxon_id),
  fungi_pa_net   %>% filter(presence > 0) %>% transmute(sample, site, node = genus)
) %>%
  distinct(sample, site, node)

# site-level vertebrate nodes -BOTH TREE AND ACOUSTIC
nodes_verts_tree <- tree_pa %>%
  distinct(sample, species) %>%
  mutate(
    site = as.integer(stringr::str_extract(sample, "(?<=D_)\\d{2}"))
  ) %>%
  transmute(site, node = species) %>%
  filter(!is.na(site)) %>%
  distinct()

nodes_verts_site_acoustic <- verts_pa %>%
  filter(presence > 0) %>%
  distinct(site, species) %>%
  transmute(site, node = species)

nodes_verts_site_union <- bind_rows(
  nodes_verts_site_acoustic,
  nodes_verts_tree
) %>%
  distinct(site, node)

# combine into network nodes
# # run with acoustics only
# nodes_soil_hybrid <- nodes_soil %>%
#   left_join(nodes_verts_site_acoustic, by = "site") %>%
#   rename(vert_node = node.y) %>%
#   rename(node = node.x) %>%
#   group_by(sample, site) %>%
#   summarise(
#     nodes = list(unique(c(node, vert_node, basal_nodes$node_id))),
#     .groups = "drop"
#   ) %>%
#   unnest(nodes) %>%
#   rename(node = nodes) %>%
#   filter(!is.na(node)) %>%
#   distinct(sample, site, node)

# run with tree roller samples included
nodes_soil_hybrid <- nodes_soil %>%
  left_join(nodes_verts_site_union, by = "site") %>%
  rename(vert_node = node.y) %>%
  rename(node = node.x) %>%
  group_by(sample, site) %>%
  summarise(
    nodes = list(unique(c(node, vert_node, basal_nodes$node_id))),
    .groups = "drop"
  ) %>%
  unnest(nodes) %>%
  rename(node = nodes) %>%
  filter(!is.na(node)) %>%
  distinct(sample, site, node)

# build edges
make_hybrid_edges <- function(sample_id, nodes_df, edges_metaweb) {
  
  nodes_s <- nodes_df %>%
    filter(sample == sample_id) %>%
    distinct(node)
  
  edges_metaweb %>%
    semi_join(nodes_s, by = c("from" = "node")) %>%
    semi_join(nodes_s, by = c("to"   = "node")) %>%
    distinct(from, to, interaction, rule) %>%
    mutate(sample = sample_id)
}
soil_samples <- sort(unique(nodes_soil_hybrid$sample))

edges_soil_hybrid <- bind_rows(
  lapply(
    soil_samples,
    make_hybrid_edges,
    nodes_df = nodes_soil_hybrid,
    edges_metaweb = edges_metaweb
  )
)

# compute metrics per soil sample (hybrid networks)
calc_soil_sample_metrics <- function(sample_id,
                                     edges_df,
                                     nodes_df,
                                     basal_nodes,
                                     n_sims = 200) {
  
  nodes_s <- nodes_df %>%
    filter(sample == sample_id) %>%
    distinct(node) %>%
    mutate(is_basal = node %in% basal_nodes$node_id)
  
  edges_s <- edges_df %>%
    filter(sample == sample_id) %>%
    distinct(from, to)
  
  nodes_s <- nodes_s %>%
    filter(node %in% c(edges_s$from, edges_s$to))
  
  if (nrow(edges_s) < 1 || nrow(nodes_s) < 2) {
    return(tibble(
      sample = sample_id,
      S = nrow(nodes_s),
      L = nrow(edges_s),
      link_density = NA_real_,
      connectance  = NA_real_,
      modularity   = NA_real_,
      robustness_R50 = NA_real_
    ))
  }
  
  # resource -> consumer orientation
  g_rc <- igraph::graph_from_data_frame(
    edges_s %>% transmute(from = to, to = from),
    directed = TRUE,
    vertices = nodes_s %>% transmute(name = node, is_basal)
  )
  
  S <- vcount(g_rc)
  L <- ecount(g_rc)
  
  tibble(
    sample = sample_id,
    S = S,
    L = L,
    link_density = L / S,
    connectance  = if (S > 1) L / (S^2) else NA_real_,
    modularity   = modularity_louvain(g_rc),
    robustness_R50 = robustness_R50(g_rc, basal_nodes$node_id, n_sims)
  )
}

metrics_soil_hybrid <- bind_rows(
  lapply(
    soil_samples,
    calc_soil_sample_metrics,
    edges_df = edges_soil_hybrid,
    nodes_df = nodes_soil_hybrid,
    basal_nodes = basal_nodes,
    n_sims = 200
  )
)

# plotting function
plot_soil_hybrid_network <- function(sample_id) {
  
  ## ---- edges for this sample ----
  edges_s <- edges_soil_hybrid %>%
    dplyr::filter(sample == sample_id)
  
  ## ---- infer realised feeding roles from edges ----
  predation_roles <- edges_s %>%
    dplyr::mutate(
      prey_type = dplyr::case_when(
        to %in% verts_traits$scientificName ~ "vertebrate",
        to %in% inverts_pa$taxon_id          ~ "invertebrate",
        to %in% fungi_pa$genus               ~ "fungi",
        to %in% basal_nodes$node_id          ~ "basal",
        TRUE                                 ~ "other"
      )
    ) %>%
    dplyr::distinct(from, prey_type)
  
  node_roles <- predation_roles %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(
      eats_vertebrates = any(prey_type == "vertebrate"),
      eats_inverts     = any(prey_type == "invertebrate"),
      eats_fungi       = any(prey_type == "fungi"),
      .groups = "drop"
    )
  
  ## ---- nodes for this sample ----
  nodes_s <- nodes_soil_hybrid %>%
    dplyr::filter(sample == sample_id) %>%
    dplyr::left_join(node_roles, by = c("node" = "from")) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        node %in% basal_nodes$node_id          ~ "basal",
        node %in% fungi_pa$genus               ~ "fungi",
        node %in% inverts_pa$taxon_id          ~ "inverts",
        node %in% verts_traits$scientificName ~ "verts",
        TRUE                                   ~ NA_character_
      ),
      
      ## ---- trophic vertical positioning ----
      y = dplyr::case_when(
        type == "basal" ~ 1.0,
        type == "fungi" ~ 2.0,
        
        # invertebrates
        type == "inverts" & eats_vertebrates ~ 3.6,
        type == "inverts" & eats_inverts     ~ 3.3,
        type == "inverts"                    ~ 3.0,
        
        # vertebrates
        type == "verts" & eats_vertebrates   ~ 4.6,
        type == "verts" & eats_inverts       ~ 4.2,
        type == "verts"                      ~ 3.9,
        
        TRUE ~ NA_real_
      )
    )
  
  ## ---- gentle adaptive horizontal spacing ----
  nodes_s <- nodes_s %>%
    dplyr::group_by(type, y) %>%
    dplyr::mutate(n_row = dplyr::n()) %>%
    dplyr::ungroup()
  
  max_row <- max(nodes_s$n_row, na.rm = TRUE)
  
  min_width <- 0.45  # 👈 adjust if needed (0.35–0.55 sweet spot)
  
  nodes_s <- nodes_s %>%
    dplyr::group_by(type, y) %>%
    dplyr::mutate(
      row_width =
        min_width +
        (1 - min_width) * (n_row / max_row),
      
      x = seq(
        from = -dplyr::first(row_width),
        to   =  dplyr::first(row_width),
        length.out = dplyr::first(n_row)
      )
    ) %>%
    dplyr::ungroup()
  nodes_s <- nodes_s %>%
    dplyr::filter(node %in% c(edges_s$from, edges_s$to))
  ## ---- build igraph ----
  g <- igraph::graph_from_data_frame(
    edges_s %>% dplyr::select(from, to),
    directed = TRUE,
    vertices = nodes_s %>%
      dplyr::select(name = node, x, y, type)
  )
  
  ## ---- plot ----
  p <- ggraph::ggraph(g, layout = "manual", x = x, y = y) +
    ggraph::geom_edge_link(
      alpha = 0.25,
      colour = "grey50",
      width = 0.3
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(shape = type),
      size = 2.4
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        basal   = 15,
        fungi   = 16,
        inverts = 17,
        verts   = 18
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.8, 4.8),
      breaks = c(1, 2, 3.0, 3.3, 4.2, 4.6),
      labels = c(
        "Basal resources",
        "Fungi",
        "Invertebrates",
        "Predatory invertebrates",
        "Invertivorous vertebrates",
        "Apex vertebrates"
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste("Soil sample", sample_id))
  
  return(p)
}

p_05_12 <- plot_soil_hybrid_network("05_12")

#ggsave(file.path("Plots", "07_10_network.png"), plot = p_07_10, width = 6, height = 5, dpi = 300)

## Boxplots of metrics ----
ylgn <- brewer.pal(9, "YlGn")
age_cols_9 <- rev(colorRampPalette(ylgn[2:9])(9))
age_cols <- c("#002e23", age_cols_9)  # site 1 darkest, site 10 lightest
age_cols_named <- setNames(age_cols, as.character(1:10))

# prep metrics dataframe
metrics_plot <- metrics_soil_hybrid %>%
  mutate(
    site = as.integer(substr(sample, 1, 2)),
    woodland_age = recode(site,
                          `1`  = 250, `2`  = 250,
                          `3`  = 100,
                          `4`  = 50,
                          `5`  = 22,
                          `6`  = 14,
                          `7`  = 11,
                          `8`  = 8,
                          `9`  = 6,
                          `10` = 4
    ),
    
    box_group = case_when(
      woodland_age == 250 & site == 2 ~ "left",
      woodland_age == 250 & site == 1 ~ "right",
      TRUE                            ~ "all"
    ),
    
    box_group = factor(box_group, levels = c("left", "right", "all")),
    site = factor(site, levels = rev(1:10))
  )

plot_metric_by_age <- function(df, metric, ylab = metric) {
  
  ggplot(
    df,
    aes(
      x = factor(woodland_age),
      y = .data[[metric]],
      fill = site,
      group = interaction(woodland_age, box_group)
    )
  ) +
    geom_boxplot(
      alpha = 0.8,
      outlier.shape = NA,
      position = position_dodge2(width = 0.9, preserve = "single")
    ) +
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width  = 0.9
      ),
      alpha = 0.4,
      size = 1,
      colour = "black"
    )+
    scale_fill_manual(values = age_cols_named) +
    labs(
      x = "Woodland age (years)",
      y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
}

p_S <- plot_metric_by_age(metrics_plot, "S", "Number of nodes (S)")
p_L <- plot_metric_by_age(metrics_plot, "L", "Number of links (L)")
p_ld <- plot_metric_by_age(metrics_plot, "link_density", "Linkage density (L/S)")
p_c  <- plot_metric_by_age(metrics_plot, "connectance", "Connectance (L/S²)")
p_m  <- plot_metric_by_age(metrics_plot, "modularity", "Modularity")
p_r  <- plot_metric_by_age(metrics_plot, "robustness_R50", "Robustness (R50)")

site_metrics <- (p_S | p_ld) /
  (p_c | p_r)
site_metrics

ggsave(file.path("Plots", "network_metrics_verts.png"), plot = site_metrics, width = 8, height = 6, dpi = 300)



# Study-wide indicative figure ----

functional_nodes <- tibble::tribble(
  ~group, ~domain, ~trophic_level,
  "Plants","basal",1,"Soil organic matter","basal",1,"Bacteria","basal",1,
  "Mycorrhizal fungi","fungi",2,"Saprotrophic fungi","fungi",2,"Endophytic fungi","fungi",2,"Plant pathogenic fungi","fungi",2,"Parasitic fungi","fungi",2,
  "Herbivorous invertebrates","inverts",3,"Detritivorous invertebrates","inverts",3,"Bacterivorous invertebrates","inverts",3,
  "Fungivorous invertebrates","inverts",3,"Omnivorous invertebrates","inverts",3,"Predatory invertebrates","inverts",3.3,
  "Herbivorous vertebrates","verts",4,"Omnivorous vertebrates","verts",4.2,"Invertivorous vertebrates","verts",4.4,"Vertebrate predators","verts",4.7
)

fungi_groups <- fungi_traits %>%
  mutate(group = case_when(
    broad_lifestyle %in% c("ecm","am") ~ "Mycorrhizal fungi",
    broad_lifestyle %in% c("sapro_soil","sapro_litter","sapro_wood","sapro_unspec") ~ "Saprotrophic fungi",
    broad_lifestyle == "endophyte" ~ "Endophytic fungi",
    broad_lifestyle == "plant_pathogen" ~ "Plant pathogenic fungi",
    broad_lifestyle == "other_parasites" ~ "Parasitic fungi",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(group)) %>% distinct(genus, group)

invert_groups <- inverts_traits %>%
  rename(taxon_id = X) %>%
  mutate(group = case_when(
    feeding_guild == "h" ~ "Herbivorous invertebrates",
    feeding_guild == "d" ~ "Detritivorous invertebrates",
    feeding_guild == "b" ~ "Bacterivorous invertebrates",
    feeding_guild == "f" ~ "Fungivorous invertebrates",
    feeding_guild == "o" ~ "Omnivorous invertebrates",
    feeding_guild == "c" ~ "Predatory invertebrates",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(group)) %>% distinct(taxon_id, group)

vert_groups <- verts_traits %>%
  select(scientificName, diet_component, tl_test) %>%
  separate_rows(diet_component, sep = ",") %>%
  mutate(diet_component = str_trim(diet_component)) %>%
  group_by(scientificName) %>%
  summarise(
    eats_plants = any(tl_test %in% c("F","H", "G")),
    eats_inverts = any(diet_component == "soil_inverts"),
    eats_verts = any(tl_test == "V"),
    eats_all = any(tl_test == "O"),
    .groups = "drop"
  ) %>%
  mutate(group = case_when(
    eats_verts ~ "Vertebrate predators",
    eats_all ~ "Omnivorous vertebrates",
    eats_plants ~ "Herbivorous vertebrates",
    eats_inverts ~ "Invertivorous vertebrates",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(group)) %>% distinct(scientificName, group)

node_sizes <- bind_rows(
  fungi_groups %>% transmute(taxon = genus, group),
  invert_groups %>% transmute(taxon = taxon_id, group),
  vert_groups %>% transmute(taxon = scientificName, group)
) %>% count(group, name = "n_taxa")

nodes_fig1 <- functional_nodes %>%
  left_join(node_sizes, by = "group") %>%
  mutate(n_taxa = replace_na(n_taxa, 0))

metaweb_edges <- tibble::tribble(
  ~from,~to,~interaction,
  "Mycorrhizal fungi","Plants","basal_plant","Endophytic fungi","Plants","basal_plant",
  "Plant pathogenic fungi","Plants","basal_plant","Parasitic fungi","Plants","basal_plant",
  "Saprotrophic fungi","Soil organic matter","detritivory",
  "Herbivorous invertebrates","Plants","herbivory",
  "Detritivorous invertebrates","Soil organic matter","detritivory",
  "Bacterivorous invertebrates","Bacteria","bacterivory",
  "Fungivorous invertebrates","Saprotrophic fungi","fungivory",
  "Fungivorous invertebrates","Mycorrhizal fungi","fungivory",
  "Fungivorous invertebrates","Plant pathogenic fungi","fungivory",
  "Fungivorous invertebrates","Parasitic fungi","fungivory",
  "Fungivorous invertebrates","Endophytic fungi","fungivory",
  "Predatory invertebrates","Herbivorous invertebrates","invertivory",
  "Predatory invertebrates","Detritivorous invertebrates","invertivory",
  "Predatory invertebrates","Fungivorous invertebrates","invertivory",
  "Predatory invertebrates","Bacterivorous invertebrates","invertivory",
  "Predatory invertebrates","Omnivorous invertebrates","invertivory",
  "Herbivorous vertebrates","Plants","herbivory",
  "Invertivorous vertebrates","Predatory invertebrates","invertivory",
  "Invertivorous vertebrates","Detritivorous invertebrates","invertivory",
  "Invertivorous vertebrates","Fungivorous invertebrates","invertivory",
  "Invertivorous vertebrates","Omnivorous invertebrates","invertivory",
  "Invertivorous vertebrates","Herbivorous invertebrates","invertivory",
  "Omnivorous invertebrates","Plants","herbivory",
  "Omnivorous invertebrates","Saprotrophic fungi","fungivory",
  "Omnivorous invertebrates","Mycorrhizal fungi","fungivory",
  "Omnivorous invertebrates","Plant pathogenic fungi","fungivory",
  "Omnivorous invertebrates","Parasitic fungi","fungivory",
  "Omnivorous invertebrates","Endophytic fungi","fungivory",
  "Omnivorous invertebrates","Herbivorous invertebrates","invertivory",
  "Omnivorous invertebrates","Predatory invertebrates","invertivory",
  "Omnivorous invertebrates","Detritivorous invertebrates","invertivory",
  "Omnivorous invertebrates","Bacterivorous invertebrates","invertivory",
  "Omnivorous invertebrates","Fungivorous invertebrates","invertivory",
  "Omnivorous vertebrates","Plants","herbivory",
  "Omnivorous vertebrates","Saprotrophic fungi","fungivory",
  "Omnivorous vertebrates","Mycorrhizal fungi","fungivory",
  "Omnivorous vertebrates","Plant pathogenic fungi","fungivory",
  "Omnivorous vertebrates","Parasitic fungi","fungivory",
  "Omnivorous vertebrates","Endophytic fungi","fungivory",
  "Omnivorous vertebrates","Predatory invertebrates","invertivory",
  "Omnivorous vertebrates","Vertebrate predators","vertivory",
  "Omnivorous vertebrates","Herbivorous vertebrates","vertivory",
  "Vertebrate predators","Invertivorous vertebrates","vertivory",
  "Vertebrate predators","Herbivorous vertebrates","vertivory",
  "Vertebrate predators","Omnivorous vertebrates","vertivory"
) %>% distinct()

trophic_order <- c(
  "Plants","Soil organic matter","Bacteria",
  "Saprotrophic fungi","Mycorrhizal fungi","Endophytic fungi","Plant pathogenic fungi","Parasitic fungi",
  "Detritivorous invertebrates","Bacterivorous invertebrates","Herbivorous invertebrates","Fungivorous invertebrates","Omnivorous invertebrates","Predatory invertebrates",
  "Herbivorous vertebrates","Omnivorous vertebrates","Invertivorous vertebrates","Vertebrate predators"
)

nodes_fig1 <- nodes_fig1 %>%
  mutate(group_chr = group, group_ord = factor(group, levels = trophic_order)) %>%
  arrange(group_ord) %>%
  mutate(angle = (seq_len(n())-1)*(2*pi/n())+(pi/2-0.15), x = -cos(angle), y = sin(angle)) %>%
  transmute(group = group_chr, domain, n_taxa, x, y)

metaweb_edges <- metaweb_edges %>%
  filter(from %in% nodes_fig1$group, to %in% nodes_fig1$group)

g_fig1 <- graph_from_data_frame(
  metaweb_edges,
  directed = TRUE,
  vertices = nodes_fig1 %>%
    transmute(name = group, domain, n_taxa, x, y)
)

p_fig1 <- ggraph(g_fig1, layout="manual", x=x, y=y) +
  geom_edge_fan(aes(edge_colour=interaction), width=0.75,
                arrow=arrow(length=unit(3,"mm"), ends="first", type="closed")) +
  geom_node_point(aes(size=n_taxa), shape=21, fill=NA, colour= NA, stroke=0.4) +
  geom_node_text(aes(label=name), repel=TRUE, size=3) +
  scale_size(range=c(3,10), guide="none") +
  scale_edge_colour_manual(values=c(
    basal_plant="#195228", herbivory="#02703b", detritivory="#805204",
    bacterivory="#c9cc04", fungivory="#80b50d", invertivory="#edce02", vertivory="#edaa02")) +
  coord_equal() + theme_void() + theme(legend.position="right") +
  guides(edge_colour=guide_legend(title="Interaction type"))
p_fig1

ggsave(file.path("Plots","indicative_network.png"), p_fig1, width=7.5, height=6, dpi=300)

# #_____________________________________
# # Site-Wide Co-occurrence Networks
# #_____________________________________
# ## Split presence by site
# verts_site <- verts_pa %>%
#   filter(presence > 0) %>%
#   distinct(site, species) %>%
#   rename(node = species)
# inverts_site <- inverts_pa %>%
#   filter(presence > 0) %>%
#   distinct(site, taxon_id) %>%
#   rename(node = taxon_id)
# fungi_site <- fungi_pa %>%
#   filter(presence > 0) %>%
#   semi_join(
#     edges_metaweb %>% dplyr::select(from),
#     by = c("genus" = "from")
#   ) %>%
#   distinct(site, genus) %>%
#   rename(node = genus)
# 
# # full site-level presence
# nodes_site <- bind_rows(
#   verts_site,
#   inverts_site,
#   fungi_site
# ) %>%
#   distinct(site, node)
# # add basal nodes
# basal_site <- expand_grid(
#   site = unique(nodes_site$site),
#   node = basal_nodes$node_id)
# # merge
# nodes_site <- bind_rows(
#   nodes_site,
#   basal_site
# ) %>%
#   distinct(site, node)
# 
# ## Filter by site co-occurrence
# edges_site <- edges_metaweb %>%
#   inner_join(
#     nodes_site,
#     by = c("from" = "node")
#   ) %>%
#   rename(site_from = site) %>%
#   inner_join(
#     nodes_site,
#     by = c("to" = "node")
#   ) %>%
#   rename(site_to = site) %>%
#   filter(site_from == site_to) %>%
#   transmute(
#     site = site_from,
#     from,
#     to,
#     interaction,
#     rule
#   )
# # check number of edges and nodes per site
# edges_site %>% count(site)
# nodes_site %>%
#   count(site, name = "n_nodes")
# nodes_site %>%
#   mutate(
#     type = case_when(
#       node %in% fungi_pa$genus ~ "fungi",
#       node %in% inverts_pa$taxon_id ~ "inverts",
#       node %in% verts_pa$species ~ "verts",
#       node %in% basal_nodes$node_id ~ "basal",
#       TRUE ~ "other"
#     )
#   ) %>%
#   count(site, type)
# 
# ## Visualise networks on igraph
# ## (with PreyAveragedTL positioning)
# site_id <- 10
# 
# edges_s <- edges_site %>%
#   filter(site == site_id)
# 
# node_levels <- tibble(
#   type = c("verts", "inverts", "fungi", "basal"),
#   y    = c(4, 3, 2, 1)   # top → bottom
# )
# 
# plot_site_metaweb <- function(site_id) {
#   
#   nodes_s <- nodes_site %>%
#     dplyr::filter(site == site_id) %>%
#     dplyr::mutate(
#       node = as.character(node),
#       type = dplyr::case_when(
#         node %in% verts_pa$species    ~ "verts",
#         node %in% inverts_pa$taxon_id ~ "inverts",
#         node %in% fungi_pa$genus      ~ "fungi",
#         node %in% basal_nodes$node_id ~ "basal",
#         TRUE ~ NA_character_
#       )
#     ) %>%
#     dplyr::filter(!is.na(type)) %>%
#     dplyr::mutate(
#       y = dplyr::case_when(
#         type == "verts"   ~ 4,
#         type == "inverts" ~ 3,
#         type == "fungi"   ~ 2,
#         type == "basal"   ~ 1
#       )
#     ) %>%
#     dplyr::distinct(node, .keep_all = TRUE) %>%
#     dplyr::group_by(type) %>%
#     dplyr::mutate(
#       x = if (dplyr::n() == 1) 0 else seq(-1, 1, length.out = dplyr::n())
#     ) %>%
#     dplyr::ungroup()
#   
#   edges_s <- edges_site %>%
#     dplyr::filter(site == site_id) %>%
#     dplyr::semi_join(nodes_s, by = c("from" = "node")) %>%
#     dplyr::semi_join(nodes_s, by = c("to"   = "node"))
#   
#   nodes_s <- nodes_s %>%
#     dplyr::filter(node %in% c(edges_s$from, edges_s$to))
#   
#   g <- igraph::graph_from_data_frame(
#     edges_s %>% dplyr::select(from, to),
#     directed = TRUE,
#     vertices = nodes_s %>%
#       dplyr::transmute(name = node, type, x, y)
#   )
#   
#   ggraph::ggraph(g, layout = "manual", x = x, y = y) +
#     ggraph::geom_edge_link(alpha = 0.25, colour = "grey50", width = 0.3) +
#     ggraph::geom_node_point(ggplot2::aes(shape = type), size = 2) +
#     ggplot2::scale_y_continuous(
#       limits = c(0.5, 4.5),
#       breaks = c(1, 2, 3, 4),
#       labels = c("Basal", "Fungi", "Invertebrates", "Vertebrates")
#     ) +
#     ggplot2::scale_shape_manual(
#       values = c(
#         basal   = 15,
#         fungi   = 16,
#         inverts = 17,
#         verts   = 18
#       ),
#       drop = FALSE
#     ) +
#     ggplot2::theme_void() +
#     ggplot2::ggtitle(paste("Site", site_id))
# }
# 
# plot_site_metaweb(8)
