library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(vegan)
library(RColorBrewer)
library(patchwork)
# ______________________________________________________________________________
dun_dat = read.csv('vert_matrix.csv', row.names = 1)
blank_association = read.csv('blank_association.csv')

library(dplyr)
library(tidyr)

# Prepare read matrix (keep unassigned)
read_mat <- dun_dat[, colnames(dun_dat) != "taxonomy"]

# Long format for LOD filtering
dun_l <- as.data.frame(as.table(as.matrix(read_mat)))
colnames(dun_l) <- c("species", "sample", "reads")
dun_l$reads <- as.numeric(dun_l$reads)

# Blank handling + LOD
blank_reads <- dun_l[dun_l$sample %in% blank_association$blank, ]
blank_stats <- aggregate(reads ~ species, blank_reads,
                         function(x) c(mean = mean(x), sd = sd(x)))
blank_stats <- do.call(data.frame, blank_stats)
colnames(blank_stats) <- c("species", "mean_blank_reads", "sd_blank_reads")
blank_stats$LOD <- blank_stats$mean_blank_reads + 3 * blank_stats$sd_blank_reads

non_blank_data <- dun_l[!dun_l$sample %in% blank_association$blank, ]
lod_data <- merge(non_blank_data, blank_stats, by = "species", all.x = TRUE)
lod_data$reads[!is.na(lod_data$LOD) & lod_data$reads < lod_data$LOD] <- 0

# Rebuild wide matrix
filtered_wide <- lod_data %>%
  dplyr::select(species, sample, reads) %>%
  pivot_wider(names_from = sample, values_from = reads, values_fill = 0) %>%
  as.data.frame()
rownames(filtered_wide) <- filtered_wide$species
filtered_wide$species <- NULL

# Drop control samples
blanks <- c("NEG","POS","EB","BK","FB")
filtered_wide <- filtered_wide[, !grepl(paste(blanks, collapse="|"), colnames(filtered_wide))]

# 0.1% relative abundance threshold
filtered_wide[t(t(filtered_wide) / colSums(filtered_wide)) < 0.001] <- 0
filtered_wide <- filtered_wide[rowSums(filtered_wide) > 0, ]

# Remove duplicated dilution samples
remove_samples <- c("D_01_B1d","D_01_P1","D_02_B1","D_02_B2","D_02_P1d",
                    "D_02_P2","D_03_B1","D_03_B2","D_03_P1","D_03_P2",
                    "D_04_B2","D_04_B3d","D_04_P1","D_05_B1d","D_05_P1")
filtered_wide <- filtered_wide[, !(colnames(filtered_wide) %in% remove_samples)]

# Attach taxonomy and remove non-target taxa
tax <- dun_dat$taxonomy; names(tax) <- rownames(dun_dat)
filtered_wide$taxonomy <- tax[rownames(filtered_wide)]
excluded_taxa <- c("unassigned","Actinopteri","Hyperoartia",
                   "Homo_sapiens","Castor_fiber","Canis_lupus","Bos_taurus")
filtered_wide <- filtered_wide[!grepl(paste(excluded_taxa, collapse="|"),
                                      filtered_wide$taxonomy), ]

# Taxonomic reassignment + aggregation
filtered_wide$species <- rownames(filtered_wide)
filtered_wide <- filtered_wide[filtered_wide$species != "Passeriformes", ]
taxa_map <- c("Anatidae"="Anas_crecca",
              "Columba"="Columba_palumbus",
              "Fringilla"="Fringilla_coelebs")
filtered_wide$species <- ifelse(filtered_wide$species %in% names(taxa_map),
                                taxa_map[filtered_wide$species],
                                filtered_wide$species)

num_cols <- sapply(filtered_wide, is.numeric)
filtered_wide <- aggregate(filtered_wide[, num_cols],
                           by = list(species = filtered_wide$species),
                           sum)

rownames(filtered_wide) <- filtered_wide$species
filtered_wide$species <- NULL

# Final clean + export
filtered_wide <- filtered_wide[rowSums(filtered_wide) > 0, ]
filtered_wide <- filtered_wide[, sort(colnames(filtered_wide))]
write.csv(filtered_wide, "dun_filtered.csv")

# Acoustic data ----------------------------------------------------------

acoust <- read.csv('all_acoustics.csv')

acoust$species <- gsub(" ", "_", acoust$species)

# Ensure datetime is POSIXct
acoust <- acoust %>%
  mutate(
    datetime = dmy_hms(paste(date, time)),
    site_day = paste0("site", site, "_", format(datetime, "%Y-%m-%d"))
  )

# Define the anchor time
anchor_time <- ymd_hm("2024-06-19 15:00")

# Rolling 24-hour sample window
acoust <- acoust %>%
  mutate(
    # Calculate number of 24h windows since anchor
    window_num = floor(as.numeric(difftime(datetime, anchor_time, units = "hours")) / 24) +1,
    # Only keep records after anchor_time
    window_num = ifelse(datetime < anchor_time, NA, window_num),
    # Label for site and window
    site_window = paste0("site", site, "_win", window_num)
  ) %>%
  filter(!is.na(window_num))  # Remove any records before anchor

# remove misidentified/ incorrect IDs
acoust <- acoust %>%
  filter(check %in% c("1", "2", "corrected"))

# # since audiomoths were collected in from 7:30am, remove last 24-hour period as not complete
# acoust <- acoust %>%
#   filter(window_num < 6)

# Create a presence/ absence table
presence_matrix <- acoust %>%
  distinct(species, site_window) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = site_window, values_from = present, values_fill = 0)

# Set species as rownames
presence_matrix <- as.data.frame(presence_matrix)
rownames(presence_matrix) <- presence_matrix$species
presence_matrix$species <- NULL
presence_matrix <- presence_matrix[, sort(names(presence_matrix))]

write.csv(presence_matrix, 'acoustics_pa.csv')

# Combine dataframes into a big presence absence matrix --------------------
# Convert DNA matrix to presence/absence and tidy format
dna_pa <- filtered_wide
dna_pa[dna_pa > 0] <- 1
dna_pa <- as.data.frame(dna_pa)
dna_pa$species <- rownames(dna_pa)
dna_long <- dna_pa %>%
  pivot_longer(
    cols = -species,
    names_to = "sample",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  dplyr::select(sample, species)

write.csv(dna_long, "tree_pa_long.csv", row.names = TRUE)

# Convert acoustic matrix (species × site_window) to long format
acoustic_long <- presence_matrix %>%
  as.data.frame() %>%
  mutate(species = rownames(.)) %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "present") %>%
  filter(present == 1) %>%
  dplyr::select(-present)

write.csv(acoustic_long, "acoustic_pa_long.csv", row.names = TRUE)

combined_long <- bind_rows(dna_long, acoustic_long)

# Turn back into a wide format
combined_pa <- combined_long %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = species, values_from = present, values_fill = 0)
combined_pa <- as.data.frame(combined_pa)
rownames(combined_pa) <- combined_pa$sample
combined_pa$sample <- NULL

write.csv(combined_pa, "tree_aud_pa.csv", row.names = TRUE)

# Make 3 different matrices, one to species, genus and family levels
# All acoustic data already to species-level

is_species <- function(x) grepl("^[A-Z][a-z]+_[a-z]+$", x)
is_genus   <- function(x) grepl("^[A-Z][a-z]+$", x)

# Species- level ------------------------------------------------------
species_only <- filtered_wide[is_species(rownames(filtered_wide)), ]
species_pa <- species_only
species_pa[species_pa > 0] <- 1

# Convert species_pa (species x sample) to long format
species_pa$species <- rownames(species_pa)
dna_long <- species_pa %>%
  pivot_longer(
    cols = -species,
    names_to = "sample",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  dplyr::select(sample, species)

# acoustic_long is already in long format (sample, species)
# Combine both
species_long <- bind_rows(dna_long, acoustic_long)

# Pivot back to wide format
species_pa <- species_long %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = species, values_from = present, values_fill = 0)
species_pa <- as.data.frame(species_pa)
rownames(species_pa) <- species_pa$sample
species_pa$sample <- NULL

write.csv(species_pa, "species_pa.csv", row.names = TRUE)

# Genus-level database -------------------------------------------
dna_pa <- as.data.frame(filtered_wide)
dna_pa[dna_pa > 0] <- 1
dna_pa$taxon <- rownames(dna_pa)

# Extract genus and collapse to genus level
dna_long <- dna_pa %>%
  mutate(genus = str_extract(taxon, "^[A-Z][a-z]+")) %>%
  pivot_longer(cols = -c(taxon, genus), names_to = "sample", values_to = "present") %>%
  filter(present == 1) %>%
  dplyr::select(genus, sample)

acoustic_pa <- as.data.frame(presence_matrix)
acoustic_pa[acoustic_pa > 0] <- 1
acoustic_pa$taxon <- rownames(acoustic_pa)

acoustic_long <- acoustic_pa %>%
  mutate(genus = str_extract(taxon, "^[A-Z][a-z]+")) %>%
  pivot_longer(cols = -c(taxon, genus), names_to = "sample", values_to = "present") %>%
  filter(present == 1) %>%
  dplyr::select(genus, sample)

combined_long <- bind_rows(dna_long, acoustic_long) %>%
  distinct()
combined_long <- combined_long %>%
  filter(!genus %in% c("Fringillidae", "Muscicapidae", "Phasianidae"))

genus_pa <- combined_long %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genus, values_from = present, values_fill = 0)

genus_pa <- as.data.frame(genus_pa)
rownames(genus_pa) <- genus_pa$sample
genus_pa$sample <- NULL

write.csv(genus_pa, "genus_pa.csv", row.names = TRUE)

# Family-level database

genus_family_map <- data.frame(
  genus = c("Alauda", "Anas", "Bufo", "Buteo", "Capreolus", "Cervus", "Columba", 
            "Dendrocopos", "Erithacus", "Fringilla", "Lepus", "Lissotriton", "Lyrurus", 
            "Martes", "Parus", "Phasianus", "Phylloscopus", "Prunella", "Rana", "Sorex", 
            "Sturnus", "Sus", "Troglodytes", "Turdus", "Vulpes", "Riparia", "Branta", 
            "Corvus", "Periparus", "Tringa", "Acanthis", "Actitis", "Gallinago", "Sylvia", 
            "Cyanistes", "Pyrrhula", "Garrulus", "Haematopus", "Spinus", "Certhia", 
            "Scolopax", "Pluvialis", "Saxicola", "Regulus", "Carduelis", "Ardea",
            "Motacilla", "Anser", "Aegithalos", "Anthus", "Loxia", "Muscicapa", "Strix", "Rallus"),
  family = c("Alaudidae", "Anatidae", "Bufonidae", "Accipitridae", "Cervidae", "Cervidae", 
             "Columbidae", "Picidae", "Muscicapidae", "Fringillidae", "Leporidae", "Salamandridae", 
             "Phasianidae", "Mustelidae", "Paridae", "Phasianidae", "Phylloscopidae", "Prunellidae", 
             "Ranidae", "Soricidae", "Sturnidae", "Suidae", "Troglodytidae", "Turdidae", "Canidae", 
             "Hirundinidae", "Anatidae", "Corvidae", "Paridae", "Scolopacidae", "Fringillidae", 
             "Scolopacidae", "Scolopacidae", "Sylviidae", "Paridae", "Fringillidae", "Corvidae", 
             "Haematopodidae", "Fringillidae", "Certhiidae", "Scolopacidae", "Charadriidae", 
             "Muscicapidae", "Regulidae", "Fringillidae", "Ardeidae", "Motacillidae", "Anatidae", 
             "Aegithalidae", "Motacillidae", "Fringillidae", "Muscicapidae", "Strigidae", "Rallidae")
)

mammal_genus_family_map <- data.frame(
  genus = c("Myotis", "Pipistrellus", "Plecotus", "Sorex"),
  family = c("Vespertilionidae", "Vespertilionidae",
             "Vespertilionidae", "Soricidae")
)

genus_family_map <- bind_rows(
  genus_family_map,
  mammal_genus_family_map
) %>%
  distinct(genus, .keep_all = TRUE)


genus_pa <- genus_pa[, !is.na(colnames(genus_pa)), drop = FALSE]

family_pa <- as.data.frame(
  sapply(unique(colnames(genus_pa)), function(fam) {
    if (sum(colnames(genus_pa) == fam) == 1) {
      genus_pa[, fam]
    } else {
      apply(genus_pa[, colnames(genus_pa) == fam, drop = FALSE], 1, max)
    }
  })
)
rownames(family_pa) <- rownames(genus_pa)

families_to_update <- c("Fringillidae", "Muscicapidae", "Phasianidae")

for (fam in families_to_update) {
  if (fam %in% colnames(family_pa) && fam %in% rownames(dna_pa)) {
    vals <- as.character(dna_pa[fam, ])
    vals[!vals %in% c("1")] <- "0"
    dna_pa_vec <- as.numeric(vals)
    names(dna_pa_vec) <- colnames(dna_pa)
    
    common_samples <- intersect(rownames(family_pa), names(dna_pa_vec))
    old_vals <- as.numeric(family_pa[common_samples, fam])
    new_vals <- as.numeric(dna_pa_vec[common_samples])
    family_pa[common_samples, fam] <- pmax(old_vals, new_vals)
  }
}

write.csv(family_pa, "family_pa.csv", row.names = TRUE)

#----------------------------------------------------------------------------
# Correlation of site variables ----
dun_meta <- read.csv("dun_meta.csv")
library(dplyr)

site_df <- dun_meta %>%
  group_by(site) %>%
  summarise(
    age       = first(X2024_age),
    condition = first(condition),   # WCA score
    mean_ndvi = first(mean_ndvi),
    size_ha   = first(size),
    grazing   = first(grazing)
  ) %>%
  ungroup() %>%
  mutate(
    log_size = log(size_ha)
  )

vars <- site_df %>%
  dplyr::select(age, mean_ndvi, condition, log_size)

round(cor(vars, method = "spearman"), 2)

library(GGally)

ggpairs(
  vars,
  upper = list(continuous = wrap("cor", method = "spearman", size = 4)),
  lower = list(continuous = wrap("smooth", method = "loess", se = FALSE)),
  diag  = list(continuous = "densityDiag")
) +
  theme_bw(base_size = 11)

library(ggplot2)
library(ggpubr)

plot_correlation <- function(df, x, y,
                             xlab, ylab,
                             log_x = FALSE) {
  
  p <- ggplot(df, aes_string(x = x, y = y)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
    stat_cor(method = "spearman",
             label.x.npc = "left",
             label.y.npc = "top",
             size = 3) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11)
  
  if (log_x) {
    p <- p + scale_x_continuous(trans = "log10")
  }
  
  return(p)
}

p_age_ndvi <- plot_correlation(
  site_df, "age", "mean_ndvi",
  xlab = "Site age (years)",
  ylab = "Mean NDVI",
  log_x = TRUE
)

p_age_wca <- plot_correlation(
  site_df, "age", "condition",
  xlab = "Site age (years)",
  ylab = "Woodland condition (WCA)",
  log_x = TRUE
)

p_age_size <- plot_correlation(
  site_df, "age", "log_size",
  xlab = "Site age (years)",
  ylab = "Patch size (log ha)",
  log_x = TRUE
)
p_ndvi_wca <- plot_correlation(
  site_df, "mean_ndvi", "condition",
  xlab = "Mean NDVI",
  ylab = "Woodland condition (WCA)"
)

p_ndvi_size <- plot_correlation(
  site_df, "mean_ndvi", "log_size",
  xlab = "Mean NDVI",
  ylab = "Patch size (log ha)"
)

p_wca_size <- plot_correlation(
  site_df, "condition", "log_size",
  xlab = "Woodland condition (WCA)",
  ylab = "Patch size (log ha)"
)

library(patchwork)

collinearity_plots <-
  (p_age_ndvi | plot_spacer() | plot_spacer()) /
  (p_age_wca | p_ndvi_wca | plot_spacer()) /
  (p_age_size | p_ndvi_size | p_wca_size)

collinearity_plots

m_ndvi <- lm(mean_ndvi ~ log(age) + log_size + grazing, data = site_df)
summary(m_ndvi)

m_wca <- lm(condition ~ log(age) + log_size + grazing, data = site_df)
summary(m_wca)


# Exploratory richness plots ----

library(ggplot2)

species_pa <- read.csv("species_pa.csv", row.names = 1)
genus_pa <- read.csv("genus_pa.csv", row.names = 1)
family_pa <- read.csv("family_pa.csv", row.names = 1)

save_plot <- function(plot, filename) {
  dir.create("Exploratory Plots", showWarnings = FALSE)
  # ggsave(file.path("Exploratory Plots", filename), plot = plot, width = 5.5, height = 3, dpi = 300)
}

plot_median_richness <- function(pa_matrix, meta, sample_type, site_cols, tax_label = "Taxa") {
  # Filter metadata for sample type
  meta_type <- meta %>% filter(type == sample_type)
  # Subset PA matrix to matching samples
  pa_sub <- pa_matrix[rownames(pa_matrix) %in% meta_type$sample, , drop = FALSE]
  # Calculate richness per sample
  richness_df <- data.frame(
    sample = rownames(pa_sub),
    richness = rowSums(pa_sub > 0)
  )
  # Merge with metadata
  richness_df <- merge(richness_df, meta_type, by = "sample")
  # Plot
  p <- ggplot(richness_df, aes(x = factor(site), y = richness, color = factor(site))) +
    geom_jitter(width = 0.1, height = 0, size = 1.2, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 4) +
    labs(x = "Site", y = paste("Sample", tax_label, "Richness")) +
    scale_color_manual(values = site_cols, name = "Site") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  save_plot(p, paste0("median_richness_", tolower(tax_label), "_", sample_type, ".png"))
  return(p)
}

plot_total_richness <- function(pa_matrix, meta, sample_type, site_cols, tax_label = "Taxa") {
  # Filter metadata for sample type
  meta_type <- meta %>% filter(type == sample_type)
  # Subset PA matrix to matching samples
  pa_sub <- pa_matrix[rownames(pa_matrix) %in% meta_type$sample, , drop = FALSE]
  # Add site info as a column
  pa_sub$sample <- rownames(pa_sub)
  pa_sub <- left_join(pa_sub, meta_type[, c("sample", "site")], by = "sample")
  taxa_cols <- setdiff(colnames(pa_sub), c("sample", "site"))
  
  # For each site, get the number of taxa present in any sample from that site
  site_richness <- pa_sub %>%
    group_by(site) %>%
    summarise(
      richness = sum(colSums(across(all_of(taxa_cols)) > 0) > 0),
      .groups = "drop"
    )
  
  p <- ggplot(site_richness, aes(x = factor(site), y = richness, fill = factor(site))) +
    geom_col() +
    scale_fill_manual(values = site_cols, name = "Site") +
    labs(x = "Site", y = paste("Total", tax_label, "Richness")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  save_plot(p, paste0("total_richness_", tolower(tax_label), "_", sample_type, ".png"))
  return(p)
  
  return(print(site_richness))
}

n_sites <- length(unique(na.omit(dun_meta$site)))

site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(n_sites)
plot_median_richness(combined_pa, dun_meta, sample_type = "eDNA", site_cols, "Species")
# Median richness, species-level
plot_median_richness(species_pa, dun_meta, sample_type = "eDNA", site_cols, "Species")
plot_median_richness(species_pa, dun_meta, sample_type = "acoustic", site_cols, "Species")
# Total richness, species-level
plot_total_richness(combined_pa, dun_meta, sample_type = "eDNA", site_cols, "Species")
plot_total_richness(combined_pa, dun_meta, sample_type = "acoustic", site_cols, "Species")

# Median richness, genus-level
plot_median_richness(genus_pa, dun_meta, sample_type = "eDNA", site_cols, "Genus")
plot_median_richness(genus_pa, dun_meta, sample_type = "acoustic", site_cols, "Genus")
# Total richness, genus-level
plot_total_richness(genus_pa, dun_meta, sample_type = "eDNA", site_cols, "Genus")
plot_total_richness(genus_pa, dun_meta, sample_type = "acoustic", site_cols, "Genus")

# Median richness, family-level
plot_median_richness(family_pa, dun_meta, sample_type = "eDNA", site_cols, "Family")
plot_median_richness(family_pa, dun_meta, sample_type = "acoustic", site_cols, "Family")
# Total richness, family-level
plot_total_richness(family_pa, dun_meta, sample_type = "eDNA", site_cols, "Family")
plot_total_richness(family_pa, dun_meta, sample_type = "acoustic", site_cols, "Family")

# ============================================================
# Richness ~ Woodland condition (DNA vs acoustic)
# ============================================================

# Load WCA
wca <- read.csv("wca.csv") %>%
  rename_with(tolower) %>%
  mutate(site = as.character(site))

# Build per-sample richness for species-level matrix
rich_cond_df <- species_pa %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  mutate(richness = rowSums(across(-sample) > 0)) %>%
  dplyr::select(sample, richness) %>%
  left_join(
    dun_meta %>% 
      mutate(site = as.character(site)) %>% 
      dplyr::select(sample, site, type),
    by = "sample"
  ) %>%
  left_join(
    wca %>% dplyr::select(site, condition),
    by = "site"
  ) %>%
  filter(!is.na(condition), !is.na(type))

rich_cond_df <- rich_cond_df %>%
  mutate(
    site_f = factor(
      site,
      levels = as.character(sort(as.numeric(unique(site))))
    )
  )
library(RColorBrewer)
ylgn <- brewer.pal(9, "YlGn")
age_cols_9 <- rev(colorRampPalette(ylgn[2:9])(9))
age_cols <- c("#002e23", age_cols_9)
age_cols_named <- age_cols
names(age_cols_named) <- levels(rich_cond_df$site_f)

p_rich_condition <- ggplot(
  rich_cond_df,
  aes(x = condition, y = richness)
) +
  geom_boxplot(
    aes(fill = type, group = interaction(condition, type)),
    alpha = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(shape = type, colour = site_f),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    size = 2,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      eDNA     = "grey80",
      acoustic = "grey30"
    ),
    name = "Sample type"
  ) +
  scale_shape_manual(values = c(eDNA = 16, acoustic = 15)) +
  scale_colour_manual(
    values = age_cols_named,
    guide  = "none"
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Woodland condition",
    y = "Species richness",
    shape = "Sample type"
  ) +
  theme(legend.position = "bottom")


verts_wca_rich <- p_rich_condition

rich_cond_df <- rich_cond_df %>%
  mutate(
    cond_num = as.numeric(condition),   # Poor = 1, Moderate = 2, Good = 3
    site     = factor(site),
    type     = factor(type)
  )
library(lme4)
library(lmerTest)

m_rich_cond <- lmer(
  richness ~ cond_num * type + (1 | site),
  data = rich_cond_df
)

summary(m_rich_cond)
anova(m_rich_cond)

# ============================================================
# Plotting species richness vs site age and NDVI ----
# Mixed-effects models with CI
# ============================================================

library(dplyr)
library(ggplot2)
library(lme4)

#-------------------------------------------------------------
# Build richness dataframe
#-------------------------------------------------------------
ndvi_stats <- read.csv('verts_nvdi_stats.csv')

library(dplyr)
library(stringr)
ndvi_clean <- ndvi_stats %>%
  # keep only audiomoth rows
  filter(str_detect(Note, "^Audiomoth_")) %>%
  
  # extract site number
  mutate(
    site = str_extract(Note, "(?<=Audiomoth_)\\d+"),
    site = as.character(as.integer(site))  # drop leading zeros
  ) %>%
  
  # summarise NDVI per site
  group_by(site) %>%
  summarise(
    mean_ndvi = mean(X_mean, na.rm = TRUE),
    .groups = "drop"
  )

richness_df <- data.frame(
  sample   = rownames(combined_pa),
  richness = rowSums(combined_pa > 0)
) %>%
  left_join(
    dun_meta %>%
      mutate(site = as.character(site)) %>%
      dplyr::select(sample, site, type, condition, X2024_age),
    by = "sample"
  ) %>%
  left_join(
    ndvi_clean %>% mutate(site = as.character(site)),
    by = "site"
  ) %>%
  mutate(
    log_site_age = log(X2024_age + 1),
    site = factor(site),                             # SAFE for modelling
    type = factor(type, levels = c("acoustic", "eDNA"))
  ) %>%
  filter(
    !is.na(site),
    !is.na(log_site_age)
  )

# enforce numeric site ordering (THIS is the key)
site_levels <- as.character(sort(as.numeric(unique(richness_df$site))))
richness_df <- richness_df %>%
  mutate(site = factor(site, levels = site_levels))

library(RColorBrewer)

site_cols_rich <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(length(site_levels))
names(site_cols_rich) <- site_levels
#-------------------------------------------------------------
# 1. Richness ~ log(site age)
#-------------------------------------------------------------

m_rich_age <- lmer(
  richness ~ log_site_age * type + (1 | site),
  data = richness_df
)

new_age <- expand.grid(
  log_site_age = seq(
    min(richness_df$log_site_age, na.rm = TRUE),
    max(richness_df$log_site_age, na.rm = TRUE),
    length.out = 100
  ),
  type = levels(richness_df$type)
)

X_age <- model.matrix(~ log_site_age * type, new_age)
beta_age <- fixef(m_rich_age)
V_age <- vcov(m_rich_age)

new_age$fit <- as.numeric(X_age %*% beta_age)
new_age$se  <- sqrt(diag(X_age %*% V_age %*% t(X_age)))

new_age <- new_age %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

# Plot: richness vs log(site age)
p_rich_logage <- ggplot(
  richness_df,
  aes(
    x = log_site_age,
    y = richness,
    colour = site,
    shape  = type
  )
) +
  geom_ribbon(
    data = new_age,
    aes(
      x = log_site_age,
      ymin = lwr,
      ymax = upr,
      fill = type
    ),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  geom_line(
    data = new_age,
    aes(
      x = log_site_age,
      y = fit,
      linetype = type
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  geom_point(alpha = 0.75, size = 2,
             position = position_jitter(width = 0.03, height = 0)) +
  scale_linetype_manual(
    values = c(acoustic = "solid", eDNA = "dashed"),
    name = "Sample type"
  ) +
  scale_fill_manual(
    values = c(acoustic = "black", eDNA = "grey40"),
    name = "Sample type"
  ) +
  scale_colour_manual(values = site_cols, name = "Site") +
  scale_shape_manual(
    values = c(acoustic = 15, eDNA = 16),
    name = "Sample type"
  ) +
  labs(
    x = "Log(site age + 1)",
    y = "Species richness"
  ) +
  theme_bw(base_size = 12)

p_rich_logage

qqnorm(residuals(m_rich_age))
qqline(residuals(m_rich_age))

plot(fitted(m_rich_age),
     residuals(m_rich_age))
abline(h = 0, lty = 2)

#-------------------------------------------------------------
# 2. Richness ~ NDVI
#-------------------------------------------------------------
rich_ndvi <- richness_df %>%
  filter(!is.na(mean_ndvi))

m_rich_ndvi <- lmer(
  richness ~ mean_ndvi * type + (1 | site),
  data = rich_ndvi
)

# Prediction grid
new_ndvi <- expand.grid(
  mean_ndvi = seq(
    min(rich_ndvi$mean_ndvi, na.rm = TRUE),
    max(rich_ndvi$mean_ndvi, na.rm = TRUE),
    length.out = 100
  ),
  type = levels(richness_df$type)
)

X_ndvi <- model.matrix(~ mean_ndvi * type, new_ndvi)
beta_ndvi <- fixef(m_rich_ndvi)
V_ndvi <- vcov(m_rich_ndvi)

new_ndvi$fit <- as.numeric(X_ndvi %*% beta_ndvi)
new_ndvi$se  <- sqrt(diag(X_ndvi %*% V_ndvi %*% t(X_ndvi)))

new_ndvi <- new_ndvi %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

# Plot: richness vs NDVI
p_rich_ndvi <- ggplot(
  rich_ndvi,
  aes(
    x = mean_ndvi,
    y = richness,
    colour = site,
    shape  = type
  )
) +
  geom_ribbon(
    data = new_ndvi,
    aes(
      x = mean_ndvi,
      ymin = lwr,
      ymax = upr,
      fill = type
    ),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  geom_line(
    data = new_ndvi,
    aes(
      x = mean_ndvi,
      y = fit,
      linetype = type
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  geom_point(
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.005, height = 0)) +
  scale_colour_manual(values = site_cols, name = "Site") +
  scale_shape_manual(
    values = c(acoustic = 15, eDNA = 16),
    name = "Sample type"
  ) +
  scale_linetype_manual(
    values = c(acoustic = "solid", eDNA = "dashed"),
    name = "Sample type"
  ) +
  scale_fill_manual(
    values = c(acoustic = "black", eDNA = "grey40"),
    name = "Sample type"
  ) +
  labs(
    x = "Mean NDVI",
    y = "Species richness"
  ) +
  theme_bw(base_size = 12)

p_rich_ndvi

qqnorm(residuals(m_rich_ndvi))
qqline(residuals(m_rich_ndvi))

plot(fitted(m_rich_ndvi),
     residuals(m_rich_ndvi))
abline(h = 0, lty = 2)

# ggsave('Exploratory plots/richness_age_verts_plot.png', plot = p_rich_logage, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/richness_ndvi_verts_plot.png', plot = p_rich_ndvi, height = 4, width = 5, dpi = 300, bg = "white")

rich_age_verts <- p_rich_logage
rich_ndvi_verts <- p_rich_ndvi

# Model results
fixef(m_rich_age)
fixef(m_rich_ndvi)

library(lmerTest)
anova(m_rich_age)
anova(m_rich_ndvi)

summary(m_rich_age)
summary(m_rich_ndvi)

library(performance)
r2(m_rich_age)
r2(m_rich_ndvi)

## Species accumulation -------------
library(iNEXT)
library(ggplot2)
library(dplyr)

get_method_taxa <- function(pa_matrix, meta, method) {
  method_samples <- meta$sample[meta$type == method]
  method_samples <- intersect(method_samples, rownames(pa_matrix))
  
  pa_sub <- pa_matrix[method_samples, , drop = FALSE]
  colSums(pa_sub > 0) > 0
}

make_pa_matrix_method <- function(pa_matrix, meta, method, samples = NULL) {
  
  # Samples for this method
  if (is.null(samples)) {
    samples <- meta$sample[meta$type == method]
  }
  samples <- intersect(samples, rownames(pa_matrix))
  
  pa_sub <- pa_matrix[samples, , drop = FALSE]
  
  # Presence–absence
  pa_sub[pa_sub > 0] <- 1
  
  # Keep only taxa detected by this method
  keep_taxa <- colSums(pa_sub) > 0
  pa_sub <- pa_sub[, keep_taxa, drop = FALSE]
  
  t(pa_sub)  # taxa × samples
}

# eDNA species accumulation curves
# Get valid eDNA samples
edna_samples <- intersect(
  dun_meta$sample[dun_meta$type == "eDNA"],
  rownames(combined_pa)
)
# Build matrix
edna_matrix_all <- make_pa_matrix(combined_pa, edna_samples)

# Run iNEXT (extrapolate to 100 samples)
edna_overall_out <- iNEXT(
  list("eDNA" = edna_matrix_all),
  datatype = "incidence_raw",
  endpoint = 200,
  conf = 0.95,
  se = TRUE
)

# Plot
p_edna_overall <- ggiNEXT(
  edna_overall_out,
  type = 1,
  color.var = "Order.q"
) +
  scale_colour_manual(
    values = c(
      "0" = "#2da191"
    ),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(
      "0" = scales::alpha("#2da191", 0.25)
    ),
    guide = "none"
  ) +
  labs(
    x = "Number of eDNA samples",
    y = "Species richness",
    title = "eDNA species accumulation (all sites pooled)"
  ) +
  theme_minimal()
p_edna_overall

ggsave('Exploratory plots/p_edna_overall.png', plot = p_edna_overall, height = 5, width = 6, dpi = 300, bg = "white")

# seperate curves by site
## eDNA species accumulation by site ----------------------------------
library(iNEXT)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales)

## ---------- helpers ----------

make_pa_matrix_method <- function(pa_matrix, meta, method, samples) {
  samps <- intersect(samples, rownames(pa_matrix))
  pa_sub <- pa_matrix[samps, , drop = FALSE]
  pa_sub[pa_sub > 0] <- 1
  pa_sub <- pa_sub[, colSums(pa_sub) > 0, drop = FALSE]
  t(pa_sub)  # taxa × samples
}

is_valid_incidence_raw <- function(mat) {
  ncol(mat) >= 2 &&
    nrow(mat) >= 2 &&
    sum(rowSums(mat) > 1) >= 2
}

## ---------- build site matrices ----------

edna_meta <- dun_meta %>% filter(type == "eDNA")
edna_site_list <- split(edna_meta$sample, edna_meta$site)

edna_site_mats <- lapply(edna_site_list, function(samps) {
  make_pa_matrix_method(
    pa_matrix = combined_pa,
    meta      = dun_meta,
    method    = "eDNA",
    samples   = samps
  )
})

edna_site_mats <- edna_site_mats[
  vapply(edna_site_mats, is_valid_incidence_raw, logical(1))
]

## ---------- run iNEXT ----------

edna_site_out <- iNEXT(
  edna_site_mats,
  datatype = "incidence_raw",
  endpoint = max(40),
  conf = 0.95,
  se = TRUE
)

# ## ---------- PLOTTING  ----------
# 
# # IMPORTANT: derive site levels ONLY from ggiNEXT output
# # site_levels <- sort(as.numeric(
# #   unique(ggiNEXT(edna_site_out)$data$Assemblage)
# # ))
# # 
# # ylgn <- brewer.pal(9, "YlGn")
# # 
# # site_cols <- c(
# #   "#002e23",
# #   rev(colorRampPalette(ylgn[2:9])(length(site_levels) - 1))
# # )
# # names(site_cols) <- site_levels
# 
# site_fills <- alpha(site_cols, 0.25)
# 
# p_edna_site <- ggiNEXT(
#   edna_site_out,
#   type = 1,
#   color.var = "Assemblage"
# )
# 
# p_edna_site$data$Assemblage <- factor(
#   p_edna_site$data$Assemblage,
#   levels = site_levels
# )
# 
# p_edna_site <- p_edna_site +
#   scale_colour_manual(
#     values = site_cols,
#     name = "Site"
#   ) +
#   scale_fill_manual(
#     values = site_fills,
#     name = "Site"
#   ) +
#   scale_shape_manual(values = rep(16, length(site_levels))) +
#   guides(shape = "none") +
#   labs(
#     x = "Number of eDNA samples",
#     y = "Species richness",
#     title = "eDNA species accumulation by site"
#   ) +
#   theme_minimal()
# 
# p_edna_site
# 
# ggsave(
#   "Exploratory plots/p_edna_site.png",
#   plot = p_edna_site,
#   height = 5,
#   width = 6,
#   dpi = 300,
#   bg = "white"
# )
# 
# # Acoustics species accumulation curves
# acoust_samples <- intersect(
#   dun_meta$sample[dun_meta$type == "acoustic"],
#   rownames(species_pa)
# )
# 
# acoust_matrix_all <- make_pa_matrix(species_pa, acoust_samples)
# 
# acoust_overall_out <- iNEXT(
#   list("Acoustic" = acoust_matrix_all),
#   datatype = "incidence_raw",
#   endpoint = 500,
#   conf = 0.95,
#   se = TRUE
# )
# 
# p_acoust_overall <- ggiNEXT(
#   acoust_overall_out,
#   type = 1,
#   color.var = "Order.q"
# ) +
#   scale_colour_manual(
#     values = c(
#       "0" = "#2E7D32"
#     ),
#     guide = "none"
#   ) +
#   scale_fill_manual(
#     values = c(
#       "0" = scales::alpha("#2E7D32", 0.25)
#     ),
#     guide = "none"
#   ) +
#   labs(
#     x = "Number of acoustic sampling days",
#     y = "Species richness",
#     title = "Acoustic species accumulation (all sites pooled)"
#   ) +
#   theme_minimal()
# 
# p_acoust_overall
# ggsave('Exploratory plots/p_acoust_overall.png', plot = p_acoust_overall, height = 5, width = 6, dpi = 300, bg = "white")
# 
# # seperate curves by site
# library(RColorBrewer)
# 
# acoust_meta <- dun_meta %>% filter(type == "acoustic")
# acoust_site_list <- split(acoust_meta$sample, acoust_meta$site)
# 
# acoust_site_mats <- lapply(acoust_site_list, function(samps) {
#   samps <- intersect(samps, rownames(species_pa))
#   make_pa_matrix(species_pa, samps)
# })
# 
# acoust_site_out <- iNEXT(
#   acoust_site_mats,
#   datatype = "incidence_raw",
#   endpoint = 30,
#   conf = 0.95,
#   se = TRUE
# )
# 
# site_levels <- sort(as.numeric(
#   unique(ggiNEXT(acoust_site_out)$data$Assemblage)
# ))
# 
# ylgn <- brewer.pal(9, "YlGn")
# 
# site_cols <- c(
#   "#002e23",
#   rev(colorRampPalette(ylgn[2:9])(length(site_levels) - 1))
# )
# names(site_cols) <- site_levels
# 
# site_fills <- alpha(site_cols, 0.25)
# 
# p_acoust_site <- ggiNEXT(
#   acoust_site_out,
#   type = 1,
#   color.var = "Assemblage"
# )
# 
# p_acoust_site$data$Assemblage <- factor(
#   p_acoust_site$data$Assemblage,
#   levels = site_levels
# )
# 
# p_acoust_site <- p_acoust_site +
#   scale_colour_manual(
#     values = site_cols,
#     name = "Site"
#   ) +
#   scale_fill_manual(
#     values = site_fills,
#     name = "Site"
#   ) +
#   scale_shape_manual(values = rep(16, length(site_levels))) +
#   guides(shape = "none") +
#   labs(
#     x = "Number of acoustic sampling days",
#     y = "Species richness",
#     title = "Acoustic species accumulation by site"
#   ) +
#   theme_minimal()
# 
# p_acoust_site
# 
# ggsave('Exploratory plots/p_acoust_site.png', plot = p_acoust_site, height = 5, width = 6, dpi = 300, bg = "white")
## Jaccard dissimilarity -------------

run_nmds_and_plot <- function(pa_matrix, meta, site_cols,
                              title_text, type = NULL,
                              exclude_samples = NULL) {
  
  if (!is.null(type)) {
    meta <- meta %>% filter(type == !!type)
  }
  if (!is.null(exclude_samples)) {
    meta <- meta %>% filter(!sample %in% exclude_samples)
  }
  
  shared_samples <- intersect(rownames(pa_matrix), meta$sample)
  pa_sub <- pa_matrix[shared_samples, , drop = FALSE]
  meta_sub <- meta %>% filter(sample %in% shared_samples)
  
  # NMDS
  nmds <- metaMDS(
    pa_sub,
    distance = "jaccard",
    k = 2,
    trymax = 100,
    autotransform = FALSE
  )
  
  # Extract scores
  nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
  nmds_scores$sample <- rownames(nmds_scores)
  nmds_scores <- left_join(nmds_scores, meta_sub, by = "sample")
  
  # Plot
  p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(
      aes(color = factor(site), shape = factor(type)),
      size = 3, shape = 15, alpha = 0.95
    ) +
    scale_color_manual(values = site_cols, name = "Site") +
    scale_shape_manual(values = c(16, 17), name = "Type") +
    labs(
      title = title_text,
      subtitle = paste("Stress =", round(nmds$stress, 3))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  clean_title <- gsub("[^a-zA-Z0-9]", "_", tolower(title_text))
  clean_type <- if (!is.null(type)) paste0("_", type) else ""
  filename <- paste0("nmds_", clean_title, clean_type, ".png")
  save_plot(p, filename)
  
  list(
    plot   = p,
    nmds   = nmds,
    stress = nmds$stress,
    pa_mat = pa_sub,
    meta   = meta_sub
  )
}

nmds_acoust <- run_nmds_and_plot(species_pa, dun_meta, site_cols, "NMDS: Acoustics", type = "acoustic")
nmds_edna <- run_nmds_and_plot(combined_pa, dun_meta, site_cols, "NMDS: eDNA Only", type = "eDNA")
# 
# # Genus
# run_nmds_and_plot(genus_pa, dun_meta, site_cols, "NMDS: All Types", exclude_samples = "D_07_B1")
# run_nmds_and_plot(genus_pa, dun_meta, site_cols,"NMDS: eDNA Only",
#                   type = "eDNA")

# ggsave('Exploratory plots/nmds_acoust_plot.png', plot = nmds_acoust, height = 4, width = 5, dpi = 300, bg = "white")

ggsave('Exploratory plots/nmds_edna_plot.png', plot = nmds_edna$plot, height = 4, width = 5, dpi = 300, bg = "white")

nmds_acoust$stress

# Distance matrix (same data as NMDS)
dist_jaccard <- vegdist(
  nmds_acoust$pa_mat,
  method = "jaccard",
  binary = TRUE
)

# Metadata (aligned)
meta_perm <- nmds_acoust$meta %>%
  mutate(site = factor(site))

rownames(meta_perm) <- meta_perm$sample

# PERMANOVA
adonis2(
  dist_jaccard ~ site,
  data = meta_perm,
  permutations = 999
)

# Distance matrix (same data as NMDS)
dist_jaccard <- vegdist(
  nmds_edna$pa_mat,
  method = "jaccard",
  binary = TRUE
)

# Metadata (aligned)
meta_perm <- nmds_edna$meta %>%
  mutate(site = factor(site))

rownames(meta_perm) <- meta_perm$sample

# PERMANOVA
adonis2(
  dist_jaccard ~ site,
  data = meta_perm,
  permutations = 999
)
# Functional analysis --------------------------------------------
library(traitdata)

data("pantheria")
data("amphibio")
data("amniota")

species_vec <- gsub("_", " ", colnames(species_pa))
# 
# # ---- Mammals: PanTHERIA ----
# mammal_traits <- pantheria %>%
#   dplyr::filter(scientificNameStd %in% species_vec) %>%
#   dplyr::rename(scientificName = scientificNameStd) %>%
#   dplyr::select(scientificName, dplyr::everything()) %>%
#   dplyr::mutate(database = "pantheria")
# 
# # ---- Mammals: EltonTraits ----
# mammal_traits2 <- elton_mammals %>%
#   dplyr::filter(scientificNameStd %in% species_vec) %>%
#   dplyr::rename(scientificName = scientificNameStd) %>%
#   dplyr::select(scientificName, dplyr::everything()) %>%
#   dplyr::mutate(database = "elton_mammals")
# 
# # ---- Split species names once ----
# species_split <- strsplit(species_vec, " ")
# species_df <- data.frame(
#   genus = sapply(species_split, `[`, 1),
#   species = sapply(species_split, `[`, 2),
#   stringsAsFactors = FALSE
# )
# 
# bird_traits <- amniota %>%
#   dplyr::filter(scientificNameStd %in% species_vec) %>%
#   dplyr::rename(scientificName = scientificNameStd) %>%
#   dplyr::select(scientificName, dplyr::everything()) %>%
#   dplyr::mutate(database = "amniota")
# 
# bird_traits2 <- avonet %>%
#   filter(Genus %in% species_df$genus & Species %in% species_df$species) %>%
#   mutate(scientificName = paste(Genus, Species)) %>%
#   dplyr::select(scientificName, everything())
# 
# bird_traits3 <- elton_birds %>%
#   filter(scientificNameStd %in% species_vec) %>%
#   rename(scientificName = scientificNameStd) %>%
#   dplyr::select(scientificName, everything())%>%
#   mutate(database = "elton_birds")
# 
# amphib_traits <- amphibio %>%
#   mutate(scientificName = paste(Genus, Species)) %>%
#   filter(scientificName %in% species_vec) %>%
#   dplyr::select(scientificName, everything()) %>%
#   mutate(database = "amphibio")
# 
# prefix_traits <- function(df) {
#   db <- unique(df$database)
#   stopifnot(length(db) == 1)
#   
#   df %>%
#     dplyr::rename_with(
#       ~ paste(db, ., sep = "__"),
#       -c(scientificName, database)
#     )
# }
# 
# mammal_traits_p   <- prefix_traits(mammal_traits)
# mammal_traits2_p  <- prefix_traits(mammal_traits2)
# bird_traits_p     <- prefix_traits(bird_traits)
# bird_traits3_p    <- prefix_traits(bird_traits3)
# amphib_traits_p   <- prefix_traits(amphib_traits)
# 
# all_traits <- bind_rows(
#   mammal_traits_p,
#   mammal_traits2_p,
#   bird_traits_p,
#   bird_traits3_p,
#   amphib_traits_p
# )
# 
# # Save to manually curate database
# write.csv(all_traits, 'traitdata_download.csv')

library(SYNCSA)
library(tibble)
library(stringr)

traits <- read.csv("vert_traitdata.csv")

names(site_cols) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

library(dplyr)
library(tibble)

# Clean names consistently
traits$scientificName <- str_replace_all(traits$scientificName, " ", "_")
colnames(species_pa)  <- str_trim(colnames(species_pa))

#### NOTE: FILTERING TO JUST ACOUSTIC SAMPLES DUE TO POOR RESOLUTION OF TREE ROLLERS

acoustic_samples <- dun_meta %>%
  filter(type == "acoustic") %>%
  pull(sample)

species_pa_acoustic <- species_pa[
  rownames(species_pa) %in% acoustic_samples,
  ,
  drop = FALSE
]

# Select traits used for functional diversity
trait_matrix <- traits %>%
  filter(scientificName %in% colnames(species_pa)) %>%
  arrange(factor(scientificName, levels = colnames(species_pa))) %>%
  dplyr::select(
    scientificName,
    body_mass_g,
    tl_test,
    hab_test
  ) %>%
  mutate(
    body_mass_g = log(as.numeric(body_mass_g) + 1),
    tl_test     = as.factor(tl_test),
    hab_test    = as.factor(hab_test)
  ) %>%
  column_to_rownames("scientificName")

# Remove samples with too few species
comm <- species_pa_acoustic[rowSums(species_pa_acoustic) > 2, ]

# Match species between traits and community
shared_species <- intersect(colnames(comm), rownames(trait_matrix))

comm_filt <- comm[, colSums(comm) >= 2, drop = FALSE]
trait_filt <- trait_matrix[colnames(comm_filt), , drop = FALSE]

## Rao's Quadratic Entropy and Functional Redundnacy --------------------------

rao_res <- rao.diversity(
  comm_filt,
  trait_filt,
  ord         = "metric",
  standardize = TRUE
)

rao_df <- data.frame(
  sample = rownames(comm_filt),
  RaoQ   = rao_res$FunRao,
  FRed   = rao_res$FunRedundancy
)

# Join with metadata
rao_df <- left_join(rao_df, dun_meta, by = "sample")

# Prep factors
rao_df <- rao_df %>%
  filter(type != "P") %>%
  mutate(
    type = factor(type, levels = c("acoustic", "eDNA")),
    site = factor(site, levels = names(site_cols))
  )

# Plot
raoq_plot <- ggplot(rao_df, aes(x = factor(site), y = RaoQ, color = factor(site))) +
  geom_jitter(width = 0.1, height = 0, size = 1.5, alpha = 0.7) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 3) + 
  labs(x = "Site", y = "Rao’s Q") +
  scale_color_manual(values = site_cols, name = "Site") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

fred_plot <- ggplot(rao_df, aes(x = factor(site), y = FRed, color = factor(site))) +
  geom_jitter(width = 0.1, height = 0, size = 1.5, alpha = 0.7) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 3) + 
  labs(x = "Site", y = "FRed") +
  scale_color_manual(values = site_cols, name = "Site") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

# ggsave('Exploratory plots/raoq_plot.png', plot = raoq_plot, height = 5, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/fred_plot.png', plot = fred_plot, height = 5, width = 5, dpi = 300, bg = "white")

## Functional Evenness --------------------------
library(FD)
feve_res <- dbFD(
  x         = trait_filt,
  a         = comm_filt,
  calc.FRic = FALSE,
  calc.FDiv = FALSE,
  calc.CWM  = FALSE,
  stand.x   = FALSE,
  corr      = "cailliez"
)

feve_df <- data.frame(
  sample = rownames(comm_filt),
  FEve   = feve_res$FEve
)

feve_df <- feve_df %>%
  left_join(dun_meta, by = "sample") %>%
  mutate(
    site = factor(site, levels = names(site_cols)),
    type = factor(type, levels = c("acoustic", "eDNA"))
  )

feve_plot <- ggplot(feve_df, aes(x = site, y = FEve, color = site)) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 3) +
  labs(x = "Site", y = "Functional Evenness (FEve)") +
  scale_color_manual(values = site_cols) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

#_____________________________________________________________________
# Combine functional diversity metrics
#_____________________________________________________________________
cond_cols <- c(
  "Poor"     = "#d73027",
  "Moderate" = "#fee08b",
  "Good"     = "#1a9850"
)

div_metrics <- rao_df %>%
  dplyr::select(
    sample,
    RaoQ,
    FRed,
    site,
    type
  ) %>%
  left_join(
    feve_df %>% dplyr::select(sample, FEve),
    by = "sample"
  ) %>%
  left_join(
    dun_meta %>% dplyr::select(sample, condition, X2024_age),
    by = "sample"
  ) %>%
  left_join(
    ndvi_clean,
    by = "site"
  ) %>%
  mutate(
    # log site age
    log_site_age = log(X2024_age + 1),
    # categorical woodland condition (interpretable)
    condition_cat = case_when(
      condition < 26                     ~ "Poor",
      condition >= 26 & condition <= 35  ~ "Moderate",
      condition > 35                     ~ "Good"
    ),
    condition_cat = factor(
      condition_cat,
      levels = c("Poor", "Moderate", "Good")
    ),
    site = factor(site, levels = names(site_cols)),
    type = factor(type, levels = c("acoustic", "eDNA"))
  )


# Site age plots
library(lme4)
library(emmeans)

plot_age <- function(yvar, ylab, log_response = FALSE) {
  
  resp <- if (log_response) paste0("log1p(", yvar, ")") else yvar
  
  # Fit model
  m <- lmer(
    as.formula(paste(resp, "~ log_site_age + (1 | site)")),
    data = div_metrics
  )
  
  # Prediction grid
  newdat <- distinct(div_metrics, log_site_age)|>
    arrange(log_site_age)
  
  if (!log_response) {
    ## ---- RAW SCALE (FEve) ----
    X <- model.matrix(~ log_site_age, newdat)
    beta <- fixef(m)
    V <- vcov(m)
    
    newdat$fit <- as.numeric(X %*% beta)
    newdat$se  <- sqrt(diag(X %*% V %*% t(X)))
    newdat$lwr <- newdat$fit - 1.96 * newdat$se
    newdat$upr <- newdat$fit + 1.96 * newdat$se
    
  } else {
    ## ---- LOG + BACK-TRANSFORM (RaoQ / FRed) ----
    emm <- emmeans(
      m,
      specs = ~ log_site_age,
      at    = list(log_site_age = newdat$log_site_age),
      type  = "response"
    )
    
    pred <- as.data.frame(emm)
    
    newdat$fit <- pred$response
    newdat$lwr <- pred$lower.CL
    newdat$upr <- pred$upper.CL
  }
  
  p <- ggplot(div_metrics, aes(x = log_site_age, y = .data[[yvar]])) +
    geom_ribbon(
      data = newdat,
      aes(x = log_site_age, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = newdat,
      aes(x = log_site_age, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.85,
      position = position_jitter(width = 0.03)
    ) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(x = "Log(site age + 1)", y = ylab, color = "Site") +
    theme(legend.position = "none")
  
  list(plot = p, model = m)
}

vert_age_feve <- plot_age("FEve", "Functional Evenness (FEve)")
vert_age_rao  <- plot_age("RaoQ", "Rao’s Q", log_response = TRUE)
vert_age_fred <- plot_age("FRed", "Functional Redundancy", log_response = TRUE)

qqnorm(residuals(vert_age_feve$model))
qqline(residuals(vert_age_feve$model))

plot(fitted(vert_age_feve$model),
     residuals(vert_age_feve$model))
abline(h = 0, lty = 2)

summary(vert_age_fred$model)

# NDVI plots
plot_ndvi <- function(yvar, ylab, log_response = FALSE) {
  
  resp <- if (log_response) paste0("log1p(", yvar, ")") else yvar
  
  m <- lmer(
    as.formula(paste(resp, "~ mean_ndvi + (1 | site)")),
    data = div_metrics
  )
  
  newdat <- distinct(div_metrics, mean_ndvi)|>
    arrange(mean_ndvi)
  
  if (!log_response) {
    ## ---- RAW SCALE (FEve) ----
    X <- model.matrix(~ mean_ndvi, newdat)
    beta <- fixef(m)
    V <- vcov(m)
    
    newdat$fit <- as.numeric(X %*% beta)
    newdat$se  <- sqrt(diag(X %*% V %*% t(X)))
    newdat$lwr <- newdat$fit - 1.96 * newdat$se
    newdat$upr <- newdat$fit + 1.96 * newdat$se
    
  } else {
    ## ---- LOG + BACK-TRANSFORM (RaoQ / FRed) ----
    emm <- emmeans(
      m,
      specs = ~ mean_ndvi,
      at    = list(mean_ndvi = newdat$mean_ndvi),
      type  = "response"
    )
    
    pred <- as.data.frame(emm)
    
    newdat$fit <- pred$response
    newdat$lwr <- pred$lower.CL
    newdat$upr <- pred$upper.CL
  }
  
  p <- ggplot(div_metrics, aes(x = mean_ndvi, y = .data[[yvar]])) +
    geom_ribbon(
      data = newdat,
      aes(x = mean_ndvi, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = newdat,
      aes(x = mean_ndvi, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.85,
      position = position_jitter(width = 0.005)
    ) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(x = "Mean NDVI", y = ylab, color = "Site") +
    theme(legend.position = "none")
  
  list(plot = p, model = m)
}

vert_ndvi_feve <- plot_ndvi("FEve", "Functional Evenness (FEve)")
vert_ndvi_rao  <- plot_ndvi("RaoQ", "Rao’s Q", log_response = TRUE)
vert_ndvi_fred <- plot_ndvi("FRed", "Functional Redundancy", log_response = TRUE)

qqnorm(residuals(vert_ndvi_fred$model))
qqline(residuals(vert_ndvi_fred$model))

plot(fitted(vert_ndvi_fred$model),
     residuals(vert_ndvi_fred$model))
abline(h = 0, lty = 2)

summary(vert_ndvi_rao$model)

# Model results
fixef(vert_age_rao$model)
fixef(vert_age_fred$model)
fixef(vert_age_feve$model)

library(lmerTest)
anova(vert_age_rao$model)
anova(vert_age_fred$model)
anova(vert_age_feve$model)

summary(vert_age_rao$model)
summary(vert_age_fred$model)
summary(vert_age_feve$model)

library(performance)
r2(vert_age_rao$model)
r2(vert_age_fred$model)
r2(vert_age_feve$model)


# ggsave('Exploratory plots/raoq_ndvi_plot.png', plot = p_ndvi_rao$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/feve_ndvi_plot.png', plot = p_ndvi_feve$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/fred_ndvi_plot.png', plot = p_ndvi_fred$plot, height = 4, width = 5, dpi = 300, bg = "white")
# 
# ggsave('Exploratory plots/raoq_age_plot.png', plot = p_age_rao$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/feve_age_plot.png', plot = p_age_feve$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/fred_age_plot.png', plot = p_age_fred, height = 4, width = 5, dpi = 300, bg = "white")

# Woodland condition plots
plot_condition <- function(yvar, ylab) {
  ggplot(div_metrics, aes(x = condition_cat, y = .data[[yvar]], fill = condition_cat)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8, color = "black") +
    scale_fill_manual(values = cond_cols) +
    theme_bw(base_size = 12) +
    labs(x = "Woodland condition", y = ylab) +
    theme(legend.position = "none")
}

p_cond_feve <- plot_condition("FEve", "Functional Evenness (FEve)")
p_cond_rao  <- plot_condition("RaoQ", "Rao’s Q")
p_cond_fred <- plot_condition("FRed", "Functional Redundancy")

# PCoA of traits -----
library(cluster)
library(MASS)
library(ggplot2)
library(dplyr)
library(fastDummies)
library(RColorBrewer)

# Gower distance on traits
gower_dist <- daisy(trait_filt, metric = "gower")

pcoa_res <- cmdscale(gower_dist, eig = TRUE, k = 2)

species_scores <- as.data.frame(pcoa_res$points)
colnames(species_scores) <- c("PC1", "PC2")
species_scores$species <- rownames(species_scores)

# Axis labels
var_expl <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
pc1_lab <- paste0("PC1 (", var_expl[1], "%)")
pc2_lab <- paste0("PC2 (", var_expl[2], "%)")

# Ensure column order matches
comm_mat <- comm_filt[, species_scores$species, drop = FALSE]
comm_mat <- as.matrix(comm_mat)

sample_scores <- comm_mat %*% as.matrix(species_scores[, c("PC1", "PC2")])
sample_scores <- as.data.frame(sample_scores)
colnames(sample_scores) <- c("PC1", "PC2")
sample_scores$sample <- rownames(sample_scores)

sample_scores <- sample_scores %>%
  left_join(
    dun_meta %>% dplyr::select(sample, site, condition, type, X2024_age),
    by = "sample"
  ) %>%
  mutate(
    site = factor(site),
    condition = factor(condition),
    type = factor(type, levels = c("acoustic", "eDNA")),
    x2024_age = as.numeric(X2024_age)
  )

library(fastDummies)

trait_num <- trait_filt %>%
  mutate(
    body_mass_g = as.numeric(body_mass_g),
    trophic_level = as.character(tl_test),
    hab_test = as.character(hab_test)
  ) %>%
  fastDummies::dummy_cols(
    select_columns = c("trophic_level", "hab_test"),
    remove_first_dummy = FALSE,   # 👈 KEEP ALL
    remove_selected_columns = TRUE
  ) %>%
  dplyr::select(where(is.numeric))


trait_loadings <- cor(
  trait_num,
  scale(species_scores[, c("PC1", "PC2")]),
  use = "pairwise.complete.obs"
)

trait_loadings_df <- as.data.frame(trait_loadings)
trait_loadings_df$trait <- rownames(trait_loadings_df)

trait_loadings_df <- trait_loadings_df %>%
  filter(sqrt(PC1^2 + PC2^2) > 0.3)

plot_radius <- 0.85 * min(
  diff(range(sample_scores$PC1)),
  diff(range(sample_scores$PC2))
)

trait_loadings_df$PC1_scaled <- trait_loadings_df$PC1 * plot_radius
trait_loadings_df$PC2_scaled <- trait_loadings_df$PC2 * plot_radius

old_sites   <- c("1", "2", "3")
young_sites <- c("7", "8", "10")

old_scores   <- sample_scores %>% filter(site %in% old_sites)
young_scores <- sample_scores %>% filter(site %in% young_sites)

# Kernal Density differences
library(MASS)

x_range <- range(
  c(sample_scores$PC1, trait_loadings_df$PC1_scaled),
  na.rm = TRUE
)

y_range <- range(
  c(sample_scores$PC2, trait_loadings_df$PC2_scaled),
  na.rm = TRUE
)

bw_x <- bandwidth.nrd(sample_scores$PC1)
bw_y <- bandwidth.nrd(sample_scores$PC2)

k_old <- kde2d(
  x = old_scores$PC1,
  y = old_scores$PC2,
  h = c(bw_x, bw_y),
  n = 200,
  lims = c(x_range, y_range)
)

k_young <- kde2d(
  x = young_scores$PC1,
  y = young_scores$PC2,
  h = c(bw_x, bw_y),
  n = 200,
  lims = c(x_range, y_range)
)

k_diff <- k_old
k_diff$z <- k_old$z - k_young$z

kdf <- expand.grid(PC1 = k_diff$x, PC2 = k_diff$y)
kdf$z <- as.vector(k_diff$z)

# Plot
# KD colouring - for 3 sites as reference
breaks <- seq(-0.9, 0.5, length.out = 18)
zero_bin <- findInterval(0, breaks, rightmost.closed = TRUE)
n_bins <- length(breaks) - 1

# number of bins on each side
n_neg <- zero_bin - 1
n_pos <- n_bins - zero_bin

cols_pal <- c(
  "#d95a04", "#de6408", "#e46e0d", "#ea7913", "#ee8620", "#f19331",
  "#f4a247", "#f7b464", "#f9c98b", "#fbdcb4", "#fdeee0", "#fff6f0",
  "#ffffff",
  "#e6f5ef", "#bfe6d6", "#8fd3ba", "#5fbfa0", "#2f9f78", "#1b7f5e",
  "#004d3a"
)
cols_pal <- c(
  "#c99600", "#d4a300", "#deaf16", "#e8bc36", "#f0c95a",
  "#f4d57a", "#f7dfa0", "#fae8bd", "#fdf2d8",
  "#ffffff", "#e6f8f1", "#b3e2cf",
  "#6ec7a8", "#3aa887", "#1f8668", "#0f7053", "#004d3a"
)

library(RColorBrewer)

site_cols <- colorRampPalette(
  c("#000000", "#ededed")  # black → light grey
)(
  length(unique(sample_scores$site))
)

# ensure ordering: oldest (site 1) = black
names(site_cols) <- sort(unique(sample_scores$site))

x_lims <- range(kdf$PC1, na.rm = TRUE)
y_lims <- range(kdf$PC2, na.rm = TRUE)

p_pcoa <- ggplot() +
  geom_contour_filled(
    data = kdf,
    aes(PC1, PC2, z = z),
    breaks = breaks
  ) +
  scale_fill_manual(values = cols_pal, name = "Old – Young\nDensity") +
  # geom_point(
  #   data = sample_scores,
  #   aes(PC1, PC2, color = site, shape = type),
  #   size = 2, alpha = 0.9
  # ) +
  # scale_color_manual(values = site_cols, name = "Site") +
  # for simplified, presentation plot
  geom_point(
    data = sample_scores,
    aes(PC1, PC2),
    colour = "black",
    shape = 16,
    size = 1.2,
    alpha = 0.6
  ) +
  geom_segment(
    data = trait_loadings_df,
    aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
    arrow = arrow(type = "closed", length = unit(0.18, "cm")),
    colour = "black", linewidth = 0.6
  ) +
  # geom_text(
  #   data = trait_loadings_df,
  #   aes(x = PC1_scaled, y = PC2_scaled, label = trait),
  #   size = 3.5, vjust = -0.2
  # ) +
  coord_cartesian(
    xlim = x_lims,
    ylim = y_lims,
    clip= "off",
    expand = FALSE
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme(panel.grid = element_blank())

p_pcoa

p_pcoa_vert <- p_pcoa

# ggsave('Exploratory plots/verts_pcoa_plot.png', plot = p_pcoa, height = 4.2, width = 6, dpi = 300, bg = "white")

library(dplyr)
library(ggplot2)
library(mgcv)
#____________________________________________________
# 1. Define reference mature woodland sites
#____________________________________________________
ref_sites <- c("1", "2", "3")

#____________________________________________________
# 2. Compute centroid of mature woodland in PCoA space
#____________________________________________________
centroid <- sample_scores %>%
  filter(site %in% ref_sites) %>%
  summarise(
    centroid_PC1 = mean(PC1, na.rm = TRUE),
    centroid_PC2 = mean(PC2, na.rm = TRUE)
  )

# Print to confirm
print(centroid)

#____________________________________________________
# 3. Compute functional distance of every sample to centroid
#---------------------------------------------------------
sample_scores <- sample_scores %>%
  mutate(
    dist_to_old = sqrt(
      (PC1 - centroid$centroid_PC1)^2 +
        (PC2 - centroid$centroid_PC2)^2
    ),
    x2024_age = as.numeric(as.character(x2024_age))
  )
# Identify the oldest age
last_age <- max(sample_scores$x2024_age, na.rm = TRUE)

# Plot-only grouping: split ONLY the last age by site
sample_scores <- sample_scores %>%
  mutate(
    box_group = case_when(
      x2024_age == last_age & site == "2" ~ "left",
      x2024_age == last_age & site == "1" ~ "right",
      TRUE ~ "all"
    ),
    box_group = factor(box_group, levels = c("left", "right", "all"))
  )
#____________________________________________________
# 4. Plot functional recovery as boxplots
#____________________________________________________
p_recovery_box <- ggplot(sample_scores,
                         aes(
                           x = factor(x2024_age),
                           y = dist_to_old,
                           fill = site,
                           group = interaction(x2024_age, box_group))) +
  geom_boxplot(
    alpha = 0.75,
    outlier.shape = NA,
    position = position_dodge2(width = 0.75, preserve = "single")
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    size = 2,
    alpha = 0.6,
    colour = "black"
  ) +
  scale_fill_manual(values = site_cols) +
  scale_colour_manual(values = site_cols) +
  theme_bw(base_size = 14) +
  labs(
    x = "Woodland Age (years)",
    y = "Functional Distance"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

p_recovery_box
vert_recovery <- p_recovery_box

# ggsave('Exploratory plots/verts_recovery_plot.png', plot = p_recovery, height = 3, width = 8, dpi = 300, bg = "white")

# Plot species onto the PCoA for indicative species silhouettes
library(ggrepel)
x_lims_sp <- range(species_scores$PC1, na.rm = TRUE)
y_lims_sp <- range(species_scores$PC2, na.rm = TRUE)

plot_radius_sp <- 0.9 * min(
  diff(x_lims_sp),
  diff(y_lims_sp)
)

trait_loadings_df <- trait_loadings_df %>%
  mutate(
    PC1_sp = PC1 * plot_radius_sp,
    PC2_sp = PC2 * plot_radius_sp
  )

p_species <- ggplot() +
  
  # species points
  geom_point(
    data = species_scores,
    aes(PC1, PC2),
    size = 2,
    colour = "black"
  ) +
  
  # repelled species labels with leader lines
  geom_text_repel(
    data = species_scores,
    aes(PC1, PC2, label = species),
    size = 3,
    seed = 123,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4
  ) +
  
  # trait arrows
  geom_segment(
    data = trait_loadings_df,
    aes(x = 0, y = 0, xend = PC1_sp, yend = PC2_sp),
    arrow = arrow(type = "closed", length = unit(0.12, "cm")),
    linewidth = 0.6,
    colour = "black"
  ) +
  # trait labels (also repelled)
  geom_text_repel(
    data = trait_loadings_df,
    aes(PC1_sp, PC2_sp, label = trait),
    size = 3.5,
    seed = 123,
    fontface = "italic",
    box.padding = 0.4,
    segment.color = "black",
    segment.size = 0.4
  ) +
  coord_equal(
    xlim = x_lims_sp,
    ylim = y_lims_sp,
    clip = "off",
    expand = FALSE
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank()
  ) +
  
  labs(
    title = "Species positions in vertebrate trait space",
    x = pc1_lab,
    y = pc2_lab
  )

p_species

## Null models: centroid shift & KDE difference ----

# Create old / young grouping variable
sample_scores <- sample_scores %>%
  mutate(
    age_group = case_when(
      site %in% old_sites   ~ "old",
      site %in% young_sites ~ "young",
      TRUE                  ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_group))

### Null model 1 - centroid shift ----
# Observed centroid shift
cent_old <- colMeans(
  sample_scores[sample_scores$age_group == "old", c("PC1", "PC2")]
)

cent_young <- colMeans(
  sample_scores[sample_scores$age_group == "young", c("PC1", "PC2")]
)

obs_centroid_shift <- sqrt(sum((cent_old - cent_young)^2))
obs_centroid_shift

# Permutation test
set.seed(42)
n_perm <- 1000
perm_centroid_shift <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  
  perm_labels <- sample(sample_scores$age_group)
  
  cent_old_p <- colMeans(
    sample_scores[perm_labels == "old", c("PC1", "PC2")]
  )
  
  cent_young_p <- colMeans(
    sample_scores[perm_labels == "young", c("PC1", "PC2")]
  )
  
  perm_centroid_shift[i] <- sqrt(sum((cent_old_p - cent_young_p)^2))
}

# One-tailed p-value
p_centroid <- (sum(perm_centroid_shift >= obs_centroid_shift) + 1) / (n_perm + 1)
p_centroid
max(perm_centroid_shift)
obs_centroid_shift
summary(perm_centroid_shift)

### Null model 2 - KDE Difference ----
# Observed KDE difference statistic
obs_kde_stat <- sum(abs(k_diff$z))
obs_kde_stat

# KDE permutation function
kde_diff_stat_perm <- function(scores, labels, bw_x, bw_y, x_range, y_range) {
  
  old_p   <- scores[labels == "old", ]
  young_p <- scores[labels == "young", ]
  
  k_old_p <- kde2d(
    x = old_p$PC1,
    y = old_p$PC2,
    h = c(bw_x, bw_y),
    n = 200,
    lims = c(x_range, y_range)
  )
  
  k_young_p <- kde2d(
    x = young_p$PC1,
    y = young_p$PC2,
    h = c(bw_x, bw_y),
    n = 200,
    lims = c(x_range, y_range)
  )
  
  sum(abs(k_old_p$z - k_young_p$z))
}

# Run KDE permutation test
set.seed(42)
perm_kde_stat <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  
  perm_labels <- sample(sample_scores$age_group)
  
  perm_kde_stat[i] <- kde_diff_stat_perm(
    scores  = sample_scores,
    labels  = perm_labels,
    bw_x    = bw_x,
    bw_y    = bw_y,
    x_range = x_range,
    y_range = y_range
  )
}

# One-tailed p-value
p_kde <- (sum(perm_kde_stat >= obs_kde_stat) + 1) / (n_perm + 1)
p_kde
max(perm_kde_stat)
obs_kde_stat
summary(perm_kde_stat)
