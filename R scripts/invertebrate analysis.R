#_____________________________________________________________________
# INVERTEBRATE 18S COMMUNITY DATA----
#_____________________________________________________________________

library(microeco)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(FD)
library(data.table)
library(stringr)
library(SYNCSA)
library(cluster)
library(scales)
library(reshape2)

#_____________________________________________________________________
## 1. Load and clean data ----
## remove occurrences of 2 or fewer reads, and remove ASVs with 5 or fewer reads across all samples
#_____________________________________________________________________

# Community table
inverts <- read.delim("18S_ASV_table_with_taxonomy.tsv",
                      check.names = FALSE, sep = "\t")

# Metadata
meta <- read.csv("soil_meta.csv") |>
  rename_with(tolower)

# Fix error in sample naming (same logic as ITS2)
fix_last <- c("06"="14","07"="15","08"="16")
cn <- colnames(inverts)
parts <- strsplit(cn,"-",fixed=TRUE)
first <- sapply(parts,`[`,1)
mid   <- sapply(parts,`[`,2)
last  <- sapply(parts,`[`,3)
suffix <- paste(mid,last,sep="-")
cn_new <- cn
safe_num <- suppressWarnings(as.numeric(first))
for(suf in unique(suffix)){
  idx <- which(suffix==suf)
  if(length(idx)>1){
    keep  <- idx[which.min(safe_num[idx])]
    dupes <- setdiff(idx, keep)
    for(i in dupes){
      old_last <- last[i]
      if(old_last %in% names(fix_last)){
        new_last   <- fix_last[[old_last]]
        cn_new[i]  <- paste(first[i], mid[i], new_last, sep="-")
      }
    }
  }
}
colnames(inverts) <- cn_new

# Identify correct sample columns (numeric + not tax/confidence columns)
non_sample_cols <- c("ASV_id", "taxonomy", "confidence",
                     paste0("rank_", 1:7),
                     "phylum","class","order","family","genus","species")

sample_cols <- setdiff(colnames(inverts), non_sample_cols)
sample_cols <- sample_cols[sapply(inverts[sample_cols], is.numeric)]

# Everything that is not a sample column = tax columns
tax_cols <- setdiff(colnames(inverts), sample_cols)


# Per sample threshold (2 reads)
inverts[ , sample_cols] <- lapply(inverts[ , sample_cols], function(x) {
  x[x <= 2] <- 0
  x
})

# Global threshold (5 reads)
global_totals <- rowSums(inverts[ , sample_cols])
inverts_filtered <- inverts[global_totals >= 5, ]

# Sort sample col names (drop first segment, join with "_")
new_sample_names <- sapply(sample_cols, function(x) {
  parts <- strsplit(x, "-", fixed = TRUE)[[1]]
  parts <- parts[-1]                 # drop the first segment
  paste(parts, collapse = "_")       # join remaining with _
})

# Apply new names only to sample columns
colnames(inverts_filtered)[colnames(inverts_filtered) %in% sample_cols] <- new_sample_names

# Clean taxonomy names (taxonomy + rank_1:rank_7 strings)
prefix_pattern <- "[dkpcofgs]__"

inverts_filtered[ , tax_cols] <- lapply(inverts_filtered[ , tax_cols], function(x) {
  gsub(prefix_pattern, "", x)
})

# Check read depth per sample
is_sample <- sapply(inverts_filtered, is.numeric)
sample_cols <- names(inverts_filtered)[is_sample]

read_depth <- colSums(inverts_filtered[, sample_cols])

read_depth_df <- data.frame(
  sample = names(read_depth),
  reads  = as.numeric(read_depth)
)

# remove replicated samples (from poor amplification the first time)
# Identify sample columns again
is_sample <- sapply(inverts_filtered, is.numeric)
sample_cols <- names(inverts_filtered)[is_sample]

# Extract base names (remove trailing r)
base_names <- gsub("r$", "", sample_cols, ignore.case = TRUE)

sample_info <- data.frame(
  original = sample_cols,
  base = base_names,
  depth = colSums(inverts_filtered[, sample_cols, drop = FALSE]),
  stringsAsFactors = FALSE
)

# Pick highest-depth sample for each base name
best_samples <- sample_info %>%
  group_by(base) %>%
  slice_max(order_by = depth, n = 1, with_ties = FALSE) %>%
  ungroup()

keep_samples <- best_samples$original

# Keep only these samples + taxonomy columns
inverts_filtered <- inverts_filtered[, c(setdiff(colnames(inverts_filtered), sample_cols), keep_samples)]

# Rename kept samples to remove "r"
clean_names <- gsub("r$", "", keep_samples, ignore.case = TRUE)
colnames(inverts_filtered)[match(keep_samples, colnames(inverts_filtered))] <- clean_names

# Recalculate read depths after replicate removal
is_sample <- sapply(inverts_filtered, is.numeric)
sample_cols <- names(inverts_filtered)[is_sample]

read_depth <- colSums(inverts_filtered[, sample_cols])
read_depth_df <- data.frame(
  sample = names(read_depth),
  reads = as.numeric(read_depth)
)

# remove poor-read depth samples (< 5000)
low_depth_samples <- read_depth_df$sample[read_depth_df$reads < 5000]
low_depth_samples

# filter to just metazoan taxa
inverts_filtered <- inverts_filtered %>%
  mutate(
    rank_2_clean = trimws(as.character(rank_2))
  )

metazoan_phyla <- c(
  "Arthropoda",
  "Annelida",
  "Nematozoa",
  "Rotifera",
  "Tardigrada",
  "Gastrotricha",
  "Platyhelminthes",
  "Mollusca"
)

inverts_filtered %>%
  mutate(is_metazoa = rank_2_clean %in% metazoan_phyla) %>%
  count(is_metazoa)

inverts_filtered <- inverts_filtered %>%
  filter(rank_2_clean %in% metazoan_phyla)

inverts_filtered <- inverts_filtered[, !(colnames(inverts_filtered) %in% low_depth_samples)]

write.csv(inverts_filtered, "18S_filtered.csv", row.names = FALSE)

#_____________________________________________________________________
# INVERTEBRATE COMMUNITY DATA — TRAITS + INITIAL EXPLORATION (ASV-LEVEL) ----
#_____________________________________________________________________

# Load bespoke trait database
traits_raw <- read.csv("dundreggan_invert_taxonomy_refine.csv",
                       check.names = FALSE)

# (traits already have phylum/class/order/family/genus/species without prefixes)
# Keep one row per ASV_id
traits_raw <- traits_raw %>% distinct(ASV_id, .keep_all = TRUE)

#_____________________________________________________________________
## TRAIT MERGE (ASV-level) ----
#_____________________________________________________________________
# Prepare clean trait reference
invert_traits_ref <- traits_raw %>%
  mutate(
    microhabitat   = as.character(microhabitat),
    feeding_guild  = as.character(feeding_guild),
    moisture_pref  = as.character(moisture_pref)
  )

# Merge traits onto filtered ASV table by ASV_id
merged_traits <- inverts_filtered %>%
  left_join(invert_traits_ref, by = "ASV_id")

# Remove .x columns and clean .y suffixes
merged_traits <- merged_traits %>%
  dplyr::select(-ends_with(".x")) %>%
  rename_with(~ gsub("\\.y$", "", .x), ends_with(".y"))

# LOG-TRANSFORM BODY SIZE HERE
merged_traits <- merged_traits %>%
  dplyr::select(-dplyr::ends_with(".y")) %>%
  dplyr::rename_with(~ sub("\\.x$", "", .x), dplyr::ends_with(".x"))

merged_traits <- merged_traits %>%
  mutate(log_body_size = log10(av_body_size + 1e-6))

# consolidate functional trait units ----
merged_traits <- merged_traits %>%
  mutate(
    functional_confidence = tolower(trimws(as.character(functional_confidence))),
    species = na_if(trimws(as.character(species)), ""),
    genus   = na_if(trimws(as.character(genus)), ""),
    family  = na_if(trimws(as.character(family)), ""),
    order   = na_if(trimws(as.character(order)), "")
  ) %>%
  mutate(
    functional_unit = case_when(
      functional_confidence == "s" & !is.na(species) ~ paste0("s__", species),
      functional_confidence == "g" & !is.na(genus)   ~ paste0("g__", genus),
      functional_confidence == "f" & !is.na(family)  ~ paste0("f__", family),
      functional_confidence == "o" & !is.na(order)   ~ paste0("o__", order),
      TRUE ~ NA_character_
    )
  )

write.csv(merged_traits, "invertebrate_trait_database.csv", row.names = FALSE)
#______________________________________________________________________________
## Richness ----
#______________________________________________________________________________
# Load metadata and ensure sample names align
meta <- read.csv("soil_meta.csv")
colnames(meta) <- tolower(colnames(meta))

#_____________________________________________________________________--
# Identify TRUE sample columns using regex pattern
# (numeric + underscore/hyphen + numeric)
#_____________________________________________________________________--
sample_cols <- grep("^[0-9]+[_-][0-9]+$", colnames(merged_traits), value = TRUE)

# Exclude blanks (EB, NG)
sample_cols <- sample_cols[!grepl("EB|NG", sample_cols, ignore.case = TRUE)]

# create dataframe for network analysis
meta_net <- meta %>%
  rename_with(tolower) %>%
  mutate(sample = as.character(sample)) %>%
  dplyr::select(sample, site)

invert_pa_long <- merged_traits %>%
  dplyr::select(functional_unit, all_of(sample_cols)) %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "sample",
    values_to = "reads"
  ) %>%
  mutate(
    sample   = as.character(sample),
    presence = as.integer(reads > 0)
  ) %>%
  dplyr::filter(!is.na(functional_unit)) %>%
  dplyr::group_by(sample, functional_unit) %>%
  dplyr::summarise(presence = max(presence), .groups = "drop") %>%
  dplyr::left_join(meta_net, by = "sample") %>%
  dplyr::rename(taxon_id = functional_unit)
write.csv(invert_pa_long, "invert_pa_long.csv", row.names = FALSE)
#_____________________________________________________________________--
# Prepare ASV long table (presence/absence)
#_____________________________________________________________________--
asv_long <- merged_traits %>%
  dplyr::select(
    ASV_id, family, genus,
    all_of(sample_cols)
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "reads"
  ) %>%
  mutate(
    sample   = gsub("^X", "", sample),
    presence = ifelse(reads > 0, 1, 0)
  )

#_____________________________________________________________________--
# Join metadata
#_____________________________________________________________________--
asv_long <- left_join(asv_long, meta, by = "sample")

# Ensure factors
asv_long <- asv_long %>%
  mutate(
    site      = factor(site),
    tree.type = factor(tree.type)
  )

#_____________________________________________________________________--
# ASV richness per sample
#_____________________________________________________________________--
asv_richness <- asv_long %>%
  group_by(sample, site, tree.type) %>%
  summarise(asv_richness = sum(presence), .groups = "drop")

# Ensure factors again
asv_richness <- asv_richness %>%
  mutate(
    site      = factor(site),
    tree.type = factor(tree.type)
  )

# Remove samples with zero ASVs (optional but consistent with NMDS)
empty_samples <- asv_richness$sample[asv_richness$asv_richness == 0]

if(length(empty_samples) > 0) {
  message("Removing ", length(empty_samples), " empty samples: ",
          paste(empty_samples, collapse = ", "))
  
  asv_richness <- asv_richness %>% 
    filter(asv_richness > 0)
  
  asc_long <- asv_long %>%
    filter(!sample %in% empty_samples)
}
#_____________________________________________________________________--
# Colour palettes
#_____________________________________________________________________--
site_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(
  length(unique(asv_richness$site))
)

tree_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(
  length(unique(asv_richness$tree.type))
)

#_____________________________________________________________________--
# PLOTS
#_____________________________________________________________________--
# --- ASV richness per site
ggplot(asv_richness, aes(x = site, y = asv_richness, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "ASV richness") +
  theme(legend.position = "none")

# --- ASV richness per tree type
ggplot(asv_richness, aes(x = tree.type, y = asv_richness, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Tree type", y = "ASV richness") +
  theme(legend.position = "none")
#_____________________________________________________________________
# ASV Accumulation curves ----
#_____________________________________________________________________
library(iNEXT)
library(RColorBrewer)
# Build PA matrix: samples x ASVs
invert_pa <- asv_long %>%
  dplyr::select(sample, ASV_id, presence) %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from  = ASV_id,
    values_from = presence,
    values_fill = 0
  )
make_pa_matrix <- function(pa_matrix, samples) {
  pa_sub <- pa_matrix[rownames(pa_matrix) %in% samples, , drop = FALSE]
  pa_sub[pa_sub > 0] <- 1
  t(pa_sub)  # taxa × samples
}

invert_pa_mat <- invert_pa %>%
  column_to_rownames("sample") %>%
  as.matrix()

# overall richness
# Get valid soil eDNA samples
invert_samples <- intersect(
  meta$sample,
  rownames(invert_pa_mat)
)

invert_matrix_all <- make_pa_matrix(
  invert_pa_mat,
  invert_samples
)

invert_overall_out <- iNEXT(
  list("Soil invertebrates" = invert_matrix_all),
  datatype = "incidence_raw",
  endpoint = 1250,
  conf = 0.95,
  se = TRUE
)

p_invert_overall <- ggiNEXT(
  invert_overall_out,
  type = 1,
  color.var = "Order.q"
) +
  scale_colour_manual(
    values = c("0" = "#bd8f13"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("0" = scales::alpha("#bd8f13", 0.25)),
    guide = "none"
  ) +
  labs(
    x = "Number of soil samples",
    y = "ASV richness",
    title = "Soil invertebrate species accumulation (all sites pooled)"
  ) +
  theme_minimal()

p_invert_overall

ggsave('Exploratory plots/p_invert_overall.png', plot = p_invert_overall, height = 5, width = 6, dpi = 300, bg = "white")

#------------------------------------------------------------
# Separate curves by site
#------------------------------------------------------------
invert_meta <- meta %>%
  mutate(site = factor(site))

invert_site_list <- split(invert_meta$sample, invert_meta$site)

invert_site_mats <- lapply(invert_site_list, function(samps) {
  samps <- intersect(samps, rownames(invert_pa_mat))
  make_pa_matrix(invert_pa_mat, samps)
})

invert_site_out <- iNEXT(
  invert_site_mats,
  datatype = "incidence_raw",
  endpoint = 150,
  conf = 0.95,
  se = TRUE
)

site_levels <- sort(as.numeric(
  unique(ggiNEXT(invert_site_out)$data$Assemblage)
))

ylgn <- RColorBrewer::brewer.pal(9, "YlGn")

site_cols <- c(
  "#002e23",
  rev(colorRampPalette(ylgn[2:9])(length(site_levels) - 1))
)
names(site_cols) <- site_levels

site_fills <- scales::alpha(site_cols, 0.25)

p_invert_site <- ggiNEXT(
  invert_site_out,
  type = 1,
  color.var = "Assemblage"
)

p_invert_site$data$Assemblage <- factor(
  p_invert_site$data$Assemblage,
  levels = site_levels
)

p_invert_site <- p_invert_site +
  scale_colour_manual(
    values = site_cols,
    name = "Site"
  ) +
  scale_fill_manual(
    values = site_fills,
    name = "Site"
  ) +
  scale_shape_manual(values = rep(16, length(site_levels))) +
  guides(shape = "none") +
  labs(
    x = "Number of soil samples",
    y = "ASV richness",
    title = "Soil invertebrate species accumulation by site"
  ) +
  theme_minimal()

p_invert_site

ggsave('Exploratory plots/p_invert_site.png', plot = p_invert_site, height = 5, width = 6, dpi = 300, bg = "white")
#_____________________________________________________________________
# ASV RICHNESS ~ NDVI and SITE AGE
#_____________________________________________________________________

library(dplyr)
library(ggplot2)
library(lme4)
library(RColorBrewer)
library(performance)
library(broom.mixed)

#_____________________________________________________________________
# 1. Prepare metadata
#_____________________________________________________________________
meta2 <- meta %>%
  rename_with(tolower) %>%
  mutate(sample = as.character(sample))

#_____________________________________________________________________
# 2. Build ASV richness dataframe with predictors
#_____________________________________________________________________
rich_df <- asv_richness %>%
  mutate(sample = as.character(sample)) %>%
  left_join(
    meta2 %>% dplyr::select(sample, x2024_age, mean_ndvi),
    by = "sample"
  ) %>%
  mutate(
    site    = factor(site),
    log_age = log(x2024_age + 1)
  ) %>%
  filter(
    !is.na(mean_ndvi),
    !is.na(x2024_age),
    !is.na(asv_richness)
  )

wca <- read.csv("wca.csv") %>%
  rename_with(tolower) %>%
  mutate(site = as.character(site))

wca <- wca %>%
  mutate(site = as.character(site))

rich_df <- rich_df %>%
  mutate(site = as.character(site)) %>%
  left_join(wca, by = "site") %>%
  mutate(
    condition = factor(
      condition,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High")
    )
  )
rich_df <- rich_df %>%
  mutate(
    site = factor(
      site,
      levels = as.character(sort(as.numeric(unique(site))))
    )
  )
ylgn <- brewer.pal(9, "YlGn")
age_cols_9 <- rev(colorRampPalette(ylgn[2:9])(9))
age_cols <- c("#002e23", age_cols_9)

p_asv_cond <- ggplot(
  rich_df,
  aes(x = condition, y = asv_richness)
) +
  geom_boxplot(
    fill = "grey85",
    colour = "grey30",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = site),
    width = 0.15, size = 2, alpha = 0.8
  ) +
  scale_colour_manual(values = age_cols) +
  theme_bw(base_size = 12) +
  labs(
    x = "Woodland condition",
    y = "OTU richness"
  )

invert_wca_rich <- p_asv_cond

rich_df <- rich_df %>%
  mutate(
    cond_num = as.numeric(condition)  # Low = 1, Moderate = 2, High = 3
  )
library(lme4)
library(lmerTest)

m_cond_trend <- lmer(
  asv_richness ~ cond_num + (1 | site),
  data = rich_df
)
summary(m_cond_trend)

#_____________________________________________________________________
# 3. Colour palette (consistent with other figures)
#_____________________________________________________________________
site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(
  length(unique(rich_df$site))
)
#_____________________________________________________________________
# 4. ASV richness ~ NDVI
#_____________________________________________________________________
plot_lmm_rich_ndvi <- function(dat = rich_df) {
  
  # Fit model
  m <- lmer(asv_richness ~ mean_ndvi + (1 | site), data = dat)
  
  # Prediction grid
  newdat <- data.frame(
    mean_ndvi = seq(
      min(dat$mean_ndvi, na.rm = TRUE),
      max(dat$mean_ndvi, na.rm = TRUE),
      length.out = 100
    )
  )
  
  # Fixed-effect predictions + CI
  X    <- model.matrix(~ mean_ndvi, newdat)
  beta <- fixef(m)
  V    <- vcov(m)
  
  newdat$fit <- as.numeric(X %*% beta)
  newdat$se  <- sqrt(diag(X %*% V %*% t(X)))
  newdat$lwr <- newdat$fit - 1.96 * newdat$se
  newdat$upr <- newdat$fit + 1.96 * newdat$se
  
  # Plot
  p <- ggplot(dat, aes(x = mean_ndvi, y = asv_richness)) +
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
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.8
    ) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(
      x = "Mean NDVI",
      y = "ASV richness",
      color = "Site"
    )
  
  list(
    plot     = p,
    model    = m,
    summary  = summary(m),
    anova    = anova(update(m, . ~ . - mean_ndvi), m),
    confint  = confint(m, method = "Wald"),
    r2       = performance::r2(m),
    tidy_fix = broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE)
  )
}

p_rich_ndvi <- plot_lmm_rich_ndvi()
p_rich_ndvi$plot

#_____________________________________________________________________
# 5. ASV richness ~ SITE AGE (log-transformed)
#_____________________________________________________________________
plot_lmm_rich_age <- function(dat = rich_df) {
  
  # Fit model
  m <- lmer(asv_richness ~ log_age + (1 | site), data = dat)
  
  # Prediction grid
  newdat <- data.frame(
    log_age = seq(
      min(dat$log_age, na.rm = TRUE),
      max(dat$log_age, na.rm = TRUE),
      length.out = 100
    )
  )
  
  # Fixed-effect predictions + CI
  X    <- model.matrix(~ log_age, newdat)
  beta <- fixef(m)
  V    <- vcov(m)
  
  newdat$fit <- as.numeric(X %*% beta)
  newdat$se  <- sqrt(diag(X %*% V %*% t(X)))
  newdat$lwr <- newdat$fit - 1.96 * newdat$se
  newdat$upr <- newdat$fit + 1.96 * newdat$se
  
  # Plot
  p <- ggplot(dat, aes(x = log_age, y = asv_richness)) +
    geom_ribbon(
      data = newdat,
      aes(x = log_age, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = newdat,
      aes(x = log_age, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.8,
      position = position_jitter(width = 0.03, height = 0)
    ) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(
      x = "Site age (log)",
      y = "ASV richness",
      color = "Site"
    )
  
  list(
    plot     = p,
    model    = m,
    summary  = summary(m),
    anova    = anova(update(m, . ~ . - log_age), m),
    confint  = confint(m, method = "Wald"),
    r2       = performance::r2(m),
    tidy_fix = broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE)
  )
}

p_rich_age <- plot_lmm_rich_age()
p_rich_age$plot

# ggsave('Exploratory plots/richness_age_inverts_plot.png', plot = p_rich_age$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/richness_ndvi_fungi_plot.png', plot = p_rich_ndvi$plot, height = 4, width = 5, dpi = 300, bg = "white")

rich_age_inverts  <- p_rich_age$plot
rich_ndvi_inverts <- p_rich_ndvi$plot

qqnorm(residuals(p_rich_ndvi$model))
qqline(residuals(p_rich_ndvi$model))

plot(fitted(p_rich_ndvi$model), resid(p_rich_ndvi$model),
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2)

# No transformations neccessary

# Model results
fixef(p_rich_age$model)
fixef(p_rich_ndvi$model)

library(lmerTest)
anova(p_rich_age$model)
anova(p_rich_ndvi$model)

summary(p_rich_age$model)
summary(p_rich_ndvi$model)

library(performance)
r2(p_rich_age$model)
r2(p_rich_ndvi$model)
#______________________________________________________________________________
## NMDS ----
#______________________________________________________________________________
# 1️⃣ Create sample-by-ASV matrix (presence/absence)
asv_matrix <- asv_long %>%
  dplyr::select(sample, site, tree.type, ASV_id, presence) %>%
  group_by(sample, site, tree.type, ASV_id) %>%
  summarise(presence = as.integer(any(presence > 0)), .groups = "drop") %>%
  pivot_wider(names_from = ASV_id, values_from = presence, values_fill = 0)

asv_matrix <- asv_matrix[, c("sample","site","tree.type",
                             colnames(asv_matrix)[!(colnames(asv_matrix) %in% c("sample","site","tree.type"))] )]

# 2️⃣ Extract metadata and numeric matrix
meta_nmds <- asv_matrix %>% dplyr::select(sample, site, tree.type)
asv_mat <- asv_matrix %>% dplyr::select(-site, -tree.type)
asv_mat <- as.data.frame(asv_mat)
rownames(asv_mat) <- asv_mat$sample
asv_mat <- asv_mat %>% dplyr::select(-sample)

# 4️⃣ Run NMDS with Jaccard distance
set.seed(123)
nmds <- metaMDS(asv_mat, distance = "jaccard", k = 2,
                trymax = 100, autotransform = FALSE)

nmds$stress

# 5️⃣ Prepare NMDS scores + metadata for plotting
nmds_points <- as.data.frame(nmds$points) %>%
  rownames_to_column("sample") %>%
  left_join(meta_nmds, by = "sample") %>%
  mutate(
    tree.type = ifelse(is.na(tree.type), "Unknown", as.character(tree.type)),
    tree.type = factor(tree.type),
    site      = factor(site)
  )

# 6️⃣ Define colour and shape palettes
site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(
  length(unique(na.omit(nmds_points$site)))
)

shape_vals <- c(
  "Birch"   = 16,
  "Pine"    = 17,
  "None"    = 15
)

# 7️⃣ Plot NMDS
inverts_nmds <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = site, shape = tree.type)) +
  geom_point(size = 3, alpha = 0.95) +
  scale_color_manual(values = site_cols, na.value = "grey50") +
  scale_shape_manual(values = shape_vals) +
  theme_bw(base_size = 12) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Site",
    shape = "Tree type",
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# stats
dist_jaccard <- vegdist(asv_mat, method = "jaccard", binary = TRUE)

# Create metadata for PERMANOVA
meta_perm <- meta_nmds %>%
  filter(sample %in% rownames(asv_mat)) %>%
  mutate(
    site = factor(site),
    tree.type = factor(tree.type)
  )

meta_perm <- as.data.frame(meta_perm)
rownames(meta_perm) <- meta_perm$sample

adonis2(
  dist_jaccard ~ site + tree.type,
  data = meta_perm,
  permutations = 999,
  by = "margin"
)

# ggsave('Exploratory plots/nmds_inverts_plot.png', plot = inverts_nmds, height = 4, width = 5, dpi = 300, bg = "white")
#_____________________________________________________________________
# FUNCTIONAL DIVERSITY ----
#_____________________________________________________________________
library(FD)
library(SYNCSA)
library(vegan)
library(tidyverse)

#_____________________________________________________________________ 
# 1. Build FUNCTIONAL UNIT × trait matrix
#_____________________________________________________________________ 

trait_filt <- merged_traits %>%
  filter(!is.na(functional_unit)) %>%
  group_by(functional_unit) %>%
  summarise(
    log_body_size  = mean(log_body_size, na.rm = TRUE),
    microhabitat  = first(microhabitat),
    feeding_guild = first(feeding_guild),
    moisture_pref = first(moisture_pref),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(log_body_size),
    !is.na(microhabitat),
    !is.na(feeding_guild),
    !is.na(moisture_pref)
  ) %>%
  mutate(
    microhabitat  = factor(microhabitat),
    feeding_guild = factor(feeding_guild),
    moisture_pref = factor(moisture_pref)
  ) %>%
  column_to_rownames("functional_unit")
 
# # create trait dataframe for networks
# trait_filt_networks <- merged_traits %>%
#   filter(!is.na(functional_unit)) %>%
# 
#   group_by(functional_unit) %>%
#   summarise(
#     # functional traits
#     microhabitat   = first(microhabitat),
#     feeding_guild  = first(feeding_guild),
#     moisture_pref  = first(moisture_pref),
#     size_class     = first(size_class),
# 
#     # taxonomy
#     phylum  = first(phylum),
#     class   = first(class),
#     order   = first(order),
#     family  = first(family),
#     genus   = first(genus),
#     species = first(species),
# 
#     .groups = "drop"
#   ) %>%
# 
#   filter(
#     !is.na(microhabitat),
#     !is.na(feeding_guild),
#     !is.na(moisture_pref),
#     !is.na(size_class)
#   ) %>%
# 
#   mutate(
#     microhabitat  = factor(microhabitat),
#     feeding_guild = factor(feeding_guild),
#     moisture_pref = factor(moisture_pref),
#     size_class    = factor(size_class, ordered = TRUE)
#   ) %>%
# 
#   column_to_rownames("functional_unit")
# 
# write.csv(trait_filt_networks, "invert_traits_networks.csv", row.names = TRUE)
#_____________________________________________________________________ 
# 2. Build community matrix (sample × functional unit)
#_____________________________________________________________________ 

comm_filt <- merged_traits %>%
  dplyr::select(functional_unit, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "reads"
  ) %>%
  mutate(presence = reads > 0) %>%
  filter(!is.na(functional_unit)) %>%
  group_by(sample, functional_unit) %>%
  summarise(presence = as.integer(any(presence)), .groups = "drop") %>%
  pivot_wider(
    names_from  = functional_unit,
    values_from = presence,
    values_fill = 0
  ) %>%
  column_to_rownames("sample")

#_____________________________________________________________________ 
# 4. Functional Diversity metrics: RaoQ, FRed
#_____________________________________________________________________ 

rao_results <- rao.diversity(
  comm_filt,
  trait_filt,
  ord         = "metric",
  standardize = TRUE
)

rao_df <- data.frame(
  sample = rownames(comm_filt),
  RaoQ   = rao_results$FunRao,
  FRed   = rao_results$FunRedundancy
)

#_____________________________________________________________________
# 5. Functional Evenness (FEve)
#_____________________________________________________________________
# Remove ASVs absent from all samples (required for dbFD)
present_asvs <- colnames(comm_filt)[colSums(comm_filt) > 0]

comm_filt  <- comm_filt[, present_asvs, drop = FALSE]
trait_filt <- trait_filt[present_asvs, , drop = FALSE]

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

#_____________________________________________________________________
# 6. Combine metrics + metadata
#_____________________________________________________________________ 

meta_clean <- meta %>%
  rename_with(tolower) %>%
  mutate(sample = gsub("^x", "", sample))

diversity_metrics <- rao_df %>%
  left_join(feve_df,    by = "sample") %>%
  left_join(meta_clean, by = "sample") %>%
  distinct(sample, .keep_all = TRUE)

write.csv(diversity_metrics,
          "invert_18S_functional_metrics.csv",
          row.names = FALSE)

#_____________________________________________________________________
# PLOTTING FUNCTIONAL METRICS
#_____________________________________________________________________

diversity_metrics$site      <- factor(diversity_metrics$site)
diversity_metrics$tree.type <- factor(diversity_metrics$tree.type)

site_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(
  length(unique(diversity_metrics$site))
)

#_____________________________________________________________________
# FEve
#_____________________________________________________________________
p_feve <- ggplot(diversity_metrics, aes(x = site, y = FEve, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Functional Evenness (FEve)") +
  theme(legend.position = "none")

#_____________________________________________________________________
# RaoQ
#_____________________________________________________________________
p_rao <- ggplot(diversity_metrics, aes(x = site, y = RaoQ, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Rao's Quadratic Entropy") +
  theme(legend.position = "none")

#_____________________________________________________________________
# FRed
#_____________________________________________________________________
p_fred <- ggplot(diversity_metrics, aes(x = site, y = FRed, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Functional Redundancy") +
  theme(legend.position = "none")

p_feve; p_rao; p_fred

#_____________________________________________________________________
# Pine vs Birch
#_____________________________________________________________________

div_tree <- diversity_metrics %>% filter(tree.type %in% c("Birch", "Pine"))

tree_cols <- c("Birch" = "#c7d46b", "Pine" = "#4a9d63")

p_feve_tree <- ggplot(div_tree, aes(x = tree.type, y = FEve, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(y = "FEve", x = "Tree type") +
  theme(legend.position = "none")

p_rao_tree <- ggplot(div_tree, aes(x = tree.type, y = RaoQ, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(y = "RaoQ", x = "Tree type") +
  theme(legend.position = "none")

p_fred_tree <- ggplot(div_tree, aes(x = tree.type, y = FRed, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(y = "FRed", x = "Tree type") +
  theme(legend.position = "none")

p_feve_tree
p_rao_tree
p_fred_tree

#_____________________________________________________________________
# Woodland Condition
#_____________________________________________________________________

wca <- read.csv("wca.csv") %>% rename_with(tolower)

div_cond <- diversity_metrics %>%
  mutate(site = as.character(site)) %>%
  left_join(wca %>% mutate(site = as.character(site)), by = "site") %>%
  mutate(
    condition = factor(condition, levels = c(1,2,3),
                       labels = c("Low","Moderate","High"))
  )

cond_cols <- c("Low"="#d73027","Moderate"="#fee08b","High"="#1a9850")

p_feve_cond <- ggplot(div_cond, aes(x = condition, y = FEve, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  theme_bw(base_size = 12) +
  labs(y = "FEve", x = "Woodland condition") +
  theme(legend.position = "none")

p_rao_cond <- ggplot(div_cond, aes(x = condition, y = RaoQ, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  theme_bw(base_size = 12) +
  labs(y = "RaoQ", x = "Woodland condition") +
  theme(legend.position = "none")

p_fred_cond <- ggplot(div_cond, aes(x = condition, y = FRed, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  theme_bw(base_size = 12) +
  labs(y = "FRed", x = "Woodland condition") +
  theme(legend.position = "none")

p_feve_cond
p_rao_cond
p_fred_cond

#_____________________________________________________________________
# NDVI correlations
#_____________________________________________________________________
library(lme4)

div_ndvi <- diversity_metrics

site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(
  length(unique(div_ndvi$site))
)

plot_lmm_ndvi <- function(response, ylab, log_response = FALSE) {
  
  is_rao <- log_response
  
  if (!is_rao) {
    # ----------------------------
    # Standard Gaussian LMM
    # ----------------------------
    m <- lmer(
      as.formula(paste(response, "~ mean_ndvi + (1 | site)")),
      data = div_ndvi
    )
    
    newdat <- data.frame(
      mean_ndvi = seq(
        min(div_ndvi$mean_ndvi, na.rm = TRUE),
        max(div_ndvi$mean_ndvi, na.rm = TRUE),
        length.out = 100
      )
    )
    
    X <- model.matrix(~ mean_ndvi, newdat)
    beta <- fixef(m)
    V <- vcov(m)
    
    fit <- as.numeric(X %*% beta)
    se  <- sqrt(diag(X %*% V %*% t(X)))
    
    newdat$fit <- fit
    newdat$lwr <- fit - 1.96 * se
    newdat$upr <- fit + 1.96 * se
    
  } else {
    # ----------------------------
    # RaoQ: log + varPower (fungi-style)
    # ----------------------------
    m <- lme(
      fixed   = log1p(RaoQ) ~ mean_ndvi,
      random  = ~1 | site,
      weights = varPower(form = ~ fitted(.)),
      data    = div_ndvi
    )
    
    newdat <- data.frame(
      mean_ndvi = seq(
        min(div_ndvi$mean_ndvi, na.rm = TRUE),
        max(div_ndvi$mean_ndvi, na.rm = TRUE),
        length.out = 100
      )
    )
    
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
  
  # ----------------------------
  # Plot (raw scale)
  # ----------------------------
  p <- ggplot(div_ndvi, aes(x = mean_ndvi, y = .data[[response]])) +
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
      linewidth = 0.9
    ) +
    geom_point(aes(color = site), size = 2, alpha = 0.85) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(x = "Mean NDVI", y = ylab, color = "Site")
  
  list(plot = p, model = m)
}

invert_feve_ndvi <- plot_lmm_ndvi("FEve", "Functional Evenness (FEve)")
invert_rao_ndvi  <- plot_lmm_ndvi("RaoQ", "Rao’s Q", log_response = TRUE)
invert_fred_ndvi <- plot_lmm_ndvi("FRed", "Functional Redundancy (FRed)")

qqnorm(residuals(invert_rao_ndvi$model, type = "pearson"))
qqline(residuals(invert_rao_ndvi$model, type = "pearson"))

plot(fitted(invert_rao_ndvi$model),
     residuals(invert_rao_ndvi$model, type = "pearson"))
abline(h = 0, lty = 2)

# ggsave('Exploratory plots/raoq_ndvi_plot.png', plot = p_rao_ndvi$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/feve_ndvi_plot.png', plot = p_feve_ndvi$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/fred_ndvi_plot.png', plot = p_fred_ndvi$plot, height = 4, width = 5, dpi = 300, bg = "white")

#_____________________________________________________________________
# Age correlations
div_age <- diversity_metrics %>%
  mutate(log_age = log(x2024_age + 1))

site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(
  length(unique(div_ndvi$site))
)

plot_lmm_age <- function(response, ylab, log_response = FALSE) {
  
  is_rao <- log_response
  
  if (!is_rao) {
    m <- lmer(
      as.formula(paste(response, "~ log_age + (1 | site)")),
      data = div_age
    )
    
    newdat <- data.frame(
      log_age = seq(
        min(div_age$log_age, na.rm = TRUE),
        max(div_age$log_age, na.rm = TRUE),
        length.out = 100
      )
    )
    
    X <- model.matrix(~ log_age, newdat)
    beta <- fixef(m)
    V <- vcov(m)
    
    fit <- as.numeric(X %*% beta)
    se  <- sqrt(diag(X %*% V %*% t(X)))
    
    newdat$fit <- fit
    newdat$lwr <- fit - 1.96 * se
    newdat$upr <- fit + 1.96 * se
    
  } else {
    m <- lme(
      fixed   = log1p(RaoQ) ~ log_age,
      random  = ~1 | site,
      weights = varPower(form = ~ fitted(.)),
      data    = div_age
    )
    
    newdat <- data.frame(
      log_age = seq(
        min(div_age$log_age, na.rm = TRUE),
        max(div_age$log_age, na.rm = TRUE),
        length.out = 100
      )
    )
    
    emm <- emmeans(
      m,
      specs = ~ log_age,
      at    = list(log_age = newdat$log_age),
      type  = "response"
    )
    
    pred <- as.data.frame(emm)
    
    newdat$fit <- pred$response
    newdat$lwr <- pred$lower.CL
    newdat$upr <- pred$upper.CL
  }
  
  p <- ggplot(div_age, aes(x = log_age, y = .data[[response]])) +
    geom_ribbon(
      data = newdat,
      aes(x = log_age, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = newdat,
      aes(x = log_age, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.85,
      position = position_jitter(width = 0.03)
    ) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(x = "Site age (log)", y = ylab, color = "Site")
  
  list(plot = p, model = m)
}

invert_feve_age <- plot_lmm_age("FEve", "Functional Evenness (FEve)")
invert_rao_age  <- plot_lmm_age("RaoQ", "Rao’s Q", log_response = TRUE)
invert_fred_age <- plot_lmm_age("FRed", "Functional Redundancy (FRed)")

fixef(invert_rao_age$model)
anova(invert_rao_age$model)
summary(invert_rao_age$model)

fixef(invert_feve_age$model)
anova(invert_feve_age$model)
summary(invert_feve_age$model)
r2(invert_feve_age$model)

fixef(invert_fred_age$model)
anova(invert_fred_age$model)
summary(invert_fred_age$model)
r2(invert_fred_age$model)


## Combined diagnostic figure: Age models
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

## 1️⃣ RaoQ (lme)
res_rao_age <- residuals(invert_rao_age$model, type = "normalized")
fit_rao_age <- fitted(invert_rao_age$model, level = 0)

qqnorm(res_rao_age, main = "RaoQ ~ age: QQ plot")
qqline(res_rao_age)

plot(fit_rao_age, res_rao_age,
     xlab = "Fitted values",
     ylab = "Normalised residuals",
     main = "RaoQ ~ age: Residuals vs fitted")
abline(h = 0, lty = 2)

## 2️⃣ FEve (lmer)
res_feve <- residuals(invert_feve_age$model)
fit_feve <- fitted(invert_feve_age$model)

qqnorm(res_feve, main = "FEve ~ age: QQ plot")
qqline(res_feve)

plot(fit_feve, res_feve,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "FEve ~ age: Residuals vs fitted")
abline(h = 0, lty = 2)

## 3️⃣ FRed (lmer)
res_fred <- residuals(invert_fred_age$model)
fit_fred <- fitted(invert_fred_age$model)

qqnorm(res_fred, main = "FRed ~ age: QQ plot")
qqline(res_fred)

plot(fit_fred, res_fred,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "FRed ~ age: Residuals vs fitted")
abline(h = 0, lty = 2)

## Reset plotting
par(mfrow = c(1, 1))

#_____________________________________________________________________
# FULL PCoA + KDE Old (1–3) vs Young (7,8,10) comparison ----
#_____________________________________________________________________
library(cluster)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fastDummies)
library(RColorBrewer)
library(grid)
#_____________________________________________________________________
# 1. Remove outlier sample
#_____________________________________________________________________
outlier_sample <- "10_07"
comm_filt <- comm_filt[!(rownames(comm_filt) == outlier_sample), , drop = FALSE]

#_____________________________________________________________________
# 2. Align trait rows to community columns (CRITICAL)
#    Your species IDs are the ROW NAMES of trait_filt (ASV strings)
#_____________________________________________________________________
# Keep only species that occur in comm_filt
trait_filt2 <- trait_filt[rownames(trait_filt) %in% colnames(comm_filt), , drop = FALSE]

# Reorder to match community column order
trait_filt2 <- trait_filt2[match(colnames(comm_filt), rownames(trait_filt2)), , drop = FALSE]

# Sanity check
stopifnot(identical(rownames(trait_filt2), colnames(comm_filt)))

#_____________________________________________________________________
# 3. Prepare traits as ecological gradients (numeric-only for Gower)
#_____________________________________________________________________
trait_ord <- trait_filt2 %>%
  mutate(
    log_body_size = as.numeric(log_body_size),
    
    soil_depth = case_when(
      microhabitat == "ep" ~ -1,
      microhabitat == "he" ~  0,
      microhabitat == "eu" ~  1,
      TRUE ~ NA_real_
    ),
    
    moisture_grad = case_when(
      moisture_pref == "m" ~ -1,
      moisture_pref == "h" ~  1,
      TRUE ~ NA_real_
    ),
    
    is_aquatic     = as.numeric(microhabitat == "aq"),
    is_phytotelmic = as.numeric(microhabitat == "ph")
  ) %>%
  fastDummies::dummy_cols(
    select_columns = "feeding_guild",
    remove_first_dummy = FALSE,
    remove_selected_columns = TRUE
  ) %>%
  dplyr::select(where(is.numeric))

#_____________________________________________________________________
# 4. Gower distance + PCoA (species positions)
#_____________________________________________________________________
gower_dist <- daisy(trait_ord, metric = "gower")
pcoa_res   <- cmdscale(gower_dist, eig = TRUE, k = 2)

species_scores <- as.data.frame(pcoa_res$points)
colnames(species_scores) <- c("PC1", "PC2")

# Put IDs back
species_scores$ASV_id <- colnames(comm_filt)
rownames(species_scores) <- species_scores$ASV_id

# Axis labels
var_expl <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
pc1_lab <- paste0("PC1 (", var_expl[1], "%)")
pc2_lab <- paste0("PC2 (", var_expl[2], "%)")

#_____________________________________________________________________
# 5. Community weighted centroids (sample positions)
#_____________________________________________________________________
comm_mat <- comm_filt[, species_scores$ASV_id, drop = FALSE]
comm_mat <- as.matrix(comm_mat)

sample_scores <- comm_mat %*% as.matrix(species_scores[, c("PC1", "PC2")])
sample_scores <- as.data.frame(sample_scores)
colnames(sample_scores) <- c("PC1", "PC2")
sample_scores$sample <- rownames(sample_scores)

#_____________________________________________________________________
# 6. Add metadata
#_____________________________________________________________________
meta_clean <- meta %>%
  mutate(sample = gsub("^x", "", sample)) %>%
  dplyr::select(sample, site, tree.type, x2024_age)

sample_scores <- sample_scores %>%
  mutate(sample = gsub("^x", "", sample)) %>%
  left_join(meta_clean, by = "sample") %>%
  mutate(
    site = factor(site),
    tree.type = factor(tree.type),
    x2024_age = as.numeric(x2024_age)
  )

#_____________________________________________________________________
# 7. Trait arrows (correlation of traits vs PCoA axes)
#_____________________________________________________________________
trait_loadings <- 0.95 * cor(
  trait_ord,
  scale(species_scores[, c("PC1", "PC2")]),
  use = "pairwise.complete.obs"
)

trait_loadings_df <- as.data.frame(trait_loadings)
trait_loadings_df$trait <- rownames(trait_loadings_df)

trait_loadings_df <- trait_loadings_df %>%
  mutate(
    label_clean = case_when(
      trait == "soil_depth"     ~ "Soil depth",
      trait == "moisture_grad"  ~ "Moisture affinity",
      trait == "is_aquatic"     ~ "Aquatic",
      trait == "is_phytotelmic" ~ "Phytotelmic",
      TRUE ~ gsub("_", " ", gsub("^feeding_guild_", "Guild: ", trait))
    )
  )

plot_radius <- min(
  diff(range(sample_scores$PC1, na.rm = TRUE)),
  diff(range(sample_scores$PC2, na.rm = TRUE))
)

trait_loadings_df <- trait_loadings_df %>%
  mutate(
    PC1_scaled = PC1 * plot_radius,
    PC2_scaled = PC2 * plot_radius,
    strength   = sqrt(PC1^2 + PC2^2)
  ) %>%
  filter(strength > 0.3)

#_____________________________________________________________________
# 8. KDE old vs young
#_____________________________________________________________________
old_sites   <- c("1", "2", "3")
young_sites <- c("7", "8", "10")

old_scores   <- sample_scores %>% filter(site %in% old_sites)
young_scores <- sample_scores %>% filter(site %in% young_sites)

# Limits include arrows (so your plot matches the old one)
x_lim <- range(c(sample_scores$PC1, trait_loadings_df$PC1_scaled), na.rm = TRUE)
y_lim <- range(c(sample_scores$PC2, trait_loadings_df$PC2_scaled), na.rm = TRUE)

pad_x <- 0.04 * diff(x_lim)
pad_y <- 0.04 * diff(y_lim)
x_lim <- x_lim + c(-pad_x, pad_x)
y_lim <- y_lim + c(-pad_y, pad_y)

bw_x <- bandwidth.nrd(sample_scores$PC1)
bw_y <- bandwidth.nrd(sample_scores$PC2)

k_old <- kde2d(old_scores$PC1, old_scores$PC2,
               n = 200, lims = c(x_lim, y_lim), h = c(bw_x, bw_y))
k_young <- kde2d(young_scores$PC1, young_scores$PC2,
                 n = 200, lims = c(x_lim, y_lim), h = c(bw_x, bw_y))

k_diff <- k_old
k_diff$z <- k_old$z - k_young$z

kdf <- expand.grid(PC1 = k_diff$x, PC2 = k_diff$y)
kdf$z <- as.vector(k_diff$z)

#_____________________________________________________________________
# 9. Heatmap colours + breaks
#_____________________________________________________________________
max_abs <- max(abs(kdf$z), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 18)

cols_pal <- c(
  "#cf5300", "#de6408", "#e46e0d", "#ee8620",
  "#f4a247", "#f7b464", "#fbdcb4", "#fdeee0",
  "#ffffff", "#d0f2e4", "#9ad6c0", "#71bda3",
  "#5ebf9e", "#2f9f78", "#1b7f5e", "#0c6e4d", "#004d3a"
)

cols_pal <- c(
  "#c99600", "#d4a300", "#deaf16", "#f0c95a",
  "#f4d57a", "#fae8bd", "#fdf2d8",
  "#ffffff", "#e6f8f1", "#b3e2cf", "#9ad6c0", "#71bda3",
  "#3aa887", "#1f8668", "#0f7053", "#004d3a"
)

# Site colours
site_cols <- colorRampPalette(c("#000000", "#ededed"))(
  length(unique(sample_scores$site))
)
names(site_cols) <- sort(unique(sample_scores$site))

#_____________________________________________________________________
# 10. Final plot
#_____________________________________________________________________
p_pcoa_invert <- ggplot() +
  geom_contour_filled(
    data = kdf,
    aes(PC1, PC2, z = z),
    breaks = breaks
  ) +
  scale_fill_manual(values = cols_pal, name = "Old – Young\nDensity") +
  # geom_point(
  #   data = sample_scores,
  #   aes(PC1, PC2, color = site, shape = tree.type),
  #   size = 1.8, alpha = 0.8
  # ) +
  # scale_color_manual(values = site_cols, name = "Site") +
  # scale_shape_manual(
  #   values = c(
  #     "Birch" = 16,
  #     "Pine"  = 17,
  #     "None"  = 15
  #   ),
  #   name = "Tree type"
  # ) +
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
    colour = "black", linewidth = 0.55
  ) +
  # geom_text(
  #   data = trait_loadings_df,
  #   aes(x = PC1_scaled * 1.06, y = PC2_scaled * 1.06, label = label_clean),
  #   size = 3.2, vjust = -0.2
  # ) +
  coord_cartesian(
    xlim = x_lim,
    ylim = y_lim,
    clip = "off",
    expand = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank()) +
  labs(x = pc1_lab, y = pc2_lab)

p_pcoa_invert
# ggsave('Exploratory plots/invert_pcoa_plot.png', plot = p_pcoa_invert, height = 6, width = 8, dpi = 300, bg = "white")
#_____________________________________________________________________
# 11. Taxonomy join (now works, because species_scores has ASV_id)
#_____________________________________________________________________
species_pcoa_tax <- species_scores %>%
  left_join(
    merged_traits %>%
      dplyr::select(ASV_id, phylum, class, order, family) %>%
      distinct(),
    by = "ASV_id"
  )
#______________________________________________________________________________
# FUNCTIONAL RECOVERY ANALYSIS — Using Sites 1, 2, 3 as Reference Woodland
#______________________________________________________________________________
library(dplyr)
library(ggplot2)

sample_scores_clean <- sample_scores %>%
  mutate(
    x2024_age = as.numeric(x2024_age),
    site      = as.character(site)
  )

#----------------------------------------------------
# 2. Compute centroid of mature woodland (sites 1–3)
#----------------------------------------------------
ref_sites <- c("1", "2", "3")

centroid <- sample_scores_clean %>%
  dplyr::filter(site %in% ref_sites) %>%
  dplyr::summarise(
    centroid_PC1 = mean(PC1, na.rm = TRUE),
    centroid_PC2 = mean(PC2, na.rm = TRUE)
  )

#----------------------------------------------------
# 3. Distance of each sample to mature woodland centroid
#----------------------------------------------------
sample_scores_clean <- sample_scores_clean %>%
  mutate(
    dist_to_old = sqrt(
      (PC1 - centroid$centroid_PC1)^2 +
        (PC2 - centroid$centroid_PC2)^2
    )
  )

#----------------------------------------------------
# 4. Plot-only grouping: split ONLY the oldest age by site
#----------------------------------------------------
last_age <- max(sample_scores_clean$x2024_age, na.rm = TRUE)

sample_scores_clean <- sample_scores_clean %>%
  mutate(
    x2024_age = factor(
      x2024_age,
      levels = sort(unique(x2024_age))
    ),
    
    box_group = case_when(
      x2024_age == last_age & site == "2" ~ "left",
      x2024_age == last_age & site == "1" ~ "right",
      TRUE ~ "all"
    ),
    
    box_group = factor(box_group, levels = c("left", "right", "all"))
  )
#----------------------------------------------------
# 5. Functional recovery boxplot
#----------------------------------------------------
p_recovery_box <- ggplot(
  sample_scores_clean,
  aes(
    x = x2024_age,
    y = dist_to_old,
    fill = site,
    group = interaction(x2024_age, box_group)
  )
) +
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
  theme_bw(base_size = 14) +
  labs(
    x = "Woodland age (years)",
    y = "Functional distance to mature woodland"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

p_recovery_box
invert_recovery <- p_recovery_box
# ggsave('Exploratory plots/invert_functional_dist.png', plot = p_recovery_box, height = 4, width = 8, dpi = 300, bg = "white")

#____________________________________________________
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
    x_range = x_lim,
    y_range = y_lim
  )
}

# One-tailed p-value
p_kde <- (sum(perm_kde_stat >= obs_kde_stat) + 1) / (n_perm + 1)
p_kde
max(perm_kde_stat)
obs_kde_stat
summary(perm_kde_stat)

species_pcoa_tax <- species_scores %>%
  left_join(
    merged_traits %>%
      dplyr::select(ASV_id, phylum, class, order, family) %>%
      distinct(),
    by = "ASV_id"
  )

