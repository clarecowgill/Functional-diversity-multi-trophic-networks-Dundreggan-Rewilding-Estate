#_____________________________________________________________________
# FUNGAL COMMUNITY DATA — CLEANING ----
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

data("fungi_func_FungalTraits", package = "microeco")

#_____________________________________________________________________
## 1. Load and clean data----
## remove occurances of 2 or fewer reads, and remove OTUs with 5 or fewer reads across all samples
#_____________________________________________________________________

its2 <- read.delim("OTU_table_ITS2.tsv", check.names = FALSE, sep = "\t")
meta <- read.csv("soil_meta.csv") |> rename_with(tolower)

# Fix error in sample naming
fix_last <- c("06"="14","07"="15","08"="16")
cn <- colnames(its2)
parts <- strsplit(cn,"-",fixed=TRUE)
first <- sapply(parts,`[`,1)
mid <- sapply(parts,`[`,2)
last <- sapply(parts,`[`,3)
suffix <- paste(mid,last,sep="-")
cn_new <- cn
safe_num <- suppressWarnings(as.numeric(first))
for(suf in unique(suffix)){
  idx <- which(suffix==suf)
  if(length(idx)>1){
    keep <- idx[which.min(safe_num[idx])]
    dupes <- setdiff(idx,keep)
    for(i in dupes){
      old_last <- last[i]
      if(old_last %in% names(fix_last)){
        new_last <- fix_last[[old_last]]
        cn_new[i] <- paste(first[i],mid[i],new_last,sep="-")
      }
    }
  }
}
colnames(its2) <- cn_new

# Identify numeric columns = samples
is_sample <- sapply(its2, is.numeric)
sample_cols <- names(its2)[is_sample]
tax_cols    <- names(its2)[!is_sample]

# Per sample threshold (2 reads)
its2[ , sample_cols] <- lapply(its2[ , sample_cols], function(x) {
  x[x <= 2] <- 0
  x
})

# Global threshold (5 reads)
global_totals <- rowSums(its2[ , sample_cols])
its2_filtered <- its2[global_totals >= 5, ]

# Sort sample col names
new_sample_names <- sapply(sample_cols, function(x) {
  parts <- strsplit(x, "-", fixed = TRUE)[[1]]
  parts <- parts[-1]                 # drop the first segment
  paste(parts, collapse = "_")       # join remaining with _
})

# Apply new names only to sample columns
colnames(its2_filtered)[colnames(its2_filtered) %in% sample_cols] <- new_sample_names

# Clean taxonomy names
prefix_pattern <- "[kpcofgs]__"

its2_filtered[ , tax_cols] <- lapply(its2_filtered[ , tax_cols], function(x) {
  gsub(prefix_pattern, "", x)
})

write.csv(its2_filtered, "ITS2_filtered_cont.csv", row.names = FALSE)

# Remove contaminants (consistently high reads in blanks with > family level assignment)
# ASVs 104, 1235
contaminants <- c("ASV_104", "ASV_1235")

its2_filtered <- its2_filtered[ !(its2_filtered$OTU_id %in% contaminants), ]

write.csv(its2_filtered, "ITS2_filtered.csv", row.names = FALSE)

# Check read depth per sample
# Identify numeric sample columns
is_sample <- sapply(its2_filtered, is.numeric)
sample_cols <- names(its2_filtered)[is_sample]

# Calculate total reads per sample
read_depth <- colSums(its2_filtered[, sample_cols])

# Convert to dataframe for plotting
read_depth_df <- data.frame(
  sample = names(read_depth),
  reads = as.numeric(read_depth)
)

head(read_depth_df)

# remove poor-read depth samples
low_depth_samples <- read_depth_df$sample[read_depth_df$reads < 5000]
low_depth_samples

its2_filtered <- its2_filtered[, !(colnames(its2_filtered) %in% low_depth_samples)]

# Remove '_<rank>_Incertae_sedis' suffixes from taxonomy columns
incertae_pattern <- "_[a-z]+_Incertae_sedis"

its2_filtered[, tax_cols] <- lapply(
  its2_filtered[, tax_cols],
  function(x) {
    str_remove(x, incertae_pattern)
  }
)
#_____________________________________________________________________
# FUNGAL COMMUNITY DATA — INITIAL EXPLORATION (genus-LEVEL)----
#_____________________________________________________________________

its2 <- its2_filtered

# Preview the trait database
colnames(fungi_func_FungalTraits)
head(fungi_func_FungalTraits)
# Rename GENUS → genus
fungi_func_FungalTraits <- dplyr::rename(
  fungi_func_FungalTraits,
  genus = GENUS
)
#_____________________________________________________________________
## TRAIT MERGE ----
#_____________________________________________________________________

library(dplyr)
library(stringr)

# Prepare clean FungalTraits reference
fungal_traits_ref <- fungi_func_FungalTraits %>%
  mutate(genus = as.character(str_trim(genus)))

# Merge safely (join only on genus)
merged_traits <- its2 %>%
  left_join(fungal_traits_ref, by = "genus")

# How many OTUs have a known lifestyle?
sum(!is.na(merged_traits$primary_lifestyle))
# Proportion of otus with traits
mean(!is.na(merged_traits$primary_lifestyle))
# Distribution of fungal trophic modes
table(merged_traits$primary_lifestyle, useNA = "ifany")

write.csv(merged_traits, "genus_pa_traits_pre.csv", row.names = FALSE)

# trait dataframe cleaned up and loaded back in
merged_traits <- read.csv("genus_pa_traits_pre.csv")

#______________________________________________________________________________
## Richness ----
#______________________________________________________________________________

# Load metadata and ensure sample names align
meta <- read.csv("soil_meta.csv")
colnames(meta) <- tolower(colnames(meta))  # ensure lowercase names
head(meta)

sample_cols <- merged_traits |>
  dplyr::select(where(is.numeric)) |>
  colnames()

# Prepare OTU long table (presence/absence), removing blanks ("EB" or "NG")
otu_long <- merged_traits %>%
  dplyr::select(OTU_id, genus, primary_lifestyle, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "reads"
  ) %>%
  filter(!grepl("EB|NG", sample, ignore.case = TRUE)) %>%
  mutate(
    sample   = gsub("^X", "", sample),
    presence = if_else(reads > 0, 1, 0)
  )

# Join metadata
otu_long <- left_join(otu_long, meta, by = "sample")

# Ensure factors
otu_long <- otu_long %>%
  mutate(
    site = factor(site),
    tree.type = factor(tree.type),
    primary_lifestyle = factor(primary_lifestyle)
    )

write.csv(otu_long, "fungi_genus_long.csv", row.names = TRUE)

# ---- OTU richness (presence/absence)
otu_richness <- otu_long %>%
  group_by(sample, site, tree.type) %>%
  summarise(otu_richness = sum(presence), .groups = "drop")

# ---- Functional group presence/absence per sample
func_presence <- otu_long %>%
  group_by(sample, site, tree.type, primary_lifestyle) %>%
  summarise(present = as.integer(any(presence > 0)), .groups = "drop")

# ---- Functional richness (# of guilds present per sample)
func_richness <- func_presence %>%
  group_by(sample, site, tree.type) %>%
  summarise(func_richness = sum(present), .groups = "drop")

# ---- Ensure categorical variables are factors
otu_richness <- otu_richness %>%
  mutate(site = factor(site),
         tree.type = factor(tree.type))

func_richness <- func_richness %>%
  mutate(site = factor(site),
         tree.type = factor(tree.type))

# --- Colour palette
site_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(otu_richness$site)))
tree_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(otu_richness$tree.type)))

# --- OTU richness per site
ggplot(otu_richness, aes(x = site, y = otu_richness, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "OTU richness") +
  theme(legend.position = "none")

# --- OTU richness per tree type
ggplot(otu_richness, aes(x = tree.type, y = otu_richness, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Tree type", y = "OTU richness") +
  theme(legend.position = "none")

# --- Functional group richness per site
ggplot(func_richness, aes(x = site, y = func_richness, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Functional group richness (guild count)") +
  theme(legend.position = "none")

# --- Functional group richness per tree type
ggplot(func_richness, aes(x = tree.type, y = func_richness, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Tree type", y = "Functional group richness (guild count)") +
  theme(legend.position = "none")

# Plot richness against age and NDVI ----
otu_richness_meta <- otu_richness %>%
  left_join(
    meta %>% dplyr::select(sample, site, x2024_age, mean_ndvi),
    by = "sample"
  ) %>%
  mutate(
    site = coalesce(as.character(site.x), as.character(site.y)),
    x2024_age = as.numeric(x2024_age),
    mean_ndvi = as.numeric(mean_ndvi)
  ) %>%
  dplyr::select(-site.x, -site.y)

# OTU accumulation curves ----
library(iNEXT)
library(RColorBrewer)
# Build PA matrix: samples x OTU
fungi_pa <- otu_long %>%
  dplyr::select(sample, OTU_id, presence) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from  = OTU_id,
    values_from = presence,
    values_fill = 0
  )
make_pa_matrix <- function(pa_matrix, samples) {
  pa_sub <- pa_matrix[rownames(pa_matrix) %in% samples, , drop = FALSE]
  pa_sub[pa_sub > 0] <- 1
  t(pa_sub)  # taxa × samples
}

fungi_pa_mat <- fungi_pa %>%
  tibble::column_to_rownames("sample") %>%
  as.matrix()

fungi_samples <- intersect(
  meta$sample,
  rownames(fungi_pa_mat)
)

fungi_matrix_all <- make_pa_matrix(
  fungi_pa_mat,
  fungi_samples
)

fungi_overall_out <- iNEXT(
  list("Soil fungi" = fungi_matrix_all),
  datatype = "incidence_raw",
  endpoint = 1250,
  conf = 0.95,
  se = TRUE
)

p_fungi_overall <- ggiNEXT(
  fungi_overall_out,
  type = 1,
  color.var = "Order.q"
) +
  scale_colour_manual(
    values = c("0" = "#bd4b04"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("0" = scales::alpha("#bd4b04", 0.25)),
    guide = "none"
  ) +
  labs(
    x = "Number of soil samples",
    y = "OTU richness",
    title = "Soil fungal OTU accumulation (all sites pooled)"
  ) +
  theme_minimal()

p_fungi_overall

ggsave('Exploratory plots/p_fungi_overall.png', plot = p_fungi_overall, height = 5, width = 6, dpi = 300, bg = "white")
#________
# OTU accumulation by site
fungi_meta <- meta %>%
  mutate(site = factor(site))

fungi_site_list <- split(fungi_meta$sample, fungi_meta$site)

fungi_site_mats <- lapply(fungi_site_list, function(samps) {
  samps <- intersect(samps, rownames(fungi_pa_mat))
  make_pa_matrix(fungi_pa_mat, samps)
})

fungi_site_out <- iNEXT(
  fungi_site_mats,
  datatype = "incidence_raw",
  endpoint = 150,
  conf = 0.95,
  se = TRUE
)

site_levels <- sort(as.numeric(
  unique(ggiNEXT(fungi_site_out)$data$Assemblage)
))

ylgn <- RColorBrewer::brewer.pal(9, "YlGn")

site_cols <- c(
  "#002e23",
  rev(colorRampPalette(ylgn[2:9])(length(site_levels) - 1))
)
names(site_cols) <- site_levels

site_fills <- scales::alpha(site_cols, 0.25)

p_fungi_site <- ggiNEXT(
  fungi_site_out,
  type = 1,
  color.var = "Assemblage"
)

p_fungi_site$data$Assemblage <- factor(
  p_fungi_site$data$Assemblage,
  levels = site_levels
)

p_fungi_site <- p_fungi_site +
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
    y = "OTU richness",
    title = "Soil fungal OTU accumulation by site"
  ) +
  theme_minimal()

p_fungi_site
ggsave('Exploratory plots/p_fungi_site.png', plot = p_fungi_site, height = 5, width = 6, dpi = 300, bg = "white")

# ---- Join woodland condition to OTU richness ----
wca <- read.csv("wca.csv") %>%
  rename_with(tolower) %>%
  mutate(site = as.character(site))

otu_richness_meta <- otu_richness_meta %>%
  mutate(site = as.character(site)) %>%
  left_join(wca, by = "site") %>%
  mutate(
    condition = factor(
      condition,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High")
    )
  )

# enforce numeric site ordering (THIS is the key)
site_levels <- as.character(sort(as.numeric(unique(otu_richness_meta$site))))
otu_richness_meta <- otu_richness_meta %>%
  mutate(site = factor(site, levels = site_levels))

library(RColorBrewer)

site_cols_rich <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(length(site_levels))
names(site_cols_rich) <- site_levels

library(lme4)
library(performance)
library(broom.mixed)

plot_lmm_ndvi_richness <- function(df) {
  
  df_ndvi <- df %>% filter(!is.na(mean_ndvi))
  
  m <- lmer(
    log1p(otu_richness) ~ mean_ndvi + (1 | site),
    data = df_ndvi
  )
  
  newdat <- data.frame(
    mean_ndvi = seq(
      min(df_ndvi$mean_ndvi, na.rm = TRUE),
      max(df_ndvi$mean_ndvi, na.rm = TRUE),
      length.out = 100
    )
  )
  
  X <- model.matrix(~ mean_ndvi, newdat)
  beta <- fixef(m)
  V <- vcov(m)
  
  newdat$fit_log <- as.numeric(X %*% beta)
  newdat$se_log  <- sqrt(diag(X %*% V %*% t(X)))
  
  newdat$lwr_log <- newdat$fit_log - 1.96 * newdat$se_log
  newdat$upr_log <- newdat$fit_log + 1.96 * newdat$se_log
  
  newdat$fit <- expm1(newdat$fit_log)
  newdat$lwr <- expm1(newdat$lwr_log)
  newdat$upr <- expm1(newdat$upr_log)
  
  p <- ggplot(df_ndvi, aes(mean_ndvi, otu_richness)) +
    geom_ribbon(
      data = newdat,
      aes(mean_ndvi, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = newdat,
      aes(mean_ndvi, fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_point(aes(color = site), size = 2, alpha = 0.8) +
    scale_color_manual(values = site_cols_rich) +
    theme_bw(base_size = 12) +
    labs(
      x = "Mean NDVI",
      y = "OTU richness",
      color = "Site"
    )
  
  list(
    plot    = p,
    model   = m,
    anova   = anova(update(m, . ~ . - mean_ndvi), m),
    r2      = performance::r2(m),
    tidy    = broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE),
    diag    = function() {
      par(mfrow = c(1, 2))
      plot(fitted(m), resid(m),
           xlab = "Fitted values",
           ylab = "Residuals",
           main = "Residuals vs fitted")
      abline(h = 0, lty = 2)
      
      qqnorm(resid(m), main = "Normal Q–Q")
      qqline(resid(m))
    }
  )
}

otu_ndvi <- plot_lmm_ndvi_richness(otu_richness_meta)
otu_ndvi$plot
otu_ndvi$anova
otu_ndvi$r2
otu_ndvi$diag()

#--------------------------------------------------
# LMM: OTU richness ~ site age (log) + (1|site)
#--------------------------------------------------

plot_lmm_age_richness <- function(df) {
  
  df_age <- df %>%
    mutate(log_age = log(x2024_age + 1)) %>%   # +1 avoids log(0)
    filter(!is.na(log_age))
  
  # Fit model
  m <- lmer(
    log1p(otu_richness) ~ log_age + (1 | site),
    data = df_age
  )
  
  # Prediction grid
  newdat <- data.frame(
    log_age = seq(
      min(df_age$log_age, na.rm = TRUE),
      max(df_age$log_age, na.rm = TRUE),
      length.out = 100
    )
  )
  
  X <- model.matrix(~ log_age, newdat)
  beta <- fixef(m)
  V <- vcov(m)
  
  newdat$fit_log <- as.numeric(X %*% beta)
  newdat$se_log  <- sqrt(diag(X %*% V %*% t(X)))
  
  newdat$lwr_log <- newdat$fit_log - 1.96 * newdat$se_log
  newdat$upr_log <- newdat$fit_log + 1.96 * newdat$se_log
  
  newdat$fit <- expm1(newdat$fit_log)
  newdat$lwr <- expm1(newdat$lwr_log)
  newdat$upr <- expm1(newdat$upr_log)
  
  # Plot
  p <- ggplot(df_age, aes(log_age, otu_richness)) +
    geom_ribbon(
      data = newdat,
      aes(log_age, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = newdat,
      aes(log_age, fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.8,
      position = position_jitter(width = 0.03, height = 0)
    ) +
    scale_color_manual(values = site_cols_rich) +
    theme_bw(base_size = 12) +
    labs(
      x = "Site age (log)",
      y = "OTU richness",
      color = "Site"
    )
  
  list(
    plot    = p,
    model   = m,
    anova   = anova(update(m, . ~ . - mean_ndvi), m),
    r2      = performance::r2(m),
    tidy    = broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE),
    diag    = function() {
      par(mfrow = c(1, 2))
      plot(fitted(m), resid(m),
           xlab = "Fitted values",
           ylab = "Residuals",
           main = "Residuals vs fitted")
      abline(h = 0, lty = 2)
      
      qqnorm(resid(m), main = "Normal Q–Q")
      qqline(resid(m))
    }
  )
}

# Run + view
otu_age <- plot_lmm_age_richness(otu_richness_meta)
otu_age$plot
otu_age$anova
otu_age$r2
otu_age$diag()

# ggsave('Exploratory plots/richness_age_fungi_plot.png', plot = otu_age$plot, height = 4, width = 5, dpi = 300, bg = "white")
# ggsave('Exploratory plots/richness_ndvi_fungi_plot.png', plot = otu_ndvi$plot, height = 4, width = 5, dpi = 300, bg = "white")

rich_age_fungi  <- otu_age$plot
rich_ndvi_fungi <- otu_ndvi$plot

# Model results
fixef(otu_age$model)
fixef(otu_ndvi$model)

library(lmerTest)
anova(otu_age$model)
anova(otu_ndvi$model)

summary(otu_age$model)
summary(otu_ndvi$model)

library(performance)
r2(otu_age$model)
r2(otu_ndvi$model)

library(RColorBrewer)
ylgn <- brewer.pal(9, "YlGn")
age_cols_9 <- rev(colorRampPalette(ylgn[2:9])(9))
age_cols <- c("#002e23", age_cols_9)

# --- OTU richness by woodland condition
p_rich_cond <- ggplot(
  otu_richness_meta,
  aes(x = condition, y = otu_richness)
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
  ) +
  theme(legend.position = "none")

fungi_wca_rich <- p_rich_cond

otu_richness_meta <- otu_richness_meta %>%
  mutate(
    condition = case_when(
      total < 26                    ~ "Poor",
      total >= 26 & total <= 35     ~ "Moderate",
      total > 35                    ~ "Good"
    ),
    condition = factor(
      condition,
      levels = c("Poor", "Moderate", "Good")
    )
  )

m_fungi_cond <- lmer(
  log1p(otu_richness) ~ condition + (1 | site),
  data = otu_richness_meta
)

m_fungi_cond <- lmer(
  log1p(otu_richness) ~ condition + (1 | site),
  data = otu_richness_meta
)

anova(m_fungi_cond)

#______________________________________________________________________________
# NDMS ----
#______________________________________________________________________________
# 1️⃣ Create sample-by-otu matrix (presence/absence)
otu_matrix <- otu_long %>%
  dplyr::select(sample, site, tree.type, genus, presence) %>%
  group_by(sample, site, tree.type, genus) %>%
  summarise(presence = as.integer(any(presence > 0)), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = presence, values_fill = 0)

otu_matrix <- otu_matrix[, colSums(otu_matrix != 0) > 0]

# 2️⃣ Extract metadata and numeric matrix
meta_nmds <- otu_matrix %>% dplyr::select(sample, site, tree.type)
otu_mat <- otu_matrix %>% dplyr::select(-site, -tree.type)
otu_mat <- as.data.frame(otu_mat)
rownames(otu_mat) <- otu_mat$sample
otu_mat <- otu_mat %>% dplyr::select(-sample)

# 3️⃣ Remove empty samples (no OTUus detected)
empty_samples <- rownames(otu_mat)[rowSums(otu_mat) == 0]
if(length(empty_samples) > 0) {
  message("Removing ", length(empty_samples), " empty samples: ",
          paste(empty_samples, collapse = ", "))
  otu_mat <- otu_mat[rowSums(otu_mat) > 0, ]
  meta_nmds <- meta_nmds %>% filter(!sample %in% empty_samples)
}

# 4️⃣ Run NMDS with Jaccard distance
set.seed(123)
nmds <- metaMDS(otu_mat, distance = "jaccard", k = 2,
                trymax = 100, autotransform = FALSE)

# Check stress (lower = better)
nmds$stress

# 5️⃣ Prepare NMDS scores + metadata for plotting
nmds_points <- as.data.frame(nmds$points) %>%
  rownames_to_column("sample") %>%
  left_join(meta_nmds, by = "sample") %>%
  mutate(
    # Convert NA tree types to explicit category
    tree.type = ifelse(is.na(tree.type), "Unknown", as.character(tree.type)),
    tree.type = factor(tree.type),
    site = factor(site)
  )

# 6️⃣ Define colour and shape palettes
# Colour palette for sites (Spectral)
site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(
  length(unique(na.omit(nmds_points$site)))
)
# Shape palette for 4 possible tree types
shape_vals <- c(
  "Birch"   = 16,  # filled circle
  "Pine"    = 17,  # filled triangle
  "None"    = 15,  # filled square
  "NA" = 1    # hollow circle
)

# 7️⃣ Plot NMDS
nmds_fungi <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = site, shape = tree.type)) +
  geom_point(size = 3, alpha = 0.95) +
  # stat_ellipse(
  #   aes(group = site, color = site),
  #   level = 0.95, linewidth = 0.8, na.rm = TRUE
  # ) +
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
dist_jaccard <- vegdist(otu_mat, method = "jaccard", binary = TRUE)

# Create metadata for PERMANOVA
meta_perm <- meta_nmds %>%
  filter(sample %in% rownames(otu_mat)) %>%
  mutate(
    site = factor(site),
    tree.type = factor(tree.type)
  )

# Set rownames to match distance matrix
meta_perm <- as.data.frame(meta_perm)
rownames(meta_perm) <- meta_perm$sample

adonis2(
  dist_jaccard ~ site + tree.type,
  data = meta_perm,
  permutations = 999,
  by = "margin"
)

# ggsave('Exploratory plots/nmds_fungi_plot.png', plot = nmds_fungi, height = 4, width = 5, dpi = 300, bg = "white")

### Splitting by functional group ----
run_functional_nmds_plot <- function(
    otu_long,
    functional_group,
    site_cols,
    distance = "jaccard",
    k = 2,
    trymax = 100,
    seed = 123,
    title_prefix = "NMDS (Jaccard)"
) {
  
  # 1. Filter to functional group
  otu_long_fg <- otu_long %>%
    dplyr::filter(primary_lifestyle == functional_group)
  
  if (nrow(otu_long_fg) == 0) {
    stop("No rows found for functional group: ", functional_group)
  }
  
  # 2. Build sample × genus PA matrix
  otu_matrix_fg <- otu_long_fg %>%
    dplyr::select(sample, genus, presence) %>%
    dplyr::group_by(sample, genus) %>%
    dplyr::summarise(
      presence = as.integer(any(presence > 0)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from  = genus,
      values_from = presence,
      values_fill = 0
    )
  
  otu_mat_fg <- otu_matrix_fg %>%
    tibble::column_to_rownames("sample") %>%
    as.data.frame()
  
  # Remove empty samples
  otu_mat_fg <- otu_mat_fg[rowSums(otu_mat_fg) > 0, , drop = FALSE]
  
  if (nrow(otu_mat_fg) < 3) {
    stop("Too few samples remain after filtering for NMDS.")
  }
  
  # 3. Run NMDS
  set.seed(seed)
  nmds_fg <- vegan::metaMDS(
    otu_mat_fg,
    distance = distance,
    k = k,
    trymax = trymax,
    autotransform = FALSE
  )
  
  # 4. Reattach metadata SAFELY
  meta_fg <- otu_long_fg %>%
    dplyr::select(sample, site, tree.type) %>%
    distinct()
  
  nmds_points_fg <- as.data.frame(nmds_fg$points) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(meta_fg, by = "sample") %>%
    dplyr::mutate(
      site = factor(site),
      tree.type = factor(tree.type, levels = c("Birch", "Pine", "None"))
    )
  
  # 5. Define filled shapes
  shape_vals <- c(
    "Birch" = 21,  # filled circle
    "Pine"  = 24,  # filled triangle
    "None"  = 22   # filled square
  )
 
   # 6. Build plot
  p <- ggplot(
    nmds_points_fg,
    aes(
      MDS1, MDS2,
      shape  = tree.type,
      fill   = site,
      colour = site
    )
  ) +
    geom_point(size = 3.2, alpha = 0.9, stroke = 0.6) +
    scale_fill_manual(values = site_cols) +
    scale_colour_manual(values = site_cols) +
    scale_shape_manual(values = shape_vals) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(title_prefix, " — ", functional_group),
      x = "NMDS1",
      y = "NMDS2",
      shape = "Tree type",
      fill = "Site",
      colour = "Site"
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  # 7. Return everything useful
  list(
    nmds = nmds_fg,
    points = nmds_points_fg,
    plot = p,
    stress = nmds_fg$stress,
    otu_matrix = otu_mat_fg,
    functional_group = functional_group
  )
}

ecm_res <- run_functional_nmds_plot(
  otu_long = otu_long,
  functional_group = "ectomycorrhizal",
  site_cols = site_cols
)
ecm_res$plot

sap_res <- run_functional_nmds_plot(
  otu_long,
  functional_group = "soil_saprotroph",
  site_cols = site_cols
)
sap_res$plot

path_res <- run_functional_nmds_plot(
  otu_long,
  functional_group = "plant_pathogen",
  site_cols = site_cols
)
path_res$plot

# richness by functional group ----
meta <- meta %>%
  dplyr::mutate(log_age = log(x2024_age + 1))

plot_fg_richness <- function(
    otu_long,
    functional_group,
    meta,
    predictor = c("x2024_age", "log_age", "mean_ndvi"),
    site_cols,
    xlab = NULL,
    ylab = "Genus richness",
    jitter_width = 0.05
) {
  
  predictor <- match.arg(predictor)
  
  # 1. Filter to functional group
  otu_fg <- otu_long %>%
    dplyr::filter(primary_lifestyle == functional_group)
  
  if (nrow(otu_fg) == 0) {
    stop("No rows found for functional group: ", functional_group)
  }
  
  # 2. Calculate richness per sample
  richness_fg <- otu_fg %>%
    dplyr::filter(presence > 0) %>%
    dplyr::distinct(sample, genus) %>%
    dplyr::count(sample, name = "richness")
  
  # 3. Join metadata
  dat_fg <- richness_fg %>%
    dplyr::left_join(
      meta %>% 
        dplyr::select(sample, site, dplyr::all_of(predictor)),
      by = "sample"
    ) %>%
    dplyr::filter(!is.na(.data[[predictor]])) %>%
    dplyr::mutate(site = factor(site))
  
  # 4. Fit simple model (for smooth / CI)
  form <- reformulate(predictor, response = "richness")
  m <- lm(form, data = dat_fg)
  
  newdat <- data.frame(
    seq(min(dat_fg[[predictor]], na.rm = TRUE),
        max(dat_fg[[predictor]], na.rm = TRUE),
        length.out = 100)
  )
  names(newdat) <- predictor
  
  preds <- cbind(
    newdat,
    as.data.frame(predict(m, newdat, interval = "confidence"))
  )
  
  # 5. Plot
  p <- ggplot(dat_fg, aes(x = .data[[predictor]], y = richness)) +
    geom_point(
      aes(colour = site),
      position = position_jitter(width = jitter_width),
      size = 2.2,
      alpha = 0.8
    ) +
    geom_ribbon(
      data = preds,
      aes(
        x = .data[[predictor]],
        ymin = lwr,
        ymax = upr
      ),
      inherit.aes = FALSE,
      alpha = 0.25,
      fill = "grey70"
    ) +
    geom_line(
      data = preds,
      aes(x = .data[[predictor]], y = fit),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    scale_colour_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(functional_group, " fungi"),
      x = ifelse(is.null(xlab), predictor, xlab),
      y = ylab,
      colour = "Site"
    ) +
    theme(
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
  
  list(
    data = dat_fg,
    model = m,
    plot = p,
    functional_group = functional_group,
    predictor = predictor
  )
}

p_ecm_logage <- plot_fg_richness(
  otu_long = otu_long,
  functional_group = "ectomycorrhizal",
  meta = meta,
  predictor = "log_age",
  site_cols = age_cols,
  xlab = "Site age (log)"
)
p_ecm_logage$plot
summary(p_ecm_logage$model)

p_ecm_ndvi <- plot_fg_richness(
  otu_long = otu_long,
  functional_group = "ectomycorrhizal",
  meta = meta,
  predictor = "mean_ndvi",
  site_cols = age_cols,
  xlab = "Mean NDVI"
)
p_ecm_ndvi$plot
summary(p_ecm_ndvi$model)

# woodland condition
ecm_richness_meta <- otu_long %>%
  filter(primary_lifestyle == "ectomycorrhizal", presence > 0) %>%
  distinct(sample, genus) %>%
  count(sample, name = "ecm_richness") %>%
  left_join(
    otu_richness_meta %>%
      dplyr::select(sample, site, condition, total),
    by = "sample"
  ) %>%
  filter(!is.na(condition))

m_ecm_cond <- lmer(
  log1p(ecm_richness) ~ condition + (1 | site),
  data = ecm_richness_meta
)
anova(m_ecm_cond)
performance::r2(m_ecm_cond)
broom.mixed::tidy(m_ecm_cond, effects = "fixed", conf.int = TRUE)

p_ecm_cond <- ggplot(
  ecm_richness_meta,
  aes(x = condition, y = ecm_richness)
) +
  geom_boxplot(
    fill = "grey85",
    colour = "grey30",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = site),
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  scale_colour_manual(values = age_cols) +
  theme_bw(base_size = 12) +
  labs(
    x = "Woodland condition",
    y = "ECM genus richness"
  ) +
  theme(
    legend.position = "none"
  )

p_ecm_cond
#_____________________________________________________________________
# FUNCTIONAL GROUP SPLIT — FUNGI ----
#_____________________________________________________________________
library(RColorBrewer)
library(scales)
library(ggplot2)
library(dplyr)

# Proportion of OTUs by primary lifestyle per sample
prop_lifestyle <- otu_long %>%
  filter(!is.na(primary_lifestyle)) %>%
  group_by(sample, site, tree.type, primary_lifestyle) %>%
  summarise(otu_count = sum(presence), .groups = "drop") %>%
  group_by(sample, site, tree.type) %>%
  mutate(prop = otu_count / sum(otu_count)) %>%
  ungroup()

# Pool across tree types within each site
site_pooled <- otu_long %>%
  filter(!is.na(primary_lifestyle)) %>%
  group_by(site, primary_lifestyle) %>%
  summarise(total_OTUs = sum(presence > 0), .groups = "drop") %>%
  group_by(site) %>%
  mutate(prop = total_OTUs / sum(total_OTUs)) %>%
  ungroup()

# Define order: older woodland lifestyles at bottom
lifestyle_order <- c(
  "phototroph_associated",
  "specialist_saprotroph",
  "animal_parasite",
  "plant_pathogen",
  "arbuscular_mycorrhizal",
  "mycoparasite",
  "unspecified_saprotroph",
  "foliar_endophyte",
  "lichen_parasite",
  "lichenized",
  "litter_saprotroph",
  "soil_saprotroph",
  "root_endophyte",
  "wood_saprotroph",
  "ectomycorrhizal"
)

# Apply factor order
site_pooled$primary_lifestyle <- factor(
  site_pooled$primary_lifestyle,
  levels = lifestyle_order
)

# Create colour palette (dark = ancient)
n_lifestyles <- length(levels(site_pooled$primary_lifestyle))
palette_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(23)

# Plot
ggplot(site_pooled, aes(x = site, y = prop, fill = primary_lifestyle)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = palette_cols, na.value = "grey80") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Site",
    y = "Proportion of OTUs",
    fill = "Primary lifestyle",
    title = "Fungal lifestyle composition per site (ancient woodland lifestyles shown in darker colours)"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

tree_pooled <- otu_long %>%
  filter(!is.na(primary_lifestyle)) %>%
  group_by(tree.type, primary_lifestyle) %>%
  summarise(total_OTUs = sum(presence > 0), .groups = "drop") %>%
  group_by(tree.type) %>%
  mutate(prop = total_OTUs / sum(total_OTUs)) %>%
  ungroup()

# Handle missing/blank tree types explicitly if needed
tree_pooled$tree.type[is.na(tree_pooled$tree.type)] <- "Unknown"

# Use same colour palette as above
ggplot(tree_pooled, aes(x = tree.type, y = prop, fill = primary_lifestyle)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = palette_cols, na.value = "grey80") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Tree type",
    y = "Proportion of OTUs",
    fill = "Primary lifestyle",
    title = "Fungal lifestyle composition per tree type (pooled across sites)"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

#_____________________________________________________________________
# FUNCTIONAL DIVERSITY — FUNGI (RaoQ, FRed, FEve) ----
#_____________________________________________________________________

library(vegan)
library(SYNCSA)
library(FD)
library(cluster)
library(tidyverse)

#_____________________________________________________________________
## Build Functional Trait Matrix ----
#_____________________________________________________________________

# # 1. Keep one record per genus with the two traits
# fungal_traits <- merged_traits %>%
#   dplyr::distinct(genus, primary_lifestyle, fruitbody_broad) %>%
#   dplyr::filter(!is.na(genus)) %>%
#   dplyr::mutate(
#     primary_lifestyle = as.factor(primary_lifestyle),
#     fruitbody_broad = as.factor(fruitbody_broad)
#   )
# 
# trait_matrix <- fungal_traits %>%
#   tibble::column_to_rownames("genus")

library(dplyr)
library(stringr)

library(dplyr)
library(stringr)

fungal_traits <- merged_traits %>%
  mutate(
    Fruitbody_type_template = Fruitbody_type_template %>%
      str_trim() %>%
      str_to_lower()
  ) %>%
  distinct(genus, primary_lifestyle, Fruitbody_type_template) %>%
  filter(!is.na(genus)) %>%
  mutate(
    
    # broad ecological lifestyle (UNCHANGED)
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
    ),
    
    # refined fruiting body / reproductive morphology
    fruitbody_broad = case_when(
      Fruitbody_type_template == "none" ~ "none",
      
      Fruitbody_type_template %in% c(
        "zoosporangium",
        "zygosporangium"
      ) ~ "asexual_micro",
      
      Fruitbody_type_template %in% c(
        "rust",
        "smut"
      ) ~ "plant_parasitic_spores",
      
      Fruitbody_type_template == "apothecium_(hymenium_on_surface)" ~ "asco_open",
      
      Fruitbody_type_template %in% c(
        "perithecium_(hymenium_hidden,_narrow_opening)",
        "cleistothecium_(closed,_spherical)"
      ) ~ "asco_closed",
      
      # exposed / ephemeral macroscopic fruiting bodies
      Fruitbody_type_template %in% c(
        "agaricoid",
        "clavarioid",
        "phalloid",
        "tremelloid",
        "cyphelloid"
      ) ~ "macro_epigeous",
      
      # persistent, substrate-attached or buffered macroscopic forms
      Fruitbody_type_template %in% c(
        "polyporoid",
        "corticioid",
        "gasteroid",
        "gasteroid-hypogeous"
      ) ~ "macro_persistent",
      
      Fruitbody_type_template == "other" ~ "other",
      
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    broad_lifestyle = factor(broad_lifestyle),
    fruitbody_broad = factor(fruitbody_broad)
  ) %>%
  dplyr::select(genus, broad_lifestyle, fruitbody_broad)

# final trait table used for functional analyses
fungal_traits_clean <- fungal_traits %>%
  filter(
    !is.na(broad_lifestyle),
    !is.na(fruitbody_broad),
    broad_lifestyle != "other",
    fruitbody_broad != "other"
  )

# build one for the networks which has other relvant columns
fungal_traits_networks <- merged_traits %>%
  distinct(genus, primary_lifestyle, Decay_substrate_template, Plant_pathogenic_capacity_template,
           Animal_biotrophic_capacity_template) %>%
  filter(!is.na(genus)) %>%
  mutate(
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
  mutate(
    broad_lifestyle = factor(broad_lifestyle),
  ) %>%
  dplyr::select(genus, broad_lifestyle, Decay_substrate_template, Plant_pathogenic_capacity_template,
                Animal_biotrophic_capacity_template)

# final trait table used for functional analyses
fungal_traits_networks <- fungal_traits_networks %>%
  filter(
    !is.na(broad_lifestyle),
    broad_lifestyle != "other",
  )

trait_matrix <- fungal_traits_clean %>%
  tibble::column_to_rownames("genus")

# 2. Build genus-level community matrix (presence/absence)
comm <- otu_long %>%
  dplyr::group_by(sample, genus) %>%
  dplyr::summarise(presence = as.integer(any(presence > 0)), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = genus,
    values_from = presence,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("sample") %>%
  as.data.frame()

# 3. Remove genera with NA in either trait
trait_matrix_clean <- trait_matrix %>%
  dplyr::filter(
    !is.na(broad_lifestyle),
    !is.na(fruitbody_broad)
  )

# 4. Match genera between comm and traits
shared_genera <- intersect(colnames(comm), rownames(trait_matrix_clean))

comm_filt  <- comm[, shared_genera, drop = FALSE]
trait_filt <- trait_matrix_clean[shared_genera, , drop = FALSE]

# 5. Remove samples with ≤ 2 genera
comm_filt <- comm_filt[rowSums(comm_filt) > 2, , drop = FALSE]

# 6. Re-align traits once more after sample filtering
final_genera <- colnames(comm_filt)
trait_filt   <- trait_filt[final_genera, , drop = FALSE]

# 7. Remove genera absent from all samples (belt-and-braces)
present_genera <- colnames(comm_filt)[colSums(comm_filt) > 0]

comm_filt  <- comm_filt[, present_genera, drop = FALSE]
trait_filt <- trait_filt[present_genera, , drop = FALSE]

#_____________________________________________________________________
## Rao’s Q and Functional Redundancy (FRed) ----
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
## Functional Evenness (FEve) ----
#_____________________________________________________________________

# dbFD handles mixed traits (factors + numeric) via Gower distance
feve_res <- dbFD(
  x         = trait_filt,   # traits (rows = genera)
  a         = comm_filt,    # community (rows = samples)
  calc.FRic = FALSE,
  calc.FDiv = FALSE,
  calc.CWM  = FALSE,
  stand.x   = FALSE,
  corr      = "cailliez"    # safe correction for negative eigenvalues
)

feve_df <- data.frame(
  sample = rownames(comm_filt),
  FEve   = feve_res$FEve
)

#_____________________________________________________________________
## Combine All Metrics + Metadata ----
#_____________________________________________________________________

# make sure meta is clean and sample IDs match
meta_clean <- meta %>%
  dplyr::rename_with(tolower) %>%
  dplyr::mutate(sample = gsub("^x", "", sample))

diversity_metrics <- rao_df %>%
  dplyr::left_join(feve_df,    by = "sample") %>%
  dplyr::left_join(meta_clean, by = "sample") %>%
  dplyr::distinct(sample, .keep_all = TRUE)

write.csv(diversity_metrics, "fungal_functional_metrics.csv", row.names = FALSE)

#_____________________________________________________________________
## Visual Checks by Site ----
#_____________________________________________________________________

diversity_metrics$site      <- as.factor(diversity_metrics$site)
diversity_metrics$tree.type <- as.factor(diversity_metrics$tree.type)

# Colour palette (same as before)
site_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(
  length(unique(diversity_metrics$site))
)

# --- Functional Evenness (FEve)
p_feve <- ggplot(diversity_metrics, aes(x = site, y = FEve, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Functional Evenness (FEve)") +
  theme(legend.position = "none")

# --- Rao’s Quadratic Entropy (RaoQ)
p_rao <- ggplot(diversity_metrics, aes(x = site, y = RaoQ, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Rao's Quadratic Entropy (RaoQ)") +
  theme(legend.position = "none")

# --- Functional Redundancy (FRed)
p_fred <- ggplot(diversity_metrics, aes(x = site, y = FRed, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_cols) +
  scale_color_manual(values = site_cols) +
  theme_bw(base_size = 12) +
  labs(x = "Site", y = "Functional Redundancy") +
  theme(legend.position = "none")

# Display
p_feve; p_rao; p_fred

#_____________________________________________________________________
## Pine vs Birch ----
#_____________________________________________________________________

div_tree <- diversity_metrics %>%
  dplyr::filter(tree.type %in% c("Birch", "Pine"))

# Define consistent colours for tree types
tree_cols <- c("Birch" = "#c7d46b", "Pine" = "#4a9d63")

# --- FEve by Tree Type ---
p_feve_tree <- ggplot(div_tree, aes(x = tree.type, y = FEve, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(y = "Functional Evenness (FEve)", x = "Tree type") +
  theme(legend.position = "none")

# --- RaoQ by Tree Type ---
p_rao_tree <- ggplot(div_tree, aes(x = tree.type, y = RaoQ, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(y = "Rao’s Quadratic Entropy", x = "Tree type") +
  theme(legend.position = "none")

# --- FRed by Tree Type ---
p_fred_tree <- ggplot(div_tree, aes(x = tree.type, y = FRed, fill = tree.type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = tree.type), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = tree_cols) +
  scale_color_manual(values = tree_cols) +
  theme_bw(base_size = 12) +
  labs(y = "Functional Redundancy", x = "Tree type") +
  theme(legend.position = "none")

# Display
p_feve_tree
p_rao_tree
p_fred_tree

#_____________________________________________________________________
## Plot by Woodland Condition ----
#_____________________________________________________________________

# Load woodland condition data
wca <- read.csv("wca.csv") %>%
  dplyr::rename_with(tolower)

div_cond <- diversity_metrics %>%
  dplyr::mutate(site = as.character(site)) %>%
  dplyr::left_join(wca %>% dplyr::mutate(site = as.character(site)), by = "site") %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High")
    )
  )

# Define consistent colour palette for woodland condition
cond_cols <- c("Low" = "#d73027", "Moderate" = "#fee08b", "High" = "#1a9850")

# --- FEve
p_feve_cond <- ggplot(div_cond, aes(x = condition, y = FEve, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  theme_bw(base_size = 12) +
  labs(y = "Functional Evenness (FEve)", x = "Woodland condition") +
  theme(legend.position = "none")

# --- RaoQ
p_rao_cond <- ggplot(div_cond, aes(x = condition, y = RaoQ, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  theme_bw(base_size = 12) +
  labs(y = "Rao’s Quadratic Entropy", x = "Woodland condition") +
  theme(legend.position = "none")

# --- FRed
p_fred_cond <- ggplot(div_cond, aes(x = condition, y = FRed, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  theme_bw(base_size = 12) +
  labs(y = "Functional Redundancy", x = "Woodland condition") +
  theme(legend.position = "none")

# Display
p_feve_cond
p_rao_cond
p_fred_cond

#_____________________________________________________________________
## Functional diversity ~ Age / NDVI (clean & explicit)
#_____________________________________________________________________

library(lme4)
library(nlme)
library(ggplot2)
library(dplyr)
library(performance)
library(broom.mixed)
library(RColorBrewer)
library(emmeans)

#---------------------------------------------------------------------
# 0. Setup
#---------------------------------------------------------------------

diversity_metrics <- diversity_metrics %>%
  mutate(site = factor(site))

site_cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(
  length(levels(diversity_metrics$site))
)
names(site_cols) <- levels(diversity_metrics$site)

#---------------------------------------------------------------------
# 1. Model fitting (ONE job only)
#---------------------------------------------------------------------

fit_fd_model <- function(data, response, predictor,
                         random = "site",
                         log_response = FALSE,
                         weighted = FALSE) {
  
  resp <- if (log_response) paste0("log1p(", response, ")") else response
  
  form_lmer <- as.formula(
    paste(resp, "~", predictor, "+ (1 |", random, ")")
  )
  
  if (!weighted) {
    return(lmer(form_lmer, data = data))
  }
  
  # weighted Gaussian model
  lme(
    fixed   = as.formula(paste(resp, "~", predictor)),
    random  = as.formula(paste("~ 1 |", random)),
    weights = varPower(form = ~ fitted(.)),
    data    = data
  )
}

#---------------------------------------------------------------------
# 2. Prediction
#---------------------------------------------------------------------

predict_fd <- function(model, predictor, newdata) {
  
  emm <- emmeans(
    model,
    specs = as.formula(paste("~", predictor)),
    at    = setNames(list(newdata[[predictor]]), predictor),
    type  = "response"
  )
  
  pred <- as.data.frame(emm)
  
  # emmeans uses different column names depending on model
  fit_col <- if ("response" %in% names(pred)) "response" else "emmean"
  
  tibble(
    x   = pred[[predictor]],
    fit = pred[[fit_col]],
    lwr = pred$lower.CL,
    upr = pred$upper.CL
  )
}

#---------------------------------------------------------------------
# 3. Plotting (RAW data + model overlay)
#---------------------------------------------------------------------

plot_fd <- function(data, response, predictor, ylab,
                    model, jitter_x = 0) {
  
  newdat <- tibble(
    !!predictor := seq(
      min(data[[predictor]], na.rm = TRUE),
      max(data[[predictor]], na.rm = TRUE),
      length.out = 100
    )
  )
  
  pred <- predict_fd(model, predictor, newdat)
  
  ggplot(data, aes(x = .data[[predictor]], y = .data[[response]])) +
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = lwr, ymax = upr),
      inherit.aes = FALSE,
      alpha = 0.25
    ) +
    geom_line(
      data = pred,
      aes(x = x, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = site),
      size = 2,
      alpha = 0.85,
      position = position_jitter(width = jitter_x, height = 0)
    ) +
    scale_color_manual(values = site_cols) +
    theme_bw(base_size = 12) +
    labs(x = predictor, y = ylab, color = "Site")
}

#---------------------------------------------------------------------
# 4. DATA PREP
#---------------------------------------------------------------------

dat_ndvi <- diversity_metrics %>%
  filter(!is.na(mean_ndvi))

dat_age <- diversity_metrics %>%
  mutate(log_age = log(x2024_age + 1)) %>%
  filter(!is.na(log_age))

#---------------------------------------------------------------------
# 5. MODEL CHOICES (based on diagnostics you already ran)
#---------------------------------------------------------------------

# FEve → untransformed LMM
m_feve_ndvi <- fit_fd_model(dat_ndvi, "FEve", "mean_ndvi")
m_feve_age  <- fit_fd_model(dat_age,  "FEve", "log_age")

# RaoQ → log + weighted Gaussian
m_rao_ndvi <- fit_fd_model(dat_ndvi, "RaoQ", "mean_ndvi",
                           log_response = TRUE, weighted = TRUE)
m_rao_age  <- fit_fd_model(dat_age,  "RaoQ", "log_age",
                           log_response = TRUE, weighted = TRUE)

# FRed → log + weighted Gaussian
m_fred_ndvi <- fit_fd_model(dat_ndvi, "FRed", "mean_ndvi",
                            log_response = TRUE, weighted = TRUE)
m_fred_age  <- fit_fd_model(dat_age,  "FRed", "log_age",
                            log_response = TRUE, weighted = TRUE)

qqnorm(residuals(m_fred_age))
qqline(residuals(m_fred_age))

plot(fitted(m_fred_age), resid(m_fred_age),
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2)

#---------------------------------------------------------------------
# 6. PLOTS (RAW SCALE, MODEL-BASED LINE)
#---------------------------------------------------------------------

fungi_feve_ndvi <- plot_fd(dat_ndvi, "FEve", "mean_ndvi",
                       "Functional Evenness (FEve)", m_feve_ndvi)

fungi_feve_age  <- plot_fd(dat_age, "FEve", "log_age",
                       "Functional Evenness (FEve)", m_feve_age, jitter_x = 0.03)

fungi_rao_ndvi <- plot_fd(dat_ndvi, "RaoQ", "mean_ndvi",
                      "Rao’s Quadratic Entropy", m_rao_ndvi)

fungi_rao_age  <- plot_fd(dat_age, "RaoQ", "log_age",
                      "Rao’s Quadratic Entropy", m_rao_age, jitter_x = 0.03)

fungi_fred_ndvi <- plot_fd(dat_ndvi, "FRed", "mean_ndvi",
                       "Functional Redundancy (FRed)", m_fred_ndvi)

fungi_fred_age  <- plot_fd(dat_age, "FRed", "log_age",
                       "Functional Redundancy (FRed)", m_fred_age, jitter_x = 0.03)

#---------------------------------------------------------------------
# 7. SAVE
#---------------------------------------------------------------------

# ggsave("Exploratory plots/feve_ndvi_plot.png",
#        plot = p_feve_ndvi, width = 5, height = 4, dpi = 300, bg = "white")
# 
# ggsave("Exploratory plots/feve_age_plot.png",
#        plot = p_feve_age, width = 5, height = 4, dpi = 300, bg = "white")
# 
# ggsave("Exploratory plots/raoq_ndvi_plot.png",
#        plot = p_rao_ndvi, width = 5, height = 4, dpi = 300, bg = "white")
# 
# ggsave("Exploratory plots/raoq_age_plot.png",
#        plot = p_rao_age, width = 5, height = 4, dpi = 300, bg = "white")
# 
# ggsave("Exploratory plots/fred_ndvi_plot.png",
#        plot = p_fred_ndvi, width = 5, height = 4, dpi = 300, bg = "white")
# 
# ggsave("Exploratory plots/fred_age_plot.png",
#        plot = p_fred_age, width = 5, height = 4, dpi = 300, bg = "white")

#===============================================================================
# FUNGAL TRAIT SPACE PCoA 
#===============================================================================
library(cluster)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fastDummies)
library(RColorBrewer)

#_____________________________________________________________________
# 1. Gower PCoA
#_____________________________________________________________________
gower_dist <- daisy(trait_filt, metric = "gower")
pcoa_res   <- cmdscale(gower_dist, eig = TRUE, k = 2)

genus_scores <- as.data.frame(pcoa_res$points)
colnames(genus_scores) <- c("PC1", "PC2")
genus_scores$genus <- rownames(genus_scores)

# Flip for consistent plotting
genus_scores$PC1 <- -genus_scores$PC1

var_expl <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
pc1_lab <- paste0("PC1 (", var_expl[1], "%)")
pc2_lab <- paste0("PC2 (", var_expl[2], "%)")

#_____________________________________________________________________
# 2. Sample positions (community-weighted)
comm_mat <- as.matrix(comm_filt[, rownames(trait_filt)])

sample_scores <- comm_mat %*% as.matrix(genus_scores[, c("PC1", "PC2")])
sample_scores <- as.data.frame(sample_scores)
colnames(sample_scores) <- c("PC1", "PC2")
sample_scores$sample <- rownames(sample_scores)

#_____________________________________________________________________
# 3. Join metadata
meta_clean <- meta %>%
  mutate(sample = gsub("^x", "", sample)) %>%
  dplyr::select(sample, site, x2024_age, tree.type)

sample_scores <- sample_scores %>%
  mutate(sample = gsub("^x", "", sample)) %>%
  left_join(meta_clean, by = "sample") %>%
  mutate(
    site = factor(site),
    tree.type = factor(tree.type),
    x2024_age = as.numeric(x2024_age)
  )

#_____________________________________________________________________
# 4. REMOVE OUTLIERS
outliers <- c("10_07", "09_11")

sample_scores_clean <- sample_scores %>%
  filter(!sample %in% outliers)

#_____________________________________________________________________
# 5. Trait loadings (dummy-coded)
trait_num <- trait_filt %>%
  mutate(across(everything(), as.character)) %>%
  fastDummies::dummy_cols(remove_first_dummy = FALSE) %>%
  dplyr::select(-broad_lifestyle, -fruitbody_broad)

trait_loadings <- cor(
  trait_num,
  scale(genus_scores[, c("PC1", "PC2")]),
  use = "pairwise.complete.obs"
)

trait_loadings_df <- as.data.frame(trait_loadings)
trait_loadings_df$trait <- rownames(trait_loadings_df)

trait_loadings_df <- trait_loadings_df %>%
  filter(sqrt(PC1^2 + PC2^2) > 0.3)
#_____________________________________________________________________
# 6. Scale trait arrows
plot_radius <- 0.75 * min(
  diff(range(sample_scores_clean$PC1)),
  diff(range(sample_scores_clean$PC2))
)

trait_loadings_df$PC1_scaled <- trait_loadings_df$PC1 * plot_radius
trait_loadings_df$PC2_scaled <- trait_loadings_df$PC2 * plot_radius

trait_loadings_df$label_clean <- trait_loadings_df$trait %>%
  gsub("broad_lifestyle_", "", .) %>%
  gsub("fruitbody_broad_", "", .) %>%
  gsub("_", " ", .)

#_____________________________________________________________________
# 7. KDE: Old (1–3) vs Young (7,8,10)
old_sites   <- c("1", "2", "3")
young_sites <- c("7", "8", "10")

old_scores   <- sample_scores_clean %>% filter(site %in% old_sites)
young_scores <- sample_scores_clean %>% filter(site %in% young_sites)

x_range <- range(c(sample_scores_clean$PC1, trait_loadings_df$PC1_scaled))
y_range <- range(c(sample_scores_clean$PC2, trait_loadings_df$PC2_scaled))

n_grid <- 200
bw_x <- bandwidth.nrd(sample_scores_clean$PC1)
bw_y <- bandwidth.nrd(sample_scores_clean$PC2)

k_old <- kde2d(
  old_scores$PC1, old_scores$PC2,
  n = n_grid, lims = c(x_range, y_range),
  h = c(bw_x, bw_y)
)

k_young <- kde2d(
  young_scores$PC1, young_scores$PC2,
  n = n_grid, lims = c(x_range, y_range),
  h = c(bw_x, bw_y)
)

k_diff <- k_old
k_diff$z <- k_old$z - k_young$z

kdf <- expand.grid(PC1 = k_diff$x, PC2 = k_diff$y)
kdf$z <- as.vector(k_diff$z)

#_____________________________________________________________________
# 8. Heatmap colours
breaks <- seq(-0.075, 0.055, length.out = 18)
zero_bin <- findInterval(0, breaks, rightmost.closed = TRUE)
n_bins <- length(breaks) - 1

max_abs <- max(abs(kdf$z), na.rm = TRUE)

cols_pal <- c(
  "#c99600", "#d4a300", "#deaf16", "#e8bc36", "#f0c95a",
  "#f4d57a", "#f7dfa0", "#fae8bd", "#fdf2d8",
  "#ffffff", "#e6f8f1", "#b3e2cf",
  "#6ec7a8", "#3aa887", "#1f8668", "#0f7053", "#004d3a"
)

#_____________________________________________________________________
# 9. Site colours
library(RColorBrewer)

site_cols <- colorRampPalette(c("#000000", "#ededed"))(
  length(unique(sample_scores_clean$site))
)
#_____________________________________________________________________
# 10. FINAL PLOT
x_lim <- range(kdf$PC1, na.rm = TRUE)
y_lim <- range(kdf$PC2, na.rm = TRUE)

p_pcoa_fungi <- ggplot() +
  geom_contour_filled(
    data = kdf,
    aes(PC1, PC2, z = z),
    breaks = breaks
  ) +
  scale_fill_manual(values = cols_pal, name = "Old – Young\nDensity") +
  # geom_point(
  #   data = sample_scores_clean,
  #   aes(PC1, PC2, color = site, shape = tree.type),
  #   size = 1.8, alpha = 0.8
  # ) +
  # scale_color_manual(values = site_cols, name = "Site") +
  # scale_shape_manual(
  #   values = c(
  #     Birch = 16,  # circle
  #     Pine  = 17,  # triangle
  #     None  = 15   # square
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
  #   aes(x = PC1_scaled * 1.08, y = PC2_scaled * 1.08, label = label_clean),
  #   size = 3.2, vjust = -0.2
  # ) +
  coord_cartesian(
    xlim = x_lim,
    ylim = y_lim,
    clip = "off",
    expand = FALSE
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

p_pcoa_fungi

# ggsave('Exploratory plots/fungi_pcoa_plot.png', plot = p_pcoa_fungi, height = 6, width = 8, dpi = 300, bg = "white")
#______________________________________________________________________________
# FUNCTIONAL RECOVERY ANALYSIS — Using Sites 1, 2, 3 as Reference Woodland
#______________________________________________________________________________
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
centroid <- sample_scores_clean %>%
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
sample_scores_clean <- sample_scores_clean %>%
  mutate(
    dist_to_old = sqrt(
      (PC1 - centroid$centroid_PC1)^2 +
        (PC2 - centroid$centroid_PC2)^2
    ),
    x2024_age = as.numeric(as.character(x2024_age))
  )
# Identify the oldest age
last_age <- max(sample_scores_clean$x2024_age, na.rm = TRUE)

# Plot-only grouping: split ONLY the last age by site
# and explicitly control left/right order
sample_scores_clean <- sample_scores_clean %>%
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
p_recovery_box <- ggplot(sample_scores_clean,
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
fungi_recovery <- p_recovery_box
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


# Plot genuses onto the PCoA for indicative genus silhouettes
library(ggrepel)
x_lims_sp <- range(genus_scores$PC1, na.rm = TRUE)
y_lims_sp <- range(genus_scores$PC2, na.rm = TRUE)

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
    data = genus_scores,
    aes(PC1, PC2),
    size = 2,
    colour = "black"
  ) +
  
  # repelled species labels with leader lines
  geom_text_repel(
    data = genus_scores,
    aes(PC1, PC2, label = genus),
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
