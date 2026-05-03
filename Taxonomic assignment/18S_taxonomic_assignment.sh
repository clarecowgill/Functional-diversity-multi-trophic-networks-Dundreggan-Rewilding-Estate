{\rtf1\ansi\ansicpg1252\cocoartf2868
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 CourierNewPSMT;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs24 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 # ============================================================\
# 18S QIIME 2 ASV Taxonomic Assignment\
# ============================================================\
\
# ----------------------------\
# 1. Activate QIIME 2 environment\
# ----------------------------\
\
conda activate qiime2-amplicon-2025.10\
conda config --env --set subdir osx-64\
\
\
# ----------------------------\
# 2. Generate FASTA from DADA2 ASV table\
# ----------------------------\
\
cd ~/Projects/metabarcoding/18S/data\
\
Rscript - <<'EOF'\
library(tidyverse)\
library(Biostrings)\
\
seqtab <- readRDS("06_seqtab.nochim.rds")\
\
asv_df <- as.data.frame(t(seqtab))\
asv_seqs <- rownames(asv_df)\
\
dna <- DNAStringSet(asv_seqs)\
names(dna) <- asv_seqs\
\
writeXStringSet(dna, "ASV_seqs_for_qiime.fasta")\
EOF\
\
\
# ----------------------------\
# 3. Import ASV sequences into QIIME 2\
# ----------------------------\
\
qiime tools import \\\
--type 'FeatureData[Sequence]' \\\
--input-path ASV_seqs_for_qiime.fasta \\\
--output-path ASV_seqs.qza\
\
\
# ----------------------------\
# 4. Export ASV table from R\
# ----------------------------\
\
Rscript - <<'EOF'\
library(tidyverse)\
\
seqtab <- readRDS("06_seqtab.nochim.rds")\
\
asv_df <- as.data.frame(t(seqtab))\
asv_df$ASV_id <- rownames(asv_df)\
asv_df <- asv_df %>% relocate(ASV_id)\
\
write.table(\
asv_df,\
"ASV_table_for_biom.tsv",\
sep = "\\t",\
quote = FALSE,\
row.names = FALSE\
)\
EOF\
\
\
# ----------------------------\
# 5. Convert ASV table to BIOM\
# ----------------------------\
\
biom convert \\\
-i ASV_table_for_biom.tsv \\\
-o ASV_table_full.biom \\\
--table-type="OTU table" \\\
--to-hdf5\
\
\
# ----------------------------\
# 6. Import BIOM table into QIIME 2\
# ----------------------------\
\
qiime tools import \\\
--input-path ASV_table_full.biom \\\
--type 'FeatureTable[Frequency]' \\\
--input-format BIOMV210Format \\\
--output-path ASV_table.qza\
\
\
# ----------------------------\
# 7. Download and prepare SILVA reference database\
# ----------------------------\
\
cd ~/Projects/metabarcoding/18S/silva\
\
qiime rescript get-silva-data \\\
--p-version 138.1 \\\
--p-target SSURef_NR99 \\\
--p-include-species-labels \\\
--o-silva-sequences silva-ssu.qza \\\
--o-silva-taxonomy silva-taxonomy.qza\
\
qiime rescript reverse-transcribe \\\
--i-rna-sequences silva-ssu.qza \\\
--o-dna-sequences silva-ssu-DNA.qza\
\
qiime feature-classifier extract-reads \\\
--i-sequences silva-ssu-DNA.qza \\\
--p-f-primer CCAGCASCYGCGGTAATTCC \\\
--p-r-primer ACTTTCGTTCTTGATYRATGA \\\
--p-min-length 100 \\\
--p-max-length 500 \\\
--o-reads silva-ssu-extracted.qza\
\
qiime rescript cull-seqs \\\
--i-sequences silva-ssu-extracted.qza \\\
--p-num-degenerates 5 \\\
--o-clean-sequences silva-ssu-extracted-clean.qza\
\
qiime feature-classifier fit-classifier-naive-bayes \\\
--i-reference-reads silva-ssu-extracted-clean.qza \\\
--i-reference-taxonomy silva-taxonomy.qza \\\
--o-classifier silva-18S-V4-classifier.qza\
\
\
# ----------------------------\
# 8. Assign taxonomy\
# ----------------------------\
\
cd ~/Projects/metabarcoding/18S/data\
\
qiime feature-classifier classify-sklearn \\\
--i-classifier ../silva/silva-18S-V4-classifier.qza \\\
--i-reads ASV_seqs.qza \\\
--o-classification ASV_taxonomy.qza\
\
\
# ----------------------------\
# 9. Export taxonomy and feature table\
# ----------------------------\
\
qiime tools export \\\
--input-path ASV_taxonomy.qza \\\
--output-path taxonomy\
\
qiime tools export \\\
--input-path ASV_table.qza \\\
--output-path feature_table\
\
\
# ----------------------------\
# 10. Merge taxonomy with ASV table in R\
# ----------------------------\
\
Rscript - <<'EOF'\
library(tidyverse)\
\
seqtab <- readRDS("06_seqtab.nochim.rds")\
\
asv_df <- as.data.frame(t(seqtab))\
asv_df$ASV_id <- rownames(asv_df)\
asv_df <- asv_df %>% relocate(ASV_id)\
\
tax <- read_tsv("taxonomy/taxonomy.tsv",\
show_col_types = FALSE) %>%\
rename(\
ASV_id = `Feature ID`,\
taxonomy = Taxon,\
confidence = Confidence\
)\
\
max_ranks <- tax$taxonomy %>%\
strsplit(";") %>%\
map_int(length) %>%\
max()\
\
rank_names <- paste0("rank_", seq_len(max_ranks))\
\
tax_split <- tax %>%\
mutate(split = strsplit(taxonomy, ";")) %>%\
mutate(split = map(split, ~ \{length(.x) <- max_ranks; .x\})) %>%\
unnest_wider(split, names_sep = "") %>%\
rename_with(~ rank_names[seq_along(.)], starts_with("split"))\
\
merged <- asv_df %>% left_join(tax_split, by = "ASV_id")\
\
write_tsv(merged, "18S_ASV_table_with_taxonomy.tsv")\
\
write_tsv(\
merged %>% select(ASV_id, taxonomy, starts_with("rank_"), confidence),\
"18S_taxonomy_table.tsv"\
)\
EOF}