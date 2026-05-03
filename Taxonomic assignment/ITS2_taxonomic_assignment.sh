{\rtf1\ansi\ansicpg1252\cocoartf2868
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
# ============================================================\
# ITS2 QIIME 2 OTU Taxonomic Assignment\
# ============================================================\
\
# ----------------------------\
# 1. Create and activate QIIME 2 environment\
# ----------------------------\
\
CONDA_SUBDIR=osx-64 conda env create \\\
  --name qiime2-amplicon-2025.10 \\\
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-macos-latest-conda.yml\
\
conda activate qiime2-amplicon-2025.10\
conda config --env --set subdir osx-64\
\
qiime info\
\
\
# ----------------------------\
# 2. Import ITS2 OTU representative sequences\
# ----------------------------\
\
qiime tools import \\\
  --type 'FeatureData[Sequence]' \\\
  --input-path ~/Projects/metabarcoding/ITS2/data/ASVs_clustered97_UPPER.fasta \\\
  --output-path ~/Projects/metabarcoding/ITS2/qiime2/ITS2_rep_seqs.qza\
\
\
# ----------------------------\
# 3. Download UNITE reference data\
# ----------------------------\
\
qiime rescript get-unite-data \\\
  --p-version 2025-02-19 \\\
  --p-taxon-group fungi \\\
  --p-cluster-id 97 \\\
  --o-taxonomy unite-tax-97.qza \\\
  --o-sequences unite-seqs-97.qza\
\
\
# ----------------------------\
# 4. Train UNITE classifier\
# ----------------------------\
\
qiime feature-classifier fit-classifier-naive-bayes \\\
  --i-reference-reads unite-seqs-97.qza \\\
  --i-reference-taxonomy unite-tax-97.qza \\\
  --o-classifier unite-classifier-97.qza\
\
\
# ----------------------------\
# 5. Assign taxonomy\
# ----------------------------\
\
qiime feature-classifier classify-sklearn \\\
  --i-classifier unite-classifier-97.qza \\\
  --i-reads ~/Projects/metabarcoding/ITS2/qiime2/ITS2_rep_seqs.qza \\\
  --o-classification ~/Projects/metabarcoding/ITS2/qiime2/ITS2_taxonomy.qza\
\
\
# ----------------------------\
# 6. Export taxonomy\
# ----------------------------\
\
qiime tools export \\\
  --input-path ~/Projects/metabarcoding/ITS2/qiime2/ITS2_taxonomy.qza \\\
  --output-path ~/Projects/metabarcoding/ITS2/taxonomy/\
}