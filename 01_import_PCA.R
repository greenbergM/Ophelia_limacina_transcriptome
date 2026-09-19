# =============================================================================
# 01_import_PCA.R
#
# Import Salmon quantifications, construct the DESeq2 dataset,
# perform variance-stabilizing transformation (VST), inspect sample
# relationships, remove the year/batch effect, and calculate mean
# expression for each developmental stage.
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(ggrepel)


# -----------------------------------------------------------------------------
# Transcript IDs and Salmon quantifications
# -----------------------------------------------------------------------------

gene2trans <- read.table(
  file = "uni_x20m30l25.txt"
)
# Transcript-level analysis: use transcript ID as both tx and gene identifier
gene2trans <- gene2trans[, c(1, 1)]


samples_names <- c(
  "OL21_01.out", "OL25_01.out", "OL25_02.out",

  "OL24_01.out", "OL24_02.out", "OL24_03.out",

  "OL24_06.out", "OL25_03.out", "OL25_04.out",

  "OL21_05.out", "OL21_06.out", "OL21_07.out",

  "OL24_19.out", "OL24_20.out", "OL24_09.out",

  "OL24_11.out", "OL24_12.out", "OL25_07.out",

  "OL21_09.out", "OL21_10.out", "OL21_11.out",

  "OL21_14.out", "OL21_15.out", "OL21_16.out",

  "OL24_13.out", "OL24_14.out", "OL24_15.out",

  "OL21_18.out", "OL21_19.out", "OL21_20.out"
)


files <- file.path(
  "salmon",
  samples_names,
  "filtered_quant.sf"
)


names(files) <- c(
  "egg_1", "egg_2", "egg_3",

  "8c_1", "8c_2", "8c_3",

  "16c_1", "16c_2", "16c_3",

  "32c_1", "32c_2", "32c_3",

  "72c_1", "72c_2", "72c_3",

  "78c_1", "78c_2", "78c_3",

  "gastr_1", "gastr_2", "gastr_3",

  "troch_1", "troch_2", "troch_3",

  "meta_1", "meta_2", "meta_3",

  "adult_1", "adult_2", "adult_3"
)


txi <- tximport(
  files,
  type = "salmon",
  tx2gene = gene2trans,
  countsFromAbundance = "lengthScaledTPM"
)


# -----------------------------------------------------------------------------
# Sample metadata
# -----------------------------------------------------------------------------

batches <- read.delim(
  "batches.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)


sampleTable <- data.frame(
  group = factor(
    batches$stage,
    levels = unique(batches$stage)
  ),
  year = factor(
    batches$year,
    levels = unique(batches$year)
  )
)

rownames(sampleTable) <- colnames(txi$counts)


# -----------------------------------------------------------------------------
# DESeq2 dataset and VST
# -----------------------------------------------------------------------------

dds <- DESeqDataSetFromTximport(
  txi,
  sampleTable,
  design = ~ year + group
)

vst_dds <- vst(dds)

# =============================================================================
# QC BEFORE BATCH CORRECTION
# =============================================================================


# -----------------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------------

pca_df <- plotPCA(
  vst_dds,
  intgroup = c("group", "year"),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_df, "percentVar")
)

pca_df$sample <- rownames(pca_df)

pca_uncorrected <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = group,
    shape = year
  )
) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    aes(label = sample),
    size = 3,
    show.legend = FALSE
  ) +
  xlab(
    paste0("PC1: ", percentVar[1], "%")
  ) +
  ylab(
    paste0("PC2: ", percentVar[2], "%")
  )

pca_uncorrected


# -----------------------------------------------------------------------------
# Sample distance heatmap
# -----------------------------------------------------------------------------

sampleDists <- dist(t(assay(vst_dds)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(
  sampleDistMatrix,
  annotation_col = sampleTable,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = viridis(100),
  show_colnames = FALSE,
  show_rownames = FALSE
)

# =============================================================================
# BATCH CORRECTION
# =============================================================================

assay(vst_dds) <- limma::removeBatchEffect(
  assay(vst_dds),
  batch = colData(vst_dds)$year,
  design = model.matrix(
    ~ group,
    data = as.data.frame(colData(vst_dds))
  )
)


vst_corrected <- assay(vst_dds)


# =============================================================================
# QC AFTER BATCH CORRECTION
# =============================================================================


# -----------------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------------

pcaData <- plotPCA(
  vst_dds,
  intgroup = "group",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)


# Display names for developmental stages
pcaData$group <- factor(
  pcaData$group,
  levels = c(
    "egg",
    "8c",
    "16c",
    "32c",
    "72c",
    "78c",
    "gastr",
    "troch",
    "meta",
    "adult"
  ),
  labels = c(
    "egg",
    "8 cells",
    "16 cells",
    "32 cells",
    "72 cells",
    "M-cells",
    "ciliated gastrula",
    "trochophore",
    "metatrochophore",
    "adult"
  )
)


stage_cols <- c(
  "egg" = "#E6867A",
  "8 cells" = "#D9A441",
  "16 cells" = "#B8B447",
  "32 cells" = "#72B36F",
  "72 cells" = "#3FAF8E",
  "M-cells" = "#46B7C3",
  "ciliated gastrula" = "#4B97D3",
  "trochophore" = "#6E7FDC",
  "metatrochophore" = "#9A73D9",
  "adult" = "#D474B3"
)


label_pos <- pcaData %>%
  group_by(group) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )


pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  geom_text_repel(
    data = label_pos,
    aes(PC1, PC2, label = group, color = group),
    size = 4,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 1.0,      
    point.padding = 1.2,    
    force = 3,              
    force_pull = 0.6,       
    segment.color = NA,
    segment.size = 0.5,
    min.segment.length = 0
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = stage_cols, name = "Developmental stage") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

pca_plot

ggsave(
  "./PCA_batch_corrected.png",
  pca_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Sample distance heatmap
# -----------------------------------------------------------------------------

sampleDists <- dist(
  t(assay(vst_dds))
)

sampleDistMatrix <- as.matrix(
  sampleDists
)


HMData <- as.data.frame(
  colData(vst_dds)[, c("group", "year"), drop = FALSE]
)


HMData$group <- factor(
  HMData$group,
  levels = c(
    "egg",
    "8c",
    "16c",
    "32c",
    "72c",
    "78c",
    "gastr",
    "troch",
    "meta",
    "adult"
  ),
  labels = c(
    "egg (0 h)",
    "8 cells (7 h)",
    "16 cells (8 h)",
    "32 cells (10 h)",
    "72 cells (15 h)",
    "M-cells (18 h)",
    "ciliated gastrula (24 h)",
    "trochophore (4 days)",
    "metatrochophore (1 month)",
    "adult (≥3 years)"
  )
)


HMData$year <- factor(
  HMData$year,
  levels = c(
    "2021",
    "2024",
    "2025"
  )
)


HMData <- HMData[
  colnames(sampleDistMatrix),
  ,
  drop = FALSE
]


colnames(HMData) <- c(
  "developmental stage",
  "batch(year)"
)


stage_cols_heatmap <- c(
  "egg (0 h)" = "#E6867A",
  "8 cells (7 h)" = "#D9A441",
  "16 cells (8 h)" = "#B8B447",
  "32 cells (10 h)" = "#72B36F",
  "72 cells (15 h)" = "#3FAF8E",
  "M-cells (18 h)" = "#46B7C3",
  "ciliated gastrula (24 h)" = "#4B97D3",
  "trochophore (4 days)" = "#6E7FDC",
  "metatrochophore (1 month)" = "#9A73D9",
  "adult (≥3 years)" = "#D474B3"
)


year_cols <- c(
  "2021" = "#1B9E77",
  "2024" = "#D95F02",
  "2025" = "#7570B3"
)


ann_colors <- list(
  "developmental stage" = stage_cols_heatmap,
  "batch(year)" = year_cols
)

HM_dist <- pheatmap(
  sampleDistMatrix,
  annotation_col = HMData,
  annotation_colors = ann_colors,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = magma(100),
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = FALSE,
  fontsize = 10,
  main = "",
  treeheight_row = 10,
  treeheight_col = 10
)

png(
  "./sample_distance_heatmap.png",
  width = 12 * 150,
  height = 10 * 150,
  res = 300
)

HM_dist

dev.off()

# =============================================================================
# MEAN EXPRESSION BY DEVELOPMENTAL STAGE
# =============================================================================

groups <- sub(
  "_[0-9]+$",
  "",
  colnames(vst_corrected)
)


vst_corrected_means <- sapply(
  unique(groups),
  function(g) {
    rowMeans(
      vst_corrected[, groups == g, drop = FALSE]
    )
  }
)

save(
  vst_corrected_means,
  file = "./vst_corrected_means.RData"
)