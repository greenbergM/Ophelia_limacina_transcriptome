# =============================================================================
# 03_transcription_factors.R
#
# Transcription-factor family composition and expression dynamics.
#
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)


# -----------------------------------------------------------------------------
# Transcription-factor PFAM annotations
# -----------------------------------------------------------------------------

# PFAM domains associated with transcription-factor families
pfam_interest <- read.table(
  "tf2pfam.tsv",
  header = TRUE,
  sep = "\t"
)

# Translated transcript to PFAM assignments
trans2pfam <- read.table(
  "transcript_to_tfpfam.tsv",
  header = FALSE,
  sep = "\t",
  col.names = c("Trans_ID", "PFAM")
)


# Split translated transcripts annotated with multiple PFAM domains
trans2pfam_long <- trans2pfam %>%
  separate_rows(
    PFAM,
    sep = ";"
  )


# Keep only PFAM domains associated with transcription factors
trans2pfam_filtered <- trans2pfam_long %>%
  inner_join(
    pfam_interest,
    by = c("PFAM" = "pfam")
  )


# =============================================================================
# TRANSCRIPTION-FACTOR FAMILY ABUNDANCE
# =============================================================================

fam_counts <- trans2pfam_filtered %>%
  distinct(
    Trans_ID,
    family
  ) %>%
  count(
    family,
    name = "n_trans"
  ) %>%
  mutate(
    percent = 100 * n_trans / sum(n_trans)
  ) %>%
  arrange(
    desc(percent)
  )


n_fam <- nrow(fam_counts)

my_cols <- scales::hue_pal()(
  n_fam
)


# -----------------------------------------------------------------------------
# TF-family abundance barplot
# -----------------------------------------------------------------------------

TF_barplot_linear <- ggplot(
  fam_counts,
  aes(
    x = reorder(family, percent),
    y = percent,
    fill = family
  )
) +
  geom_col(
    width = 0.8
  ) +
  geom_text(
    aes(label = n_trans),
    hjust = -0.2,
    size = 3
  ) +
  scale_fill_manual(
    values = my_cols
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  coord_flip() +
  labs(
    x = "TF family",
    y = "Percent of TF transcripts",
    title = "TF family abundance (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  ) +
  expand_limits(
    y = max(fam_counts$percent) * 1.1
  )


TF_barplot_linear


# -----------------------------------------------------------------------------
# Sqrt-scaled version used in the combined figure
# -----------------------------------------------------------------------------

TF_barplot <- ggplot(
  fam_counts,
  aes(
    x = reorder(family, percent),
    y = percent,
    fill = family
  )
) +
  geom_col(
    width = 0.8
  ) +
  geom_text(
    aes(label = n_trans),
    hjust = -0.2,
    size = 3
  ) +
  scale_fill_manual(
    values = my_cols
  ) +
  scale_y_sqrt(
    labels = function(x) paste0(round(x, 1), "%")
  ) +
  coord_flip() +
  labs(
    x = "TF family",
    y = "Percent of TF transcripts",
    title = "TF family abundance (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  ) +
  expand_limits(
    y = max(fam_counts$percent) * 1.15
  )


png(
  "barplot_tf.png",
  width = 16 * 150,
  height = 12 * 150,
  res = 300
)

TF_barplot

dev.off()


# =============================================================================
# TF-FAMILY EXPRESSION
# =============================================================================

load("./vst_corrected_means.RData")

# Convert stage-mean expression matrix to data frame
expr_df <- as.data.frame(
  vst_corrected_means
)

expr_df$Trans_ID <- rownames(
  expr_df
)


# Add TF-family annotation
expr_annot <- expr_df %>%
  inner_join(
    trans2pfam_filtered,
    by = "Trans_ID"
  )


# Convert expression data to long format
expr_long <- expr_annot %>%
  pivot_longer(
    cols = -c(
      Trans_ID,
      PFAM,
      family
    ),
    names_to = "stage",
    values_to = "expr"
  )


# 75th percentile of expression for each TF family and stage
family_summary <- expr_long %>%
  group_by(
    family,
    stage
  ) %>%
  summarise(
    q75 = quantile(
      expr,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# Heatmap matrix
# -----------------------------------------------------------------------------

heatmat <- family_summary %>%
  pivot_wider(
    names_from = stage,
    values_from = q75
  ) %>%
  column_to_rownames(
    "family"
  ) %>%
  as.matrix()


# Developmental stages included in the figure
stage_levels <- c(
  "egg",
  "8c",
  "16c",
  "32c",
  "72c",
  "78c",
  "gastr"
)

stage_names <- c(
  "egg",
  "8 cells",
  "16 cells",
  "32 cells",
  "72 cells",
  "M-cells",
  "ciliated gastrula"
)


heatmat <- heatmat[
  ,
  stage_levels,
  drop = FALSE
]


# Row-wise Z-score
heatmat_z <- t(
  scale(
    t(heatmat)
  )
)


# -----------------------------------------------------------------------------
# TF-family expression heatmap
# -----------------------------------------------------------------------------

TF_HM <- ComplexHeatmap::Heatmap(
  heatmat_z,
  name = "Z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  column_names_rot = 45,
  col = colorRampPalette(
    rev(
      brewer.pal(
        n = 11,
        name = "RdBu"
      )
    )
  )(100),
  row_title = NULL,
  column_title = "Family-level expression of TF",
  column_labels = stage_names
)


png(
  "TF_heatmap.png",
  width = 23 * 150,
  height = 15 * 150,
  res = 300
)

TF_HM

dev.off()

