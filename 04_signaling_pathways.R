# =============================================================================
# 04_signaling_pathways.R
#
# Abundance and developmental expression of major signaling pathways.
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)


# -----------------------------------------------------------------------------
# Annotation
# -----------------------------------------------------------------------------

annotation <- read.csv(
  "eggnog.txt",
  sep = "\t"
)


# -----------------------------------------------------------------------------
# Developmental signaling pathway definitions
# -----------------------------------------------------------------------------

dev_pathways <- read.table(
  "ko2pathway.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  rename(
    pathway_id = ko,
    pathway_name = pathway
  )


# =============================================================================
# TRANSCRIPT-TO-PATHWAY MAPPING
# =============================================================================

transcript_col <- "query"
pathway_col <- "KEGG_Pathway"


trans2pathway <- annotation %>%
  dplyr::select(
    Trans_ID = all_of(transcript_col),
    pathway_raw = all_of(pathway_col)
  ) %>%
  filter(
    !is.na(pathway_raw),
    pathway_raw != "",
    pathway_raw != "-"
  ) %>%
  separate_rows(
    pathway_raw,
    sep = "[,;]"
  ) %>%
  mutate(
    pathway_id = str_trim(pathway_raw),
    pathway_id = ifelse(
      str_detect(pathway_id, "^\\d+$"),
      paste0("ko", pathway_id),
      pathway_id
    )
  ) %>%
  distinct(
    Trans_ID,
    pathway_id
  )


# Keep only pathways included in the manually defined pathway table

trans2majorpathway <- trans2pathway %>%
  inner_join(
    dev_pathways,
    by = "pathway_id"
  )


# =============================================================================
# DEVELOPMENTAL SIGNALING PATHWAY ABUNDANCE
# =============================================================================

dev_counts <- trans2majorpathway %>%
  distinct(
    Trans_ID,
    pathway_name
  ) %>%
  count(
    pathway_name,
    name = "n_trans"
  ) %>%
  mutate(
    percent = 100 * n_trans / sum(n_trans)
  ) %>%
  arrange(
    desc(percent)
  )


my_cols <- scales::hue_pal()(
  nrow(dev_counts)
)


# Display labels
custom_path_names <- c(
  "Wnt" = "Wnt",
  "TGF-b" = "TGF-β",
  "Ras/MAPK" = "RAS/MAPK",
  "Hippo" = "Hippo",
  "Notch" = "Notch",
  "Hedgehog" = "Hedgehog"
)


pathway_barplot <- ggplot(
  dev_counts,
  aes(
    x = reorder(pathway_name, percent),
    y = percent,
    fill = pathway_name
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
  scale_x_discrete(
    labels = custom_path_names
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%")
  ) +
  coord_flip() +
  labs(
    x = "Signaling pathway",
    y = "Percent of pathway-annotated transcripts",
    title = "Pathway transcripts abundance (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  ) +
  expand_limits(
    y = max(dev_counts$percent) * 1.12
  )


ggsave(
  "./dev_pathway_barplot.png",
  pathway_barplot,
  width = 10,
  height = 7,
  dpi = 300
)

# =============================================================================
# MAJOR PATHWAY EXPRESSION DURING EMBRYONIC DEVELOPMENT
# =============================================================================

load("vst_corrected_means.RData")

expr_df <- as.data.frame(
  vst_corrected_means
)

expr_df$Trans_ID <- rownames(
  expr_df
)


expr_path_annot <- expr_df %>%
  inner_join(
    trans2majorpathway,
    by = "Trans_ID"
  )


pathway_summary <- expr_path_annot %>%
  pivot_longer(
    cols = -c(
      Trans_ID,
      pathway_id,
      pathway_name
    ),
    names_to = "stage",
    values_to = "expr"
  ) %>%
  group_by(
    pathway_name,
    stage
  ) %>%
  summarise(
    median_expr = median(
      expr,
      na.rm = TRUE
    ),
    n_trans = n_distinct(Trans_ID),
    .groups = "drop"
  )


row_labels_df <- pathway_summary %>%
  distinct(
    pathway_name,
    n_trans
  ) %>%
  mutate(
    label = pathway_name
  )

stage_levels_pre_troch <- c(
  "egg",
  "8c",
  "16c",
  "32c",
  "72c",
  "78c",
  "gastr"
)


stage_names_pre_troch <- c(
  "egg",
  "8 cells",
  "16 cells",
  "32 cells",
  "72 cells",
  "M-cells",
  "ciliated gastrula"
)


heatmat <- pathway_summary %>%
  left_join(
    row_labels_df,
    by = "pathway_name"
  ) %>%
  dplyr::select(
    label,
    stage,
    median_expr
  ) %>%
  pivot_wider(
    names_from = stage,
    values_from = median_expr
  ) %>%
  column_to_rownames(
    "label"
  ) %>%
  as.matrix()


heatmat <- heatmat[
  ,
  stage_levels_pre_troch,
  drop = FALSE
]


heatmat_z <- t(
  scale(
    t(heatmat)
  )
)

heatmat_z[
  is.na(heatmat_z)
] <- 0


pathway_HM <- ComplexHeatmap::Heatmap(
  heatmat_z,
  name = "Z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = custom_path_names[
    rownames(heatmat_z)
  ],
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
  column_title = "Pathway median expression",
  column_labels = stage_names_pre_troch
)

png(
  "pathway_heatmap.png",
  width = 23 * 150,
  height = 15 * 150,
  res = 300
)

pathway_HM

dev.off()

