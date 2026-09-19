# =============================================================================
# 02_GO_analysis.R
#
# GO-slim classification of transcript annotations.
# GO terms assigned by eggNOG are mapped to the generic GO-slim categories
# for Biological Process (BP), Molecular Function (MF), and
# Cellular Component (CC).
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)


# -----------------------------------------------------------------------------
# Input data
# -----------------------------------------------------------------------------

annotation <- read.csv(
  "eggnog.txt",
  sep = "\t"
)

obo_path <- "go_analysis/goslim_generic.obo"


# =============================================================================
# GO-SLIM TERMS
# =============================================================================

# Extract GO IDs included in the generic GO-slim ontology

obo <- readLines(
  obo_path,
  warn = FALSE
)

goslim_ids <- unique(
  sub(
    "^id: ",
    "",
    grep("^id: GO:", obo, value = TRUE)
  )
)


# =============================================================================
# TRANSCRIPT-TO-GO MAPPING
# =============================================================================

gene_go <- annotation %>%
  transmute(
    gene = query,
    GOs
  ) %>%
  filter(
    !is.na(GOs),
    GOs != "-",
    GOs != ""
  ) %>%
  mutate(
    GOs = gsub("\\s+", "", GOs)
  ) %>%
  separate_rows(
    GOs,
    sep = ","
  ) %>%
  rename(
    go = GOs
  ) %>%
  filter(
    grepl("^GO:\\d{7}$", go)
  ) %>%
  distinct(
    gene,
    go
  )


# =============================================================================
# REMOVE INVALID AND OBSOLETE GO TERMS
# =============================================================================

gene_go <- gene_go %>%
  mutate(
    go = trimws(
      gsub("\r", "", go)
    )
  ) %>%
  distinct(
    gene,
    go
  )


# Check which GO IDs are present in GO.db

uniq_go <- unique(gene_go$go)

exists_tbl <- tibble(
  go = uniq_go,
  exists = vapply(
    uniq_go,
    function(x) !is.null(GOTERM[[x]]),
    logical(1)
  )
)


gene_go2 <- gene_go %>%
  inner_join(
    exists_tbl %>%
      filter(exists),
    by = "go"
  ) %>%
  dplyr::select(
    gene,
    go
  )


# Remove obsolete GO terms

term_tbl <- tibble(
  go = unique(gene_go2$go),
  term = Term(
    GOTERM[unique(gene_go2$go)]
  )
)


gene_go2 <- gene_go2 %>%
  inner_join(
    term_tbl,
    by = "go"
  ) %>%
  filter(
    !is.na(term),
    !grepl("^obsolete", term)
  ) %>%
  dplyr::select(
    gene,
    go
  ) %>%
  distinct()


# =============================================================================
# MAP GO TERMS TO GO-SLIM CATEGORIES
# =============================================================================

build_gene_slim <- function(
  gene_go_df,
  slim_ids,
  ont = c("BP", "MF", "CC")
) {

  ont <- match.arg(ont)

  anc_map <- switch(
    ont,
    BP = GOBPANCESTOR,
    MF = GOMFANCESTOR,
    CC = GOCCANCESTOR
  )


  # Keep GO terms belonging to the selected ontology

  gene_go_df <- gene_go_df %>%
    mutate(
      ont = Ontology(GOTERM[go])
    ) %>%
    filter(
      ont == !!ont
    )


  go_unique <- unique(
    gene_go_df$go
  )


  # For each GO term, find ancestors that belong to the GO-slim set

  go2slim_list <- lapply(
    go_unique,
    function(g) {

      anc <- anc_map[[g]]

      all_terms <- unique(
        c(g, anc)
      )

      intersect(
        all_terms,
        slim_ids
      )
    }
  )

  names(go2slim_list) <- go_unique


  go2slim_df <- tibble(
    go = names(go2slim_list),
    slim = unname(go2slim_list)
  ) %>%
    unnest(slim)


  gene_go_df %>%
    inner_join(
      go2slim_df,
      by = "go"
    ) %>%
    distinct(
      gene,
      slim
    )
}


gene_slim_bp <- build_gene_slim(
  gene_go2,
  goslim_ids,
  "BP"
)

gene_slim_mf <- build_gene_slim(
  gene_go2,
  goslim_ids,
  "MF"
)

gene_slim_cc <- build_gene_slim(
  gene_go2,
  goslim_ids,
  "CC"
)


# =============================================================================
# GO-SLIM CATEGORY ABUNDANCE
# =============================================================================

summarise_slim <- function(
  gene_slim,
  ont,
  denom
) {

  gene_slim %>%
    count(
      slim,
      name = "n_genes"
    ) %>%
    mutate(
      percent = 100 * n_genes / denom,
      ontology = ont,
      term = Term(GOTERM[slim])
    ) %>%
    filter(
      !is.na(term)
    ) %>%
    arrange(
      desc(percent)
    )
}


# Percentages are calculated relative to the number of transcripts
# represented in each ontology.

bp_sum <- summarise_slim(
  gene_slim_bp,
  "BP",
  length(unique(gene_slim_bp$gene))
)

mf_sum <- summarise_slim(
  gene_slim_mf,
  "MF",
  length(unique(gene_slim_mf$gene))
)

cc_sum <- summarise_slim(
  gene_slim_cc,
  "CC",
  length(unique(gene_slim_cc$gene))
)


go_slim_dist <- bind_rows(
  bp_sum,
  mf_sum,
  cc_sum
)


# =============================================================================
# PLOT TOP GO-SLIM CATEGORIES
# =============================================================================

topN <- 10

ont_order <- c(
  "BP",
  "MF",
  "CC"
)


plot_df <- go_slim_dist %>%
  mutate(
    ontology = factor(
      ontology,
      levels = ont_order
    )
  ) %>%
  group_by(
    ontology
  ) %>%
  slice_max(
    order_by = percent,
    n = topN
  ) %>%
  ungroup()


go_plot <- ggplot(
  plot_df,
  aes(
    x = reorder(term, percent),
    y = percent,
    fill = ontology
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  facet_grid(
    . ~ ontology,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      BP = "#8DA0CB",
      MF = "#FC8D62",
      CC = "#66C2A5"
    )
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = NULL,
    y = "Annotated transcripts (%)",
    title = NULL
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    strip.text.x = element_text(
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey92",
      color = "black"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1
    )
  )


go_plot

png(
  "./terms_transcriptome.png",
  width = 18 * 150,
  height = 13 * 150,
  res = 300
)

go_plot

dev.off()