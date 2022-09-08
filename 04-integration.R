# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


source(file = "analysis/2022-07-20/script/utils.R")

# Load data ---------------------------------------------------------------

all_expr_count_matrix_integration <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr_count_matrix_integration.rds.gz"
  )
)

all_expr_count_matrix_integration %>% 
  dplyr::select(
    -dplyr::contains(match = "SS"),
    -dplyr::contains(match = "SD")
  ) ->
  expr_data

matchnames <- c("BD2" = "B2B", "BD6" = "B6B", "S2" = "S2S", "S4" = "S4S", "BS2" = "B2S", "BS6" = "B6S", "SB2" = "S2B", "SB4" = "S4B")

tibble::tibble(
  barcode = colnames(expr_data)[-1],
) %>% 
  dplyr::mutate(
    new_barcode = gsub(pattern = "BS2", replacement = "B2S", barcode),
  ) %>% 
  dplyr::mutate(
    new_barcode = gsub(pattern = "BD2", replacement = "B2B", new_barcode),
    new_barcode = gsub(pattern = "BD6", replacement = "B6B", new_barcode),
    new_barcode = gsub(pattern = "S2", replacement = "S2S", new_barcode),
    new_barcode = gsub(pattern = "S4", replacement = "S4S", new_barcode),
    new_barcode = gsub(pattern = "BS6", replacement = "B6S", new_barcode),
    new_barcode = gsub(pattern = "SB2", replacement = "S2B", new_barcode),
    new_barcode = gsub(pattern = "SB4", replacement = "S4B", new_barcode),
  ) %>% 
  dplyr::mutate(
    group = gsub(pattern = "_\\d", replacement = "", x = new_barcode)
  ) %>% 
  dplyr::mutate(
    seq = ifelse(group == "BC", "Brain", group)
  ) %>% 
  dplyr::mutate(
    seq = ifelse(seq == "SC", "Skull", seq),
    seq = ifelse(grepl(pattern = "S$", x = group), "Skull", seq),
    seq = ifelse(grepl(pattern = "B$", x = group), "Brain", seq)
  ) %>% 
  dplyr::mutate(
    case = ifelse(group == "BC", "Brain Control", group)
  ) %>% 
  dplyr::mutate(
    case = ifelse(case == "SC", "Skull Control", case),
    case = ifelse(grepl(pattern = "B2", x = case), "Brain2", case),
    case = ifelse(grepl(pattern = "B6", x = case), "Brain6", case),
    case = ifelse(grepl(pattern = "S2", x = case), "Skull2", case),
    case = ifelse(grepl(pattern = "S4", x = case), "Skull4", case)
  ) %>% 
  dplyr::mutate(
    batch = ifelse(group %in% c("BC", "SC", "B2B", "B6B","S2S", "S4S"), 1, 2)
  ) ->
  barcode_group
  
colnames(expr_data) <- c("GeneName", barcode_group$new_barcode)

readr::write_tsv(
  x = expr_data,
  file = file.path(
    "analysis/2022-07-20/result/integration/immune",
    "all_expr_count_matrix_integration_for_immune_raw_count.tsv"
  )
)

# Meta --------------------------------------------------------------------


# PCA ---------------------------------------------------------------------
colors <- RColorBrewer::brewer.pal(8, "Set2")

expr_data %>% 
  tibble::column_to_rownames("GeneName") %>% 
  as.matrix() ->
  expr_data.matrix

barcode_group %>% 
  # dplyr::select(-barcode) %>% 
  dplyr::mutate(barcode = new_barcode) %>% 
  tibble::column_to_rownames("new_barcode") ->
  expr_data.matrix.meta

brain_skull_pca <- PCAtools::pca(
  mat = expr_data.matrix,
  metadata = expr_data.matrix.meta,
  removeVar = 0.2
)

PCAtools::biplot(
  pcaobj = brain_skull_pca,
  colby = "seq",
  colLegendTitle = "Seq",
  shape = "case",
  shapeLegendTitle = "Case",
  legendPosition = "right",
  encircle = TRUE, 
  encircleFill = TRUE,
  lab = NULL
) ->
  brain_skull_pca_plot

ggsave(
  filename = "pca-brain-skull.pdf",
  plot = brain_skull_pca_plot,
  path = "analysis/2022-07-20/result/integration",
  device = "pdf",
  width = 7, height = 5
)

pdf(
  file = "analysis/2022-07-20/result/integration/rle-brain-skull-before.pdf", 
  width = 6, 
  height = 5
)
EDASeq::newSeqExpressionSet(
  counts = expr_data.matrix,
  phenoData = Biobase::AnnotatedDataFrame(
    data.frame(
      condition = factor(barcode_group$seq),
      row.names = barcode_group$new_barcode
    )
  )
) %>% 
  EDASeq::plotRLE(
  outline = FALSE,
  k = 2, 
  labels = TRUE,
  col = colors[factor(barcode_group$seq)],
  xaxt = "t", 
  las=2,
  ylab = 'RLE',
  ylim = c(-4,4), 
  main = 'Before removing batch effect'
)
dev.off()

# Batch effect removal ----------------------------------------------------


expr_data.matrix.be <- sva::ComBat_seq(
  counts = expr_data.matrix, 
  batch = barcode_group$batch
)

brain_skull_pca.be <- PCAtools::pca(
  mat = expr_data.matrix.be,
  metadata = expr_data.matrix.meta,
  removeVar = 0.2
)

PCAtools::biplot(
  pcaobj = brain_skull_pca.be,
  colby = "seq",
  colLegendTitle = "Seq",
  shape = "case",
  shapeLegendTitle = "Case",
  legendPosition = "right",
  encircle = TRUE, 
  encircleFill = TRUE,
  lab = NULL
) ->
  brain_skull_pca.be_plot

ggsave(
  filename = "pca-brain-skull-be.pdf",
  plot = brain_skull_pca.be_plot,
  path = "analysis/2022-07-20/result/integration",
  device = "pdf",
  width = 7, height = 5
)


pdf(
  file = "analysis/2022-07-20/result/integration/rle-brain-skull-after.pdf", 
  width = 6, 
  height = 5
)
EDASeq::newSeqExpressionSet(
  counts = expr_data.matrix.be,
  phenoData = Biobase::AnnotatedDataFrame(
    data.frame(
      condition = factor(barcode_group$seq),
      row.names = barcode_group$new_barcode
    )
  )
) %>% 
  EDASeq::plotRLE(
    outline = FALSE,
    k = 2, 
    labels = TRUE,
    col = colors[factor(barcode_group$seq)],
    xaxt = "t", 
    las=2,
    ylab = 'RLE',
    ylim = c(-4,4), 
    main = 'After removing batch effect'
  )
dev.off()

# new expr ----------------------------------------------------------------


se <- SummarizedExperiment::SummarizedExperiment(
  assays = expr_data.matrix.be,
  colData = expr_data.matrix.meta
)
readr::write_rds(
  x = se,
  file = file.path(
    "analysis/2022-07-20/rda",
    "count_matrix.se.rds.gz"
  )
)


# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/04-integration.rda")
