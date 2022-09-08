# @DESCRIPTION: all data

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


source(file = "analysis/2022-07-20/script/utils.R")

# Load data ---------------------------------------------------------------

all_expr <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr.rds.gz"
  )
)
msig_df <- msigdbr::msigdbr(species = "Mus musculus")

# select data -------------------------------------------------------------

purrr::map(
  .x = all_expr$data,
  .f = function(.x) {
    .x$GeneName
  }
) %>% 
  purrr::reduce(
    .f = intersect
  ) ->
  intersect_genes

purrr::map(
  .x = all_expr$data,
  .f = function(.x) {
    .x %>% 
      dplyr::filter(GeneName %in% intersect_genes) %>% 
      dplyr::select(1, 5:10) %>% 
      dplyr::distinct()
  }
) ->
  selected_samples

selected_samples %>% 
  purrr::reduce(
    .f = function(.x, .y) {
      .diff_names <- setdiff(colnames(.y), colnames(.x))
      .yy <- .y %>% 
        dplyr::select(1, .diff_names)
      .x %>% 
        dplyr::inner_join(.yy, by = "GeneName")
    }
  ) ->
  all_expr_count_matrix

readr::write_rds(
  x = all_expr_count_matrix,
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr_count_matrix.rds.gz"
  )
)



all_expr_count_matrix %>% 
  dplyr::select(
    1,
    dplyr::contains(match = "UVB"), 
    dplyr::contains(match = "UVS"),
    dplyr::contains(match = "BC"),
    dplyr::contains(match = "SC")
  ) ->
  all_expr_count_matrix_uv

readr::write_rds(
  x = all_expr_count_matrix_uv,
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr_count_matrix_uv.rds.gz"
  )
)

all_expr_count_matrix %>% 
  dplyr::select(
    1,
    dplyr::contains(match = "Q"),
    dplyr::contains(match = "SC")
  ) ->
  all_expr_count_matrix_qsc

readr::write_rds(
  x = all_expr_count_matrix_qsc,
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr_count_matrix_qsc.rds.gz"
  )
)

colnames(all_expr_count_matrix)
all_expr %>% 
  dplyr::pull()

all_expr_count_matrix %>% 
  dplyr::select(
    1, 
    dplyr::contains(match = "BC"),
    dplyr::contains(match = "SC"),
    dplyr::contains(match = "BD"),
    dplyr::contains(match = "SD"),
    dplyr::contains(match = "BS"),
    dplyr::contains(match = "SB"),
    dplyr::contains(match = "SS"),
    dplyr::contains(match = "S2"),
    dplyr::contains(match = "S4")
  ) ->
  all_expr_count_matrix_integration


readr::write_rds(
  x = all_expr_count_matrix_integration,
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr_count_matrix_integration.rds.gz"
  )
)

# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/03-all-count-matrix.rda")
