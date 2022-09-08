# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)


# Load data ---------------------------------------------------------------


a <- list.files("analysis/2021-11-05/raw", pattern = "All.xls", full.names = TRUE, recursive = TRUE)

b <- list.files("analysis/2022-02-08-bulk-rna-seq/raw", pattern = "All.xls", full.names = TRUE, recursive = TRUE)


c(a, b) %>% 
  tibble::enframe() %>% 
  dplyr::mutate(filename = basename(value)) %>% 
  dplyr::distinct(filename, .keep_all = T) ->
  filenames
  

# function ----------------------------------------------------------------

fn_loaddata <- function(.x, .y) {
  .d <- readr::read_tsv(file = .x)
  .d
}


# Add ---------------------------------------------------------------------



filenames %>% 
  dplyr::mutate(
    data = purrr::map2(
      .x = value, 
      .y= filename,
      .f = fn_loaddata
    )
  ) ->
  all_expr

all_expr %>%
  dplyr::mutate(
    vs = gsub(
      pattern = ".All.xls",
      replacement = "",
      x = filename
    )
  ) %>% 
  dplyr::mutate(
    vs = gsub(
      pattern = "VS",
      replacement = "_",
      x = vs
    )
  ) %>% 
  dplyr::select(
    vs,
    data, 
    path = value
  ) ->
  all_expr_update


# save --------------------------------------------------------------------

readr::write_rds(
  x = all_expr_update,
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr.rds.gz"
  )
)


# save image --------------------------------------------------------------


save.image(
  file = "analysis/2022-07-20/rda/01-loaddata.rda"
)
load("analysis/2022-07-20/rda/01-loaddata.rda")
