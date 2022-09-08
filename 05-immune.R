# @DESCRIPTION: de

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(SummarizedExperiment)


# Load data ---------------------------------------------------------------

se <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "count_matrix.se.rds.gz"
  )
)

immune <- readr::read_tsv(
  file = "analysis/2022-07-20/result/integration/immune/all_expr_count_matrix_integration_ImmuCellAI_mouse_abundance_result.txt"
) %>% 
  dplyr::rename(barcode = `...1`) %>% 
  tidyr::gather(key = "cell_type", value = "score", - barcode) %>% 
  dplyr::inner_join(
    colData(se) %>% as.data.frame(),
    by = "barcode"
  )

readr::write_tsv(
  x = immune,
  file = "analysis/2022-07-20/result/integration/immune/se_immune.tsv"
)

immune %>% 
  dplyr::select(barcode, cell_type, score) %>% 
  tidyr::spread(key = cell_type, value = score) %>% 
  dplyr::arrange(barcode) %>% 
  writexl::write_xlsx(
    path = file.path(
      "analysis/2022-07-20/result/integration/immune",
      "tables5_immune_cell_abundance.xlsx"
    )
  )

readr::write_rds(
  x = immune,
  file = file.path(
    "analysis/2022-07-20/rda",
    "immune.rds.gz"
  )
)


# function ----------------------------------------------------------------

fn_plot_immune <- function(.immune, .levels = c("BC", "B2B", "B6B")) {
  .title <- paste0(.levels, collapse = "_")
  .immune %>% 
    dplyr::filter(cell_type != "Infiltration_score") %>% 
    dplyr::group_by(cell_type) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(f = purrr::map_lgl(.x = data, .f = function(.x){all(.x$score <0.01)})) %>% 
    # dplyr::filter(!f) %>% 
    dplyr::select(-f) %>% 
    dplyr::mutate(test = purrr::map(.x = data, .f = function(.x) {
      .xx <- kruskal.test(score ~ group, data = .x)
      .p <- broom::tidy(.xx)$p.value
      .m <- mean(.x$score)
      tibble::tibble(pval = .p, m = .m)
    })) %>% 
    dplyr::ungroup() %>% 
    tidyr::unnest(cols = test) %>% 
    dplyr::arrange(-m) %>% 
    # dplyr::mutate(cell_type = factor(x = cell_type, levels = cell_type)) %>% 
    dplyr::mutate(sig = ifelse(test = pval < 0.05, yes = "sig", no = "non-sig")) %>% 
    tidyr::unnest(cols = data) %>% 
    dplyr::mutate(group = factor(x = group, levels = .levels)) ->
    .immune_cell_for_plot
  
  .immune_cell_for_plot %>% 
    dplyr::select(cell_type, pval, m, sig) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(mark = dplyr::case_when(
      pval < 0.05 ~ "**",
      pval < 0.1 ~ "*",
      TRUE ~ ""
    )) %>% 
    dplyr::mutate(score = m + 0.2) ->
    .immune_cell_for_plot_mark
  
  
  .immune_cell_for_plot %>% 
    ggplot(aes(x = cell_type, y = score, color = group)) +
    geom_boxplot() +
    scale_color_manual(values = c("#00F5FF", "#CD0000",  "#191970"), name = "Group") +
    # scale_x_discrete(limit = .immune_cell_for_plot_mark$cell_type) +
    labs(x = "Immune cell types", y = "Immune infiltration score") +
    # annotate(
    #   geom = "text",
    #   x = .immune_cell_for_plot_mark$cell_type,
    #   y = .immune_cell_for_plot_mark$score,
    #   label = .immune_cell_for_plot_mark$mark
    # ) +
    theme(
      panel.background = element_rect(fill = NA),
      
      panel.grid = element_blank(),
      
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 16),
      axis.line.y.left = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      
      legend.position = c(0.9, 0.7),
      legend.background = element_blank(),
      legend.key = element_blank(),
    ) +
    labs(
      title = .title
    ) ->
    .immune_cell_plot
  list(
    plot = .immune_cell_plot,
    cell_type = .immune_cell_for_plot_mark
  )
  
}

# Plots -------------------------------------------------------------------
#BB
immune_BB_plot <- fn_plot_immune(
  .immune = immune %>% dplyr::filter(group %in% c("BC", "B2B", "B6B")), 
  .levels = c("BC", "B2B", "B6B")
  )

ggsave(
  filename = "immune_BB_plot.pdf",
  plot = immune_BB_plot$plot +
    scale_x_discrete(
      limit = immune_BB_plot$cell_type$cell_type
    ) +
    annotate(
      geom = "text",
      x = immune_BB_plot$cell_type$cell_type,
      y = immune_BB_plot$cell_type$score,
      label = immune_BB_plot$cell_type$mark
    ),
  path = "analysis/2022-07-20/result/integration/immune",
  width = 10, height = 4
)

# ss


immune_SS_plot <- fn_plot_immune(
  .immune = immune %>% dplyr::filter(group %in% c("SC", "S2S", "S4S")), 
  .levels = c("SC", "S2S", "S4S")
)


immune_SS_plot$cell_type %>% 
  dplyr::slice(match(immune_BB_plot$cell_type$cell_type, cell_type)) ->
  immune_SS_plot$cell_type


ggsave(
  filename = "immune_SS_plot.pdf",
  plot = immune_SS_plot$plot + 
    scale_x_discrete(
      limit = immune_SS_plot$cell_type$cell_type
    ) +
    annotate(
      geom = "text",
      x = immune_SS_plot$cell_type$cell_type,
      y = immune_SS_plot$cell_type$score,
      label = immune_SS_plot$cell_type$mark
    ),
  path = "analysis/2022-07-20/result/integration/immune",
  width = 10, height = 4
)

# BS


# immune %>% 
#   dplyr::filter(group %in% c("SC", "B2S", "B6S")) ->
#   immune_BS
# 
# immune_BS_plot <- fn_plot_immune(.immune = immune_BS, .levels = c("SC", "B2S", "B6S"))
# 


immune_BS_plot <- fn_plot_immune(
  .immune = immune %>% dplyr::filter(group %in% c("SC", "B2S", "B6S")), 
  .levels = c("SC", "B2S", "B6S")
)


immune_BS_plot$cell_type %>% 
  dplyr::slice(match(immune_BB_plot$cell_type$cell_type, cell_type)) ->
  immune_BS_plot$cell_type


ggsave(
  filename = "immune_BS_plot.pdf",
  plot = immune_BS_plot$plot + 
    scale_x_discrete(
      limit = immune_BS_plot$cell_type$cell_type
    ) +
    annotate(
      geom = "text",
      x = immune_BS_plot$cell_type$cell_type,
      y = immune_BS_plot$cell_type$score,
      label = immune_BS_plot$cell_type$mark
    ),
  path = "analysis/2022-07-20/result/integration/immune",
  width = 10, height = 4
)



# immune %>% 
#   dplyr::filter(group %in% c("BC", "S2B", "S4B")) ->
#   immune_SB
# 
# immune_SB_plot <- fn_plot_immune(.immune = immune_SB, .levels = c("BC", "S2B", "S4B"))

immune_SB_plot <- fn_plot_immune(
  .immune = immune %>% dplyr::filter(group %in% c("BC", "S2B", "S4B")), 
  .levels = c("BC", "S2B", "S4B")
)

immune_SB_plot$cell_type %>% 
  dplyr::slice(match(immune_BB_plot$cell_type$cell_type, cell_type)) ->
  immune_SB_plot$cell_type

ggsave(
  filename = "immune_SB_plot.pdf",
  plot = immune_SB_plot$plot + 
    scale_x_discrete(
      limit = immune_SB_plot$cell_type$cell_type
    ) +
    annotate(
      geom = "text",
      x = immune_SB_plot$cell_type$cell_type,
      y = immune_SB_plot$cell_type$score,
      label = immune_SB_plot$cell_type$mark
    ),
  path = "analysis/2022-07-20/result/integration/immune",
  width = 10, height = 4
)


# ggsave(
#   filename = "immune_SB_plot.pdf",
#   plot = immune_SB_plot,
#   path = "analysis/2022-07-20/result/integration/immune",
#   width = 7, height = 4
# )


# Individual --------------------------------------------------------------
bb <- c("BC", "B2B", "B6B")
ss <- c("SC", "S2S", "S4S")
bs <- c("SC", "B2S", "B6S")
sb <- c("BC", "S2B", "S4B")


colData(se) %>% as.data.frame() -> se_coldata

list(bb = bb, ss = ss, bs = bs, sb = sb) %>% 
  purrr::map(
    .f = function(.x) {
      se_coldata %>% 
        dplyr::filter(group %in% .x) ->
        .xx
      
      assay(se)[, .xx$barcode] %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "GeneName")
    }
  ) %>% 
  tibble::enframe(name = "group", value = "count") ->
  count_individual


purrr::pmap(.l = count_individual, .f = function(group, count) {
  .filename <- glue::glue("{group}_count.tsv")
  readr::write_tsv(
    x = count,
    file = file.path(
      "analysis/2022-07-20/result/integration/immune/individual",
      .filename
    )
  )
})



## load individual immune abundance and plots.

list(bb = bb, ss = ss, bs = bs, sb = sb) %>% 
  tibble::enframe(name = "group", value = "subgroup") %>% 
  dplyr::mutate(
    immune = purrr::map(
      .x = group,
      .f = function(.x) {
        # .x <- "bb"
        
        .immune <- readr::read_tsv(
          file = glue::glue("analysis/2022-07-20/result/integration/immune/individual/{.x}_count_ImmuCellAI_mouse_abundance_result.txt")
        ) %>% 
          dplyr::rename(barcode = `...1`) %>% 
          tidyr::gather(key = "cell_type", value = "score", - barcode) %>% 
          dplyr::inner_join(
            colData(se) %>% as.data.frame(),
            by = "barcode"
          )
      }
    )
  ) ->
  individual_count_immune

# individual_count_immune %>% 
#   dplyr::select(immune) %>% 
#   tidyr::unnest(cols = immune) %>% 
#   dplyr::select(barcode, cell_type, score, group) %>% 
#   dplyr::filter(cell_type == "B_cell") %>% 
#   dplyr::filter(group == "BC")
# 
# immune %>% 
#   dplyr::select(barcode, cell_type, score, group) %>% 
#   dplyr::filter(cell_type == "B_cell") %>% 
#   dplyr::filter(group == "BC")
  

readr::write_rds(
  x = individual_count_immune,
  file = file.path(
    "analysis/2022-07-20/rda",
    "individual_count_immune.rds.gz"
  )
)


individual_count_immune %>% 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(group, subgroup, immune), 
      .f = function(.x, .y, .z) {
        
        .p <- fn_plot_immune(.immune = .z, .levels = .y)
        .filename <- glue::glue("immune_{.x}_plot.pdf")
        ggsave(
          filename = .filename,
          plot = .p,
          path = "analysis/2022-07-20/result/integration/immune/individual/",
          width = 7, 
          height = 4
        )
        .p
    })
  ) ->
  individual_count_immune_plot

# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/05-immune.rda")
load(file = "analysis/2022-07-20/rda/05-immune.rda")