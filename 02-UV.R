# @DESCRIPTION: UV

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)


# Load data ---------------------------------------------------------------

all_expr <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "all_expr.rds.gz"
  )
)
msig_df <- msigdbr::msigdbr(species = "Mus musculus")

# function ----------------------------------------------------------------
source(file = "analysis/2022-07-20/script/utils.R")

# UVB ---------------------------------------------------------------------


all_expr %>% 
  dplyr::filter(grepl(pattern = "UVB", x = vs)) ->
 UVB

UVB %>% 
  dplyr::mutate(
    dfilter = purrr::map2(
      .x = data,
      .y = vs,
      .f = fn_selectdata,
      outpath = "analysis/2022-07-20/result/UV"
    )
  ) ->
  UVB_filter



UVB_filter$dfilter %>% fn_union_genes -> UVB_union_genes
UVB_filter$data %>% 
  fn_combine_data(.name = "UVB0_BC") ->
  UVB_filter_raw_count

UVB_filter_raw_count %>% 
  readr::write_tsv(
    file = "analysis/2022-07-20/result/UV/immune/BC_UVB0_UVB1_for_immune_raw_count.tsv"
  )

UVB_filter_raw_count %>% 
  dplyr::filter(GeneName %in% UVB_union_genes) ->
  UVB_filter_deg_union
  
UVB_filter_deg_union_heatmap <- fn_heatmap(UVB_filter_deg_union, .order = c("BC", "UVB0", "UVB1"))


ComplexHeatmap::row_order(UVB_filter_deg_union_heatmap) %>% 
  purrr::map(
    .f = function(.x) {
      UVB_filter_deg_union[.x, ]
    }
  ) ->
  UVB_filter_deg_union_heatmap_upgenes

names(UVB_filter_deg_union_heatmap_upgenes) <- c("BC", "UVB0", "UVB1")



# Immune ------------------------------------------------------------------

UVB_immune <- readr::read_tsv(
  file = "analysis/2022-07-20/result/UV/immune/ImmuCellAI_mouse_abundance_result_UVB.txt"
) %>% 
  dplyr::rename(sample = `...1`) %>% 
  tidyr::gather(key = "cell_type", value = "score", - sample) %>% 
  dplyr::mutate(group = gsub(pattern = "_\\d", replacement = "", x = sample))


UVB_immune %>% 
  dplyr::filter(cell_type != "Infiltration_score") %>% 
  dplyr::group_by(cell_type) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(f = purrr::map_lgl(.x = data, .f = function(.x){all(.x$score <0.01)})) %>% 
  dplyr::filter(!f) %>% 
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
  dplyr::mutate(cell_type = factor(x = cell_type, levels = cell_type)) %>% 
  dplyr::mutate(sig = ifelse(test = pval < 0.05, yes = "sig", no = "non-sig")) %>% 
  tidyr::unnest(cols = data) %>% 
  dplyr::mutate(group = factor(x = group, levels =c("BC", "UVB0", "UVB1"))) ->
  UVB_immune_cell_for_plot

UVB_immune_cell_for_plot %>% 
  dplyr::select(cell_type, pval, m, sig) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(mark = dplyr::case_when(
    pval < 0.05 ~ "**",
    pval < 0.1 ~ "*",
    TRUE ~ ""
  )) %>% 
  dplyr::mutate(score = m + 0.2) ->
  UVB_immune_cell_for_plot_mark


UVB_immune_cell_for_plot %>% 
  ggplot(aes(x = cell_type, y = score, color = group)) +
  geom_boxplot() +
  scale_color_manual(values = c("#00F5FF", "#CD0000",  "#191970"), name = "Group") +
  labs(x = "Immune cell types", y = "Immune infiltration score") +
  annotate(
    geom = "text",
    x = UVB_immune_cell_for_plot_mark$cell_type,
    y = UVB_immune_cell_for_plot_mark$score,
    label = UVB_immune_cell_for_plot_mark$mark
  ) +
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
  ) ->
  UVB_immune_cell_plot;UVB_immune_cell_plot

ggsave(
  filename = "UVB-immune-cell-distribution.pdf",
  plot = UVB_immune_cell_plot,
  path = "analysis/2022-07-20/result/UV/immune/",
  width = 7, height = 4
)


# Save heatmap ------------------------------------------------------------

pdf(file = "analysis/2022-07-20/result/UV/UVB_BC_deg_heatmap.pdf", width = 8, height = 10)
ComplexHeatmap::draw(object = UVB_filter_deg_union_heatmap)
dev.off()


UVB_filter %>% 
  dplyr::mutate(up = c("BC", "UVB0", "UVB1")) %>% 
  dplyr::left_join(
    UVB_filter_deg_union_heatmap_upgenes %>% 
      tibble::enframe(name = "up", value = "upgenes"), 
    by = "up"
  ) ->
  UVB_filter_upgenes

UVB_filter_upgenes %>% 
  dplyr::mutate(
    upd = purrr::map2(
      .x = data,
      .y = upgenes,
      .f = fn_upgenes
    )
  ) ->
  UVB_filter_upgenes_upd

UVB_filter_upgenes_upd %>% 
  dplyr::select(up, upd) %>% 
  tibble::deframe() %>% 
  writexl::write_xlsx(path = "analysis/2022-07-20/result/UV/UV-up-genes.xlsx")

UVB_filter_upgenes_upd %>% 
  dplyr::mutate(
    gobp = purrr::pmap(
      .l = list(data, upgenes, up),
      .f = fn_gobp,
      outpath = "analysis/2022-07-20/result/UV"
    )
  ) ->
  UVB_filter_upgenes_upd_gobp

# UVB_filter_upgenes_upd_gobp$gobp %>% 
#   unlist()


# UVS ---------------------------------------------------------------------


all_expr %>% 
  dplyr::filter(grepl(pattern = "U_", x = vs)) ->
  UVS

UVS %>% 
  dplyr::mutate(
    dfilter = purrr::map2(
      .x = data,
      .y = vs,
      .f = fn_selectdata,
      outpath = "analysis/2022-07-20/result/UV"
    )
  ) ->
  UVS_filter


UVS_filter$dfilter %>% fn_union_genes -> UVS_union_genes
UVS_filter$data %>% 
  fn_combine_data(.name = "UVS0_SC") ->
  UVS_filter_raw_count 

UVS_filter_raw_count %>% 
  readr::write_tsv(
    file = "analysis/2022-07-20/result/UV/immune/SC_UVS0_UVS1_for_immune_raw_count.tsv"
  )

UVS_filter_raw_count %>% 
  dplyr::filter(GeneName %in% UVS_union_genes) ->
  UVS_filter_deg_union

UVS_filter_deg_union_heatmap <- fn_heatmap(UVS_filter_deg_union, .order = c("SC", "UVS0", "UVS1"))


ComplexHeatmap::row_order(UVS_filter_deg_union_heatmap) %>% 
  purrr::map(
    .f = function(.x) {
      UVS_filter_deg_union[.x, ]
    }
  ) ->
  UVS_filter_deg_union_heatmap_upgenes

names(UVS_filter_deg_union_heatmap_upgenes) <- c("SC", "UVS0", "UVS1")

# Save heatmap ------------------------------------------------------------
pdf(file = "analysis/2022-07-20/result/UV/UVS_BC_deg_heatmap.pdf", width = 8, height = 10)
ComplexHeatmap::draw(object = UVS_filter_deg_union_heatmap)
dev.off()


UVS_filter %>% 
  dplyr::mutate(up = c("SC", "UVS0", "UVS1")) %>% 
  dplyr::left_join(
    UVS_filter_deg_union_heatmap_upgenes %>% 
      tibble::enframe(name = "up", value = "upgenes"), 
    by = "up"
  ) ->
  UVS_filter_upgenes

UVS_filter_upgenes %>% 
  dplyr::mutate(
    upd = purrr::map2(
      .x = data,
      .y = upgenes,
      .f = fn_upgenes
    )
  ) ->
  UVS_filter_upgenes_upd

UVS_filter_upgenes_upd %>% 
  dplyr::select(up, upd) %>% 
  tibble::deframe() %>% 
  writexl::write_xlsx(path = "analysis/2022-07-20/result/UV/UV-up-genes.xlsx")

UVS_filter_upgenes_upd %>% 
  dplyr::mutate(
    gobp = purrr::pmap(
      .l = list(data, upgenes, up),
      .f = fn_gobp,
      outpath = "analysis/2022-07-20/result/UV"
    )
  ) ->
  UVS_filter_upgenes_upd_gobp




# Immune ------------------------------------------------------------------

UVS_immune <- readr::read_tsv(
  file = "analysis/2022-07-20/result/UV/immune/ImmuCellAI_mouse_abundance_result_UVS.txt"
) %>% 
  dplyr::rename(sample = `...1`) %>% 
  tidyr::gather(key = "cell_type", value = "score", - sample) %>% 
  dplyr::mutate(group = gsub(pattern = "_\\d", replacement = "", x = sample))


UVS_immune %>% 
  dplyr::filter(cell_type != "Infiltration_score") %>% 
  dplyr::group_by(cell_type) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(f = purrr::map_lgl(.x = data, .f = function(.x){all(.x$score <0.01)})) %>% 
  dplyr::filter(!f) %>% 
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
  dplyr::mutate(cell_type = factor(x = cell_type, levels = cell_type)) %>% 
  dplyr::mutate(sig = ifelse(test = pval < 0.05, yes = "sig", no = "non-sig")) %>% 
  tidyr::unnest(cols = data) %>% 
  dplyr::mutate(group = factor(x = group, levels =c("SC", "UVS0", "UVS1"))) ->
  UVS_immune_cell_for_plot

UVS_immune_cell_for_plot %>% 
  dplyr::select(cell_type, pval, m, sig) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(mark = dplyr::case_when(
    pval < 0.05 ~ "**",
    pval < 0.1 ~ "*",
    TRUE ~ ""
  )) %>% 
  dplyr::mutate(score = m + 0.2) ->
  UVS_immune_cell_for_plot_mark


UVS_immune_cell_for_plot %>% 
  ggplot(aes(x = cell_type, y = score, color = group)) +
  geom_boxplot() +
  scale_color_manual(values = c("#00F5FF", "#CD0000",  "#191970"), name = "Group") +
  labs(x = "Immune cell types", y = "Immune infiltration score") +
  annotate(
    geom = "text",
    x = UVS_immune_cell_for_plot_mark$cell_type,
    y = UVS_immune_cell_for_plot_mark$score,
    label = UVS_immune_cell_for_plot_mark$mark
  ) +
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
  ) ->
  UVS_immune_cell_plot;UVS_immune_cell_plot

ggsave(
  filename = "UVS-immune-cell-distribution.pdf",
  plot = UVS_immune_cell_plot,
  path = "analysis/2022-07-20/result/UV/immune/",
  width = 7, height = 4
)



# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/02-uv.rda")
load(file = "analysis/2022-07-20/rda/02-uv.rda")
