# @DESCRIPTION: immune new

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(DESeq2)
library(ComplexHeatmap)
library(ggupset)
library(fgsea)


# Load data ---------------------------------------------------------------

immune <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "immune.rds.gz"
  )
)

se <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "count_matrix.se.rds.gz"
  )
)

se_group_de_enrichment <- readr::read_rds(
  file = "analysis/2022-07-20/rda/07-enrichment.rds.gz"
) %>%
  tidyr::unnest(cols = enrichment)

se_gsva <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "se_hallmark_GSVA.rds.gz"
  )
)

hallmark_deg <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "hallmark_deg.se.rds.gz"
  )
)

gene_set_celldeath <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "gene_set_celldeath.rds.gz"
  )
)

color_group <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "color_group.rds.gz"
  )
)

msig_df <- msigdbr::msigdbr(species = "Mus musculus")

# Function ----------------------------------------------------------------


# DEG immune correlation --------------------------------------------------

se_group_de_enrichment %>%
  dplyr::select(bs, seq, se, vs, des_color) %>%
  dplyr::mutate(
    immune = purrr::map2(
      .x = se,
      .y = des_color,
      .f = function(.se, .des) {
        .des %>%
          dplyr::filter(color != "grey") ->
          .des_diff

        .dds <- DESeqDataSet(.se, design = ~ group)
        .dds <- estimateSizeFactors(.dds)
        .dds <- estimateDispersions(.dds)
        .dds <- nbinomWaldTest(.dds)

        .vsd <-  vst(.dds, blind=FALSE)

        assay(.vsd[.des_diff$GeneName,]) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "GeneName") %>%
          tidyr::gather(key = barcode, value = expr, - GeneName) %>%
          dplyr::inner_join(immune, by = "barcode") %>%
          dplyr::select(GeneName, cell_type, expr, score, group) %>%
          dplyr::group_by(GeneName, cell_type) %>%
          tidyr::nest() %>%
          dplyr::ungroup() ->
          .d
        future::plan(future::multisession, workers = 10)
        .d %>%
          dplyr::mutate(
            corr = furrr::future_map(
              .x = data,
              .f = function(.x) {
                cor.test(.x$expr, .x$score) %>%
                  broom::tidy() %>%
                  dplyr::select(corr = estimate, pval = p.value)
              }
            )
          ) %>%
          dplyr::select(-data) %>%
          tidyr::unnest(cols = corr)
      }
    )
  ) ->
  se_group_de_enrichment_immune_corr



# two group ---------------------------------------------------------------


se_group_de_enrichment_immune_corr %>% 
  dplyr::rename(immune_corr = immune) %>% 
  # head(1) %>% 
  dplyr::mutate(
    immune_biplot = purrr::map2(
      .x = vs,
      .y = seq,
      .f = function(.vs, .seq, .immune_data) {
        # .vs <- "B2B_vs_BC"
        # .seq <- "Brain"
        print(.vs)
        
        .vss <- strsplit(.vs, split = "_vs_")[[1]]
        .vs1 <- .vss[[1]]
        .vs2 <- .vss[[2]]
        
        .immune_data %>%
          dplyr::filter(group %in% .vss) %>% 
          dplyr::filter(cell_type != "Infiltration_score") %>%
          dplyr::group_by(cell_type) %>% 
          tidyr::nest() %>% 
          dplyr::mutate(f = purrr::map_lgl(
            .x = data, .f = function(.x){all(.x$score <0.01)})
            ) %>%
          dplyr::filter(!f) %>%
          dplyr::select(-f) %>%
          dplyr::mutate(
            test = purrr::map(
              .x = data,
              .f = function(.x) {
                .xx <- t.test(score ~ group, data = .x)
                .p <- broom::tidy(.xx)$p.value
                .e <- broom::tidy(.xx)$estimate
                .m <- mean(.x$score)
                tibble::tibble(pval = .p, m = .m, e = .e)
          })) %>%
          dplyr::ungroup() %>%
          tidyr::unnest(cols = test) %>% 
          dplyr::arrange(-m) %>%
          dplyr::mutate(
            cell_type = factor(
              x = cell_type,
              levels = cell_type
              )
            ) %>%
          dplyr::mutate(
            sig = ifelse(test = pval < 0.05, yes = "sig", no = "non-sig")
            ) %>%
          tidyr::unnest(cols = data) %>%
          dplyr::mutate(group = factor(x = group, levels = .vss)) ->
          .immune_cell_for_plot

        .immune_cell_for_plot %>%
          dplyr::select(cell_type, pval, m, e, sig) %>%
          dplyr::distinct() %>%
          dplyr::mutate(mark = dplyr::case_when(
            pval < 0.05 ~ "**",
            pval < 0.1 ~ "*",
            TRUE ~ ""
          )) %>%
          dplyr::mutate(score = m + 0.2) ->
          .immune_cell_for_plot_mark

        .title <- glue::glue("{.seq} - {.vs}")

        .immune_cell_for_plot %>%
          ggplot(aes(x = cell_type, y = score, color = group)) +
          geom_boxplot() +
          scale_color_manual(values = c("#00F5FF", "#CD0000") %>% rev(), name = "Group") +
          labs(x = "Immune cell types", y = "Immune infiltration score") +
          annotate(
            geom = "text",
            x = .immune_cell_for_plot_mark$cell_type,
            y = .immune_cell_for_plot_mark$score,
            label = .immune_cell_for_plot_mark$mark
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
          ) +
          labs(
            title = .title
          ) ->
          .immune_cell_plot
        ggsave(
          filename = glue::glue("{.seq}_{.vs}_immune_plot.pdf"),
          plot = .immune_cell_plot,
          path = "analysis/2022-07-20/result/integration/immune/twogroup",
          width = 7, height = 4
        )

        .immune_cell_for_plot_mark %>%
          # dplyr::filter(sig == "sig") %>%
          dplyr::select(cell_type, pval, m, e) %>%
          dplyr::arrange(pval)
      },
      .immune_data = immune
    )
  ) ->
  se_group_de_enrichment_immune_corr_plot

readr::write_rds(
  x = se_group_de_enrichment_immune_corr_plot,
  file = file.path(
    "analysis/2022-07-20/rda",
    "se_group_de_enrichment_immune_corr_plot.rds.gz"
  )
)

# se_group_de_enrichment_immune_corr_plot %>% 
#   dplyr::select(seq, vs, immune_corr, immune_biplot) %>% 
#   dplyr::mutate(
#     a = purrr::map2(
#       .x = immune_corr,
#       .y = immune_biplot,
#       .f = function(.corr, .bi) {
#         # .corr <- .d$immune_corr[[1]]
#         # .bi <- .d$immune_biplot[[1]]
#         
#         .corr %>% 
#           dplyr::filter(abs(corr) > 0.7, pval < 0.05) %>% 
#           dplyr::filter(cell_type %in% .bi$cell_type)
#         
#         
#         
#       }
#     )
#   )

# Immune bubble plot ------------------------------------------------------


se_group_de_enrichment_immune_corr_plot %>% 
  dplyr::select(seq, vs, immune_biplot) %>% 
  tidyr::unnest(cols = immune_biplot) %>% 
  dplyr::filter(pval < 0.05) %>% 
  dplyr::mutate(e = dplyr::case_when(
    e > 0.15 ~ 0.15,
    e < -0.15 ~ -0.15,
    TRUE ~ e
  )) %>% 
  dplyr::mutate(pval = pval) %>% 
  dplyr::mutate(pval_log = -log10(pval)) %>% 
  dplyr::mutate(pval_log = ifelse(pval_log > 4, 4, pval_log)) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", replacement = " vs. ", x = vs)) -> 
  se_group_de_enrichment_immune_corr_plot_for_bubble

se_group_de_enrichment_immune_corr_plot_for_bubble %>% 
  dplyr::mutate(
    pattern = purrr::map2_dbl(
      .x = e,
      .y = pval, 
      .f = function(.x, .y) {
        if ((.x > 0.1) && (.y < 0.05)) {
          return(1)
        } else if ((.x < -0.1) && (.y < 0.05)) {
          return(-1)
        } else {
          return(0)
        }
      }
    )
  ) %>% 
  dplyr::select(vs, cell_type, pattern) %>% 
  tidyr::spread(key = vs, value = pattern) %>% 
  dplyr::mutate_if(
    .predicate = is.numeric, 
    .fun = dplyr::funs(ifelse(is.na(.), 0, .))
  ) ->
  se_group_de_enrichment_immune_corr_plot_for_bubble_pattern

se_group_de_enrichment_immune_corr_plot_for_bubble %>% 
  dplyr::select(seq, vs) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(color = ifelse(seq == "Brain", "red", "black")) ->
  seq_color

se_group_de_enrichment_immune_corr_plot_for_bubble_pattern %>% 
  dplyr::rowwise() %>%
  dplyr::do(
    symbol = .$cell_type,
    rank =  unlist(.[-1], use.names = F) %>% abs %>% sum(),
    up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
    down = (unlist(.[-1], use.names = F) == -1) %>% sum()
  ) %>%
  dplyr::ungroup() %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(up_p = up / 15, down_p = down / 15, none = 1 - up_p - down_p) %>% 
  dplyr::arrange(-rank, -up) ->
  celltype_rank

se_group_de_enrichment_immune_corr_plot_for_bubble_pattern %>% 
  dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(.))) %>%
  tidyr::gather(key = vs, value = rank) %>%
  dplyr::left_join(seq_color, by = "vs") %>% 
  dplyr::arrange(color, rank) ->
  group_rank


se_group_de_enrichment_immune_corr_plot_for_bubble %>% 
  ggplot(aes(x = cell_type, y = vs, color = e, size = pval_log)) +
  geom_point() +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15),
    labels = c("-0.15", "-0.1", "-0.05", "0", "0.05",  "0.1", "0.15"),
    name = "Diff"
  ) +
  scale_size_continuous(
    name = "P value",
    limit = c(-log10(0.05), 4),
    range = c(1, 7),
    breaks = c(-log10(0.05), 2, 3, 4),
    labels = c("0.05", latex2exp::TeX("$10^{-2}$"), latex2exp::TeX("$10^{-3}$"), latex2exp::TeX("$< 10^{-4}$"))
  ) +
  scale_x_discrete(limit = celltype_rank$symbol) +
  scale_y_discrete(limit = group_rank$vs) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.text.y = element_text(color = group_rank$color, size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.margin = unit(c(10, 10, 10, 40), units = "points")
  ) +
  coord_cartesian(clip = 'off', xlim = c(1, nrow(se_group_de_enrichment_immune_corr_plot_for_bubble_pattern))) +
  annotate(geom = "segment",  x = -3.1, xend = -3.1, y = 9, yend = 16) +
  annotate(geom = "text", x = -3.5, y = 13, label = "Brain", size = 6, angle = 90) +
  annotate(geom = "segment",  x = -3.1, xend = -3.1, y = 1, yend = 8) +
  annotate(geom = "text", x = -3.5, y = 4, label = "Skull", size = 6, angle = 90)  ->
  immune_bubble_plot

ggsave(
  plot = immune_bubble_plot,
  filename = "immune_bubble_plot.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/immune",
  width = 12, 
  height = 6
)


# immune correlates with genes --------------------------------------------


se_group_de_enrichment_immune_corr_plot %>% 
  # dplyr::filter(grepl(pattern = "SC$|BC$", x = vs)) %>%
  # dplyr::filter(
  #   vs %in% c("B2B_vs_BC", "B6B_vs_B2B", "S2B_vs_BC", "S4B_vs_S2B", "S2S_vs_SC", "S4S_vs_S2S", "B2S_vs_SC", "B6S_vs_B2S")
  # ) %>%
  # dplyr::filter(
  #   vs %in% c("B2B_vs_BC", "S2S_vs_SC", "S2B_vs_BC",  "B6S_vs_SC")
  # ) %>%
  dplyr::filter(
    vs %in% c("B2B_vs_BC", "B6B_vs_B2B", "B6B_vs_BC", "S2S_vs_SC", "S4S_vs_S2S", "S4S_vs_SC")
  ) %>%
  dplyr::select(seq, vs, des_color, immune_corr, immune_biplot) %>% 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .ic = immune_corr,
        .ib = immune_biplot,
        .des = des_color
      ),
      .f = function(.ic, .ib, .des) {
        # .ic <- se_group_de_enrichment_immune_corr_plot$immune_corr[[15]]
        # .ib <- se_group_de_enrichment_immune_corr_plot$immune_biplot[[15]]
        # .des <- se_group_de_enrichment_immune_corr_plot$des_color[[15]]
        
        .ib %>%
          # dplyr::filter(cell_type %in% c("Macrophage", "NK", "B_cell")) %>% 
          dplyr::filter(m > 0.04, pval < 0.05) %>% 
          dplyr::mutate(imm_color = ifelse(e > 0, "red", "green")) %>% 
          dplyr::select(cell_type, imm_color) ->
          .ib_filter
        
        .des %>% 
          dplyr::filter(color != "grey") %>% 
          dplyr::select(GeneName, deg_color =  color) ->
          .des_filter
        
        .ic %>% 
          # dplyr::filter(cell_type %in% c("Macrophage", "NK", "B_cell")) %>% 
          dplyr::filter(cell_type %in% .ib_filter$cell_type) %>% 
          dplyr::filter(pval < 0.05, abs(corr) > 0.9) %>% 
          dplyr::mutate(corr_color = ifelse(corr > 0, "red", "green")) %>% 
          dplyr::select(GeneName, cell_type, corr_color) ->
          .ic_filter
        
        # .ic_filter %>% 
        #   dplyr::inner_join(.des_filter, by = "GeneName") %>% 
        #   dplyr::inner_join(.ib_filter, by = "cell_type")
        
        .ic_filter %>% 
          dplyr::filter(cell_type %in% c("Macrophage", "NK", "B_cell")) %>% 
          dplyr::group_by(cell_type) %>% 
          tidyr::nest() %>% 
          dplyr::ungroup() ->
          .b
        
        .ic_filter %>% 
          dplyr::group_by(cell_type, corr_color) %>% 
          dplyr::count() %>% 
          dplyr::ungroup() %>% 
          tidyr::spread(key = corr_color, value = n) -> 
          .a
        
        tibble::tibble(
          a = list(.a),
          b = list(.b)
        )
        
        
          
      }
    )
  ) %>% 
  tidyr::unnest(cols = a) ->
  se_group_de_enrichment_immune_corr_plot_filter_gene

{
  
  se_group_de_enrichment_immune_corr_plot_filter_gene %>% 
    dplyr::filter( vs %in% c("B2B_vs_BC", "B6B_vs_B2B", "B6B_vs_BC")) %>%
    # dplyr::filter( vs %in% c("S2S_vs_SC", "S4S_vs_S2S", "S4S_vs_SC")) %>%
    # dplyr::mutate(vs = factor(vs, levels = c("B2B_vs_BC", "B6B_vs_B2B", "B6B_vs_BC"))) %>% 
    dplyr::select(seq, vs, a) %>% 
    tidyr::unnest(cols = a) %>% 
    tidyr::gather(key = "type", value = "n", -c(seq, vs, cell_type)) %>% 
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>% 
    dplyr::mutate(np = ifelse(type == "red", "Positive", "Negative")) %>% 
    dplyr::mutate(np = factor(np, levels = c("Positive", "Negative"))) %>% 
    dplyr::mutate(vs = gsub(pattern = "_vs_", replacement = " vs. ", x = vs)) ->
    se_group_de_enrichment_immune_corr_plot_filter_gene_for_plot
  
  se_group_de_enrichment_immune_corr_plot_filter_gene_for_plot %>% 
    dplyr::group_by(cell_type) %>% 
    dplyr::summarise(m = sum(n)) %>% 
    dplyr::arrange(-m) ->
    cell_type_rank_pie
  
  se_group_de_enrichment_immune_corr_plot_filter_gene_for_plot %>% 
    dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_rank_pie$cell_type)) %>% 
    ggplot(aes(x = "", y = n, fill = np)) +
    geom_bar(stat = "identity") +
    coord_polar("y", start=0) +
    geom_text(
      aes(label = n),
      position = position_stack(vjust = 0.5)
    ) +
    scale_fill_brewer(
      name = "Correlation",
      palette = "Set1", 
      # direction = -1
    ) +
    theme(
      panel.background = element_rect(fill = "#F4F4F4", color = NA),
      panel.spacing = unit(0, units = "mm"),
      
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      
      # axis.line.x.bottom = element_line(color = "black"),
      # axis.line.y.left = element_line(color = "black"),
      
      strip.background = element_rect(fill = "white", color = NA),
      strip.text.x = element_text(angle = 90, size = 12, hjust = 1),
      strip.text.y = element_text(angle = 90, size = 10),
      legend.position = "top"
    ) +
    facet_grid(
      rows = dplyr::vars(vs), 
      cols = dplyr::vars(cell_type),
      switch = "both",
    ) ->
    brain_immune_grid_pie;brain_immune_grid_pie
  
  ggsave(
    filename = "Brain_Immune_gene_correlation.pdf",
    plot = brain_immune_grid_pie,
    device = "pdf",
    path = "analysis/2022-07-20/result/integration/immune",
    width = 10,
    height = 5
  )
}

{
  

  se_group_de_enrichment_immune_corr_plot_filter_gene %>% 
    # dplyr::filter( vs %in% c("B2B_vs_BC", "B6B_vs_B2B", "B6B_vs_BC")) %>%
    dplyr::filter( vs %in% c("S2S_vs_SC", "S4S_vs_S2S", "S4S_vs_SC")) %>%
    dplyr::select(seq, vs, a) %>% 
    tidyr::unnest(cols = a) %>% 
    tidyr::gather(key = "type", value = "n", -c(seq, vs, cell_type)) %>% 
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>% 
    dplyr::mutate(np = ifelse(type == "red", "Positive", "Negative")) %>% 
    dplyr::mutate(np = factor(np, levels = c("Positive", "Negative"))) %>% 
    dplyr::mutate(vs = gsub(pattern = "_vs_", replacement = " vs. ", x = vs)) ->
    se_group_de_enrichment_immune_corr_plot_filter_gene_for_plot_s
  
  se_group_de_enrichment_immune_corr_plot_filter_gene_for_plot_s %>% 
    dplyr::group_by(cell_type) %>% 
    dplyr::summarise(m = sum(n)) %>% 
    dplyr::arrange(-m) ->
    cell_type_rank_pie_s
  
  se_group_de_enrichment_immune_corr_plot_filter_gene_for_plot_s %>% 
    dplyr::mutate(cell_type = factor(cell_type, levels = cell_type_rank_pie_s$cell_type)) %>% 
    ggplot(aes(x = "", y = n, fill = np)) +
    geom_bar(stat = "identity") +
    coord_polar("y", start=0) +
    geom_text(
      aes(label = n),
      position = position_stack(vjust = 0.5)
    ) +
    scale_fill_brewer(
      name = "Correlation",
      palette = "Set1", 
      # direction = -1
    ) +
    theme(
      panel.background = element_rect(fill = "#F4F4F4", color = NA),
      panel.spacing = unit(0, units = "mm"),
      
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      
      # axis.line.x.bottom = element_line(color = "black"),
      # axis.line.y.left = element_line(color = "black"),
      
      strip.background = element_rect(fill = "white", color = NA),
      strip.text.x = element_text(angle = 90, size = 12, hjust = 1),
      strip.text.y = element_text(angle = 90, size = 10),
      legend.position = "top"
    ) +
    facet_grid(
      rows = dplyr::vars(vs), 
      cols = dplyr::vars(cell_type),
      switch = "both",
    ) ->
    skull_immune_grid_pie;skull_immune_grid_pie
  
  ggsave(
    filename = "Skull_Immune_gene_correlation.pdf",
    plot = skull_immune_grid_pie,
    device = "pdf",
    path = "analysis/2022-07-20/result/integration/immune",
    width = 10,
    height = 5
  )
}



# Filter genes ------------------------------------------------------------

hallmark_deg %>% 
  dplyr::filter(vs %in% c("B2B_vs_BC", "B6B_vs_B2B", "B6B_vs_BC")) %>% 
  dplyr::select(vs, GeneName, log2FC, FDR, color) %>%
  dplyr::group_by(vs) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::inner_join(
    se_group_de_enrichment_immune_corr_plot_filter_gene %>% 
      dplyr::filter(vs %in% c("B2B_vs_BC", "B6B_vs_B2B", "B6B_vs_BC")) %>% 
      dplyr::select(vs, b) ,
    by = "vs"
  ) ->
  se_group_de_enrichment_immune_corr_plot_filter_gene_hallmark_deg

se_group_de_enrichment_immune_corr_plot_filter_gene_hallmark_deg %>% 
  dplyr::filter(vs %in% c("B2B_vs_BC", "B6B_vs_BC")) %>% 
  dplyr::mutate(
    a = purrr::map2(
      .x = data,
      .y = b,
      .f = function(.x, .y) {
        # .x <- se_group_de_enrichment_immune_corr_plot_filter_gene_hallmark_deg$data[[1]]
        # .y <- se_group_de_enrichment_immune_corr_plot_filter_gene_hallmark_deg$b[[1]]
        .x %>% 
          dplyr::distinct() %>% 
          dplyr::arrange(-FDR, -abs(log2FC)) ->
          .xx
        
        .y %>% 
          tidyr::unnest(cols = data) %>% 
          dplyr::rename(color = corr_color) ->
          .yy
        
        .yy %>% 
          dplyr::inner_join(
            .xx,
            by = c("GeneName", "color")
          ) %>% 
          dplyr::arrange(-FDR, -abs(log2FC)) 
        
      }
    )
  ) ->
  se_group_de_enrichment_immune_corr_plot_filter_gene_hallmark_deg_inter

se_group_de_enrichment_immune_corr_plot_filter_gene_hallmark_deg_inter %>% 
  dplyr::select(vs, a) %>% 
  tidyr::spread(key = vs, value = a) %>% 
  dplyr::mutate(
    a = purrr::map2(
      .x = B2B_vs_BC,
      .y = B6B_vs_BC,
      .f = function(.x, .y) {
        .x <- d$B2B_vs_BC[[1]]
        .y <- d$B6B_vs_BC[[1]]
        
        .x %>% 
          dplyr::inner_join(
            .y,
            by = c("cell_type", "GeneName")
          ) %>% 
          dplyr::arrange(-FDR.x, -FDR.y)
        
      }
    )
  ) %>% 
  dplyr::select(a) %>% 
  tidyr::unnest(cols = a) %>% 
  dplyr::rename(
    color.x.B2B_vs_BC = color.x,
    color.y.B6B_vs_BC = color.y
  ) ->
  hallmark_deg_immune

hallmark_deg_immune %>% 
  dplyr::slice(1:50) %>% 
  dplyr::pull(GeneName) %>% 
  readr::write_lines(file = "analysis/2022-07-20/result/integration/candidate_genes.txt", sep = " ")
  






# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/09-immune2.rda")
load("analysis/2022-07-20/rda/09-immune2.rda")
