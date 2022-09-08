# @DESCRIPTION: heatmap

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

se_group_de_enrichment <- readr::read_rds(
  file = "analysis/2022-07-20/rda/07-enrichment.rds.gz"
) %>% 
  tidyr::unnest(cols = enrichment)

se <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "count_matrix.se.rds.gz"
  )
)

msig_df <- msigdbr::msigdbr(species = "Mus musculus") 

msig_df %>% 
  dplyr::filter(gs_cat == "H") %>% 
  dplyr::select(gs_name, gene_symbol) ->
  msig_df_h

msig_df %>% 
  dplyr::filter(grepl(
    pattern = "ferroptosis",
    x = gs_name,
    ignore.case = T
  )) %>% 
  dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol) ->
  # dplyr::group_by(gs_cat, gs_subcat) %>% 
  # dplyr::count()
  msig_df_h_ferroptosis

msig_df %>% 
  dplyr::filter(grepl(
    pattern = "apoptosis",
    x = gs_name,
    ignore.case = T
  )) %>% 
  dplyr::filter(gs_subcat == "GO:BP") %>% 
  dplyr::select(gs_cat, gs_subcat , gs_name, gene_symbol) ->
  # dplyr::group_by(gs_cat, gs_subcat) %>%
  # dplyr::count() 
  msig_df_h_apoptosis

msig_df %>% 
  dplyr::filter(grepl(
    pattern = "Pyroptosis",
    x = gs_name,
    ignore.case = T
  )) %>% 
  dplyr::filter(gs_subcat == "GO:BP") %>%
  dplyr::select(gs_cat, gs_subcat , gs_name, gene_symbol) ->
  # dplyr::group_by(gs_cat, gs_subcat) %>% 
  # dplyr::count()
  msig_df_h_pyroptosis
  
msig_df %>% 
  dplyr::filter(grepl(
    pattern = "NECROPTOTIC",
    x = gs_name,
    ignore.case = T
  )) %>% 
  dplyr::filter(gs_subcat == "GO:BP") %>%
  dplyr::select(gs_cat, gs_subcat , gs_name, gene_symbol) ->
  # dplyr::group_by(gs_cat, gs_subcat) %>%
  # dplyr::count()
  msig_df_h_necroptosis

colors <- tibble::tibble(
  group = unique(se$group),
  color = c("black", "grey", RColorBrewer::brewer.pal(8, name = "Paired"))
)
readr::write_rds(
  x = colors,
  file = file.path(
    "analysis/2022-07-20/rda",
    "color_group.rds.gz"
  )
)

# Function ----------------------------------------------------------------

fn_heatmap_se <- function(.se, .type = "Brain") {
  # .se <- se_brain_deg
  
  
  .matrix <- assay(.se)
  .meta <- colData(.se)

  .matrix_scale <- .matrix %>% 
    apply(1, scale) %>%
    t()
  
  colnames(.matrix_scale) <- .meta$barcode
  
  colors %>% 
    dplyr::filter(group %in% .meta$group) %>% 
    tibble::deframe() ->
    .cluster_col
  
  hma_top = ComplexHeatmap::HeatmapAnnotation(
    df = as.data.frame(.meta) %>% 
      dplyr::select(Group = group),
    gap = unit(c(2,2), "mm"),
    col = list(Group = .cluster_col),
    which = "column"
  )
  
  ComplexHeatmap::Heatmap(
    # data and color
    matrix = .matrix_scale,
    col =  circlize::colorRamp2(
      breaks = c(-1.1, 0, 1.1), 
      colors = c("blue", "white", "red"), 
      space = "RGB"
    ),
    name = "Normalized counts",
    na_col = 'grey', 
    color_space = 'LAB', 
    rect_gp = gpar(col = NA),
    border = NA, 
    cell_fun = NULL, 
    layer_fun = NULL,
    jitter = FALSE,
    
    # title
    # row_title = 'Selected genes', # OC44
    row_title = 'Gene name',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    column_title = 'Samples',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 10),
    column_title_rot = 0,
    
    # clustering of row
    cluster_rows = T,
    cluster_row_slices = T,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D",
    row_dend_side = 'left',
    row_dend_width = unit(10, 'mm'),
    show_row_dend = F,
    row_dend_reorder = T,
    row_dend_gp = gpar(),
    
    # clustering of column
    cluster_columns = F,
    # cluster_column_slices = T,
    # clustering_distance_columns = "pearson",
    # clustering_method_columns = "ward.D",
    # column_dend_side = 'top',
    # column_dend_height = unit(10, 'mm'),
    # show_column_dend = F, # OC521
    # column_dend_gp = gpar(),
    # column_dend_reorder = F,
    
    row_order = NULL,
    # column_order = c(7,8,9,4,5,6,1,2,3),
    
    # row labels
    row_labels = rownames(.matrix_scale),
    row_names_side = 'right',
    show_row_names = F,
    row_names_max_width = unit(6, 'cm'),
    row_names_gp = gpar(fontsize = 3),
    row_names_rot = 0,
    row_names_centered = FALSE,
    
    # column labels
    # column_labels = colnames(.matrix_scale),
    column_names_side = 'bottom',
    show_column_names = T,
    column_names_max_height = unit(6, 'cm'),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_centered = FALSE,
    
    # annotation
    top_annotation = hma_top,
    bottom_annotation = NULL,
    # left_annotation = hma_right,
    # right_annotation = hma_left,
    left_annotation = NULL,
    right_annotation = NULL,
    
    
    # kmeans cluster number
    # row cluster is 1
    # column cluster is 2 with 10 repeats
    km = 1,
    split = NULL,
    row_km = 5,
    row_km_repeats = 5,
    row_split = NULL,
    # column_km = 3,
    # column_km_repeats = 10,
    column_split = rep(c("A", "B", "C", "D", "E"), each = 3), ###AAAAA
    # gap = unit(1, 'mm'),
    row_gap = unit(1, 'mm'),
    column_gap = unit(1, 'mm'),
    
    show_heatmap_legend = T,
    heatmap_legend_param = list(title = 'Row Z-score'),
    
    # raster_device = 'tiff',
    raster_quality = 2,
    raster_device_param = list(),
    
    raster_resize = F,
    
    post_fun = NULL
  ) -> 
    .heatmap
  
  .row_order <- ComplexHeatmap::row_order(.heatmap)
  
  .row_order %>% 
    tibble::enframe() %>% 
    dplyr::mutate(
      p = purrr::map2(
        .x = name,
        .y = value,
        .f = function(.rn, .rr) {
          # .rn <- .u$name[[1]]
          # .rr <- .u$value[[1]]
          
          .rr_genename <- rownames(.matrix_scale)[.rr]
          
          .go_bp <- clusterProfiler::enrichGO(
            gene = .rr_genename,
            universe = rownames(se),
            keyType = "SYMBOL",
            OrgDb = org.Mm.eg.db::org.Mm.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
          )
          
          .go_bp %>%
            tibble::as_tibble()  %>% 
            dplyr::mutate(
              Description = stringr::str_wrap(
                stringr::str_to_sentence(string = Description), width = 60
              )
            ) %>%
            dplyr::mutate(adjp = -log10(p.adjust)) %>%
            dplyr::select(ID, Description, adjp, Count) %>%
            head(20) %>%
            dplyr::arrange(adjp, Count) %>%
            dplyr::mutate(Description = factor(Description, levels = Description)) ->
            .go_bp_for_plot
          
          .go_bp_for_plot %>%
            ggplot(aes(x = Description, y = adjp)) +
            geom_col(fill = "#AE1700", color = NA, width = 0.7) +
            geom_text(aes(label = Count), hjust = 4, color = "white", size = 5) +
            labs(
              title = glue::glue("{.type} {.rn}; n=({length(.rr_genename)})"),
              y = "-log10(Adj. P value)"
            ) +
            scale_y_continuous(expand = c(0, 0.02)) +
            coord_flip() +
            theme(
              panel.background = element_rect(fill = NA),
              panel.grid = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(color = "black", size = 13, hjust = 1),
              axis.ticks.length.y = unit(3, units = "mm"),
              axis.text.x = element_text(color = "black", size = 12)
            ) ->
            .go_bp_plot
          
          .filename_prefix <- glue::glue("{.type}_heatmap_{.rn}_n=({length(.rr_genename)})")
          
          writexl::write_xlsx(
            x = as.data.frame(.go_bp), 
            path = file.path(
              "analysis/2022-07-20/result/integration/heatmap",
              glue::glue("{.filename_prefix}.xlsx")
            )
          )
          ggsave(
            plot = .go_bp_plot,
            filename = glue::glue("{.filename_prefix}.pdf"),
            device = "pdf",
            path = "analysis/2022-07-20/result/integration/heatmap",
            width = 10, 
            height = 6.5
          )
          
        }
      )
    )
    
  
  .heatmap
}

fn_heatmap_se_g <- function(.se, .g) {
  # .se <- se_brain_deg
  
  .matrix <- assay(.se)
  .meta <- colData(.se)
  
  .matrix_scale <- .matrix %>% 
    apply(1, scale) %>%
    t()
  
  colnames(.matrix_scale) <- .meta$barcode
  
  colors %>% 
    dplyr::filter(group %in% .meta$group) %>% 
    tibble::deframe() ->
    .cluster_col
  
  hma_top = ComplexHeatmap::HeatmapAnnotation(
    df = as.data.frame(.meta) %>% 
      dplyr::select(Group = group),
    gap = unit(c(2,2), "mm"),
    col = list(Group = .cluster_col),
    which = "column"
  )
  
  hma_index_right <- match(.g$g, rownames(.matrix_scale))
  hma_right <- ComplexHeatmap::rowAnnotation(
    link = anno_mark(
      at = hma_index_right,
      labels = rownames(.matrix_scale)[hma_index_right],
      which = "row", 
      side = "left",
      lines_gp = gpar(
        lwd = 0.5,
        col = .g$color
        ),
      labels_gp = gpar(
        fontsize = 7,
        col = .g$color
        ),
      padding = unit(0.5, "mm"),
      link_width = unit(5, "mm"),
    )
  )
  
  # hma_index_left <- match(BD2_mark, rownames(BD6_BD2_BC_deg_mat))
  # hma_left <- ComplexHeatmap::rowAnnotation(
  #   link = anno_mark(
  #     at = hma_index_left,
  #     labels = rownames(BD6_BD2_BC_deg_mat)[hma_index_left],
  #     which = "row",
  #     side = "right",
  #     lines_gp = gpar(lwd = 0.5, col = B6_B2_BC_color["B2"]),
  #     labels_gp = gpar(fontsize = 10, col = B6_B2_BC_color["B2"]),
  #     padding = unit(0.5, "mm"),
  #     link_width = unit(5, "mm"),
  #   )
  # )
  
  
  
  ComplexHeatmap::Heatmap(
    # data and color
    matrix = .matrix_scale,
    col =  circlize::colorRamp2(
      breaks = c(-1.1, 0, 1.1), 
      colors = c("blue", "white", "red"), 
      space = "RGB"
    ),
    name = "Normalized counts",
    na_col = 'grey', 
    color_space = 'LAB', 
    rect_gp = gpar(col = NA),
    border = NA, 
    cell_fun = NULL, 
    layer_fun = NULL,
    jitter = FALSE,
    
    # title
    # row_title = 'Selected genes', # OC44
    row_title = 'Gene name',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    column_title = 'Samples',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 10),
    column_title_rot = 0,
    
    # clustering of row
    cluster_rows = T,
    cluster_row_slices = T,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D",
    row_dend_side = 'left',
    row_dend_width = unit(10, 'mm'),
    show_row_dend = F,
    row_dend_reorder = T,
    row_dend_gp = gpar(),
    
    # clustering of column
    cluster_columns = F,
    # cluster_column_slices = T,
    # clustering_distance_columns = "pearson",
    # clustering_method_columns = "ward.D",
    # column_dend_side = 'top',
    # column_dend_height = unit(10, 'mm'),
    # show_column_dend = F, # OC521
    # column_dend_gp = gpar(),
    # column_dend_reorder = F,
    
    row_order = NULL,
    # column_order = c(7,8,9,4,5,6,1,2,3),
    
    # row labels
    row_labels = rownames(.matrix_scale),
    row_names_side = 'right',
    show_row_names = F,
    row_names_max_width = unit(6, 'cm'),
    row_names_gp = gpar(fontsize = 3),
    row_names_rot = 0,
    row_names_centered = FALSE,
    
    # column labels
    # column_labels = colnames(.matrix_scale),
    column_names_side = 'bottom',
    show_column_names = T,
    column_names_max_height = unit(6, 'cm'),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_centered = FALSE,
    
    # annotation
    top_annotation = hma_top,
    bottom_annotation = NULL,
    left_annotation = hma_right,
    # right_annotation = hma_right,
    
    
    # kmeans cluster number
    # row cluster is 1
    # column cluster is 2 with 10 repeats
    km = 1,
    split = NULL,
    row_km = 5,
    row_km_repeats = 5,
    row_split = NULL,
    # column_km = 3,
    # column_km_repeats = 10,
    column_split = rep(c("A", "B", "C", "D", "E"), each = 3),
    # gap = unit(1, 'mm'),
    row_gap = unit(1, 'mm'),
    column_gap = unit(1, 'mm'),
    
    show_heatmap_legend = T,
    heatmap_legend_param = list(title = 'Row Z-score'),
    
    # raster_device = 'tiff',
    raster_quality = 2,
    raster_device_param = list(),
    
    raster_resize = F,
    
    post_fun = NULL
  ) -> 
    .heatmap
  .heatmap
}

fn_heatmap_se_celldeath <- function(.se, .g) {
  # .se <- se_brain_deg
  
  .g %>% 
    dplyr::mutate(celldeath = factor(celldeath)) ->
    .gg
  
  .gg_color <- .gg$color
  
  names(.gg_color) <- gsub(pattern = " Pathway", replacement = "", x = names(.gg$color))
  
  
  
  .matrix <- assay(.se)
  
  .meta <- colData(.se) %>% 
    as.data.frame() %>% 
    dplyr::inner_join(
      colors,
      by = "group"
    ) %>% 
    dplyr::mutate(group = factor(group)) 
  
  .color <- .meta$color
  names(.color) <- .meta$group
  .meta <- .meta %>% 
    dplyr::mutate(color = .color)
  
  .matrix_scale <- .matrix %>% 
    apply(1, scale) 
  
  rownames(.matrix_scale) <- .meta$barcode
  
  
  
  ComplexHeatmap::Heatmap(
    matrix = .matrix_scale,
    col = circlize::colorRamp2(
      breaks = c(-1.3, 0, 1.3), 
      colors = c("#00fefe", "#000000", "#fe0000")
    ),
    name = "Normalized counts",
    # rect_gp = gpar(col = "white", lwd = 0.05),
    
    # column
    column_title = NULL,
    column_split = length(unique(.gg$color)),
    clustering_distance_columns = "pearson",
    cluster_columns = cluster_within_group(
      mat = .matrix_scale,
      factor = .gg %>%
        dplyr::select(GeneName, celldeath) %>%
        tibble::deframe()
    ),
    column_names_gp = gpar(fontsize = 8),
    bottom_annotation = HeatmapAnnotation(
      Pathway = names(.gg_color),
      col = list(Pathway = .gg_color),
      height = unit(2, "mm")
    ),
    show_column_names = TRUE,
    
    
    
    # row
    row_title = "Samples",
    row_gap = unit(0.05, "mm"),
    # row_dend_reorder = TRUE,
    # row_split = 10,
    clustering_distance_rows = "pearson",
    cluster_rows = cluster_within_group(
      mat = t(.matrix_scale),
      factor = .meta %>%
        dplyr::select(barcode, seq)  %>%
        tibble::deframe()
    ),
    # row_names_gp = gpar(fill = .meta$color),
    right_annotation = rowAnnotation(
      Group = names(.meta$color),
      col = list(Group = .meta$color),
      annotation_legend_param = list(
        Group = list(
          nrow = 2,
          title_position = "leftcenter"
          )
      )
    ),
    show_row_names = FALSE,
    
    # legend
    heatmap_legend_param = list(
      legend_width = unit(3, "cm"),
      direction = "horizontal"
      # title_position = "leftcenter-rot"
    )
  ) 
}

fn_gsea_running_score <- function(pathway, stats, gseaParam = 1) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- tibble::tibble(x = c(0, xs, n + 1), y = c(0, ys, 0))
  
  list(es = toPlot, pos = tibble::tibble(x = pathway))
}

fn_gseaplot <- function(es, pos, rank_genes = rank_genes, pathway = "H") {
  CPCOLS <- c("#ffffff", RColorBrewer::brewer.pal(9, "Set1"))
  ggplot() +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic(11) +
    theme(
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      legend.position = c(0.8, 0.8),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), 
      # axis.text.y = element_text(size = 12,colour = "black"),
      # axis.title = element_text(size=16),
      plot.margin = margin(
        t = 0.2, r = 0.2, b = 0, l = 0.2,
        unit = "cm"
      )
    ) + 
    labs(x = NULL, y = NULL) ->
    pproto
  
  pproto + 
    geom_line(data = es, aes(x = x, y = y), color = CPCOLS[4]) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    labs(x = NULL, y = "Enrichment score", title = glue::glue("{pathway}")) ->
    pup
  
  pos %>% dplyr::mutate(ymin = 0, ymax = 1) -> pos
  pos %>% 
    ggplot(aes(x = x)) +
    geom_linerange(aes(ymin = ymin, ymax = ymax)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic(11) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank(),
    ) +
    labs(x = NULL, y = NULL) ->
    pmiddle
  
  v <- seq(1, nrow(pos), length.out = 9)
  vt <- ifelse(is.na(match(seq_along(rank_genes), pos$x) ), 0, 1)
  inv <- findInterval(rev(cumsum(vt)), v)
  if (min(inv) == 0) {
    inv <- inv + 1
  }
  col <- c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
  ymin <- min(pmiddle$data$ymin)
  yy <- max(pmiddle$data$ymax - pmiddle$data$ymin) * 0.5
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  
  d <- data.frame(
    ymin = ymin, ymax = yy, xmin = xmin,
    xmax = xmax, col = col[unique(inv)]
  )
  
  pmiddle +
    geom_rect(
      data = d,
      aes(xmin = xmin, xmax = xmax, ymin =ymin, ymax = ymax, fill = I(col)),
      alpha = 0.9, inherit.aes = FALSE
    ) ->
    pmid
  
  df <- rank_genes %>% 
    tibble::enframe(name = 'x', value = 'y') %>% 
    dplyr::arrange(-y) %>% 
    dplyr::mutate(x = seq_along(y))
  
  pproto +
    geom_segment(data = df, aes(x = x, xend = x, y = y, yend = 0), color = 'grey') +
    theme(
      plot.margin = margin(
        t = -0.1, r = 0.2, b = 0.2,
        l = 0.2, unit = "cm"
      ),
      axis.line.x = element_line(),
      # axis.text.x = element_text(size = 12,colour = "black"),
      axis.ticks.x = element_line(),
    ) +
    labs(
      x = "Rank",
      y = "Ranked List Metric"
    ) ->
    pdown
  
  plotlist <- list(pup, pmid, pdown)
  
  cowplot::plot_grid(
    plotlist = plotlist,
    ncol = 1,
    align = 'v',
    rel_heights = c(1.5, 0.2, 1)
  )
}

fn_plot_single_gsea <- function(.x, .y, .d, .h) {
  # .x <- d$seq[[1]]
  # .y <- d$vs[[1]]
  # .d <- d$des_color[[1]]
  # .h <- d$H[[1]]
  
  .yy <- gsub(pattern = "_vs_", replacement = " vs. ", x = .y)
  
  .d %>% 
    dplyr::arrange(-log2FC) %>%
    dplyr::select(GeneName, log2FC) %>%
    tidyr::drop_na() %>%
    tibble::deframe() ->
    geneList
  
  .h %>% 
    dplyr::mutate(
      p = purrr::pmap(
        .l = list(
          des = Description,
          nes = NES,
          p.adj = p.adjust
        ),
        .f = function(des, nes, p.adj) {
          msig_df_h %>% 
            dplyr::filter(
              gs_name == des
            ) %>% 
            dplyr::pull(gene_symbol) ->
            .gs_genes
          .title <- glue::glue("{.x} - {.yy}\n{des}\nNES={round(nes,3)}\np.adjust={signif(p.adj,4)}")
          .es <- fn_gsea_running_score(pathway = .gs_genes, stats = geneList)
          .p <- fn_gseaplot(es = .es$es, pos = .es$pos, rank_genes = geneList, pathway = .title)
          ggsave(
            plot = .p,
            filename = glue::glue("{.x}_{.y}_{des}.pdf"),
            path = "analysis/2022-07-20/result/integration/go/gsea_single_plot",
            width = 5, 
            height = 5
          )
          .p
        }
      )
    ) 
  
}

fn_gsva <- function(.x, .y, .d) {
  # .x <- d$gs_name[[1]]
  # .y <- d$data[[1]]
  # .d <- se_count_matrix
  print(.x)
  .es <- GSVA::gsva(
    expr = .d, 
    gset.idx.list = list("geneset" = .y$gene_symbol),
    method = "ssgsea", 
    # parallel.sz = 10, 
    verbose=FALSE
  )
  
  as.data.frame(.es)
}

fn_gsva_plot <- function(.x, .y, .d) {
  # .x <- se_gsva$gs_name[[23]]
  # .y <- se_gsva$gsva[[23]]
  # 
  print(.x)
  
  .y %>% 
    tidyr::gather(key = "barcode", value = "gsva") %>% 
    dplyr::inner_join(.d, by = "barcode") ->
    .yy
  
  bb <- c("BC", "B2B", "B6B")
  ss <- c("SC", "S2S", "S4S")
  bs <- c("SC", "B2S", "B6S")
  sb <- c("BC", "S2B", "S4B")
  
  list(bb = bb, ss = ss, bs = bs, sb = sb) %>% 
    tibble::enframe() %>% 
    dplyr::mutate(
      p = purrr::map(
        .x = value,
        .f = function(.m) {
          print(.m)
          
          my_comparisons <-  list( .m[c(1, 2)], .m[c(2, 3)], .m[c(1, 3)] )
          
          .t.test <- my_comparisons %>% 
            tibble::enframe() %>% 
            dplyr::mutate(name = purrr::map_chr(.x = value, .f = function(.x) {
              paste0(.x %>% rev(), collapse = "_vs_")
            })) %>% 
            dplyr::mutate(
              pval = purrr::map_dbl(
                .x = value,
                .f = function(.mm) {
                  .t <- t.test(gsva ~ group, data = .yy %>% dplyr::filter(group %in% .mm))
                  .t$p.value
                }
              )
            ) %>% 
            dplyr::select(vs = name, pval_t_test = pval) 
          
          .yy %>% 
            dplyr::filter(group %in% .m) %>% 
            ggpubr::ggviolin(
              x = "group", 
              y = "gsva",  
              fill = "group", 
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "boxplot", 
              add.params = list(fill = "white")
            ) +
            ggpubr::stat_compare_means(
              comparisons = my_comparisons,
              method = "t.test",
            ) +
            theme(legend.position = "none") +
            labs(
              title = .x,
              x = "Group",
              y = "GSVA Score"
            ) ->
            .p
          
          tibble::tibble(
            pval_t_test = list(.t.test),
            p = list(.p)
          )
        }
      )
    ) %>% 
    dplyr::mutate(vs = purrr::map_chr(.x = value, .f = function(.x) {paste0(.x, collapse = "_vs_")})) %>% 
    dplyr::select(-value) %>% 
    tidyr::unnest(cols = p) ->
    .aa
  
  .aa %>% 
    dplyr::mutate(
      a = purrr::map2(
        .x = vs,
        .y = p,
        .f = function(.m, .p) {
          # .m <- "BC_vs_B2B_vs_B6B"
          .filename <- glue::glue("{.x}_{.m}.pdf")
          
          ggsave(
            filename = .filename,
            plot = .p,
            device = "pdf",
            path = "analysis/2022-07-20/result/integration/go/gsva_single_plot",
            width = 5, 
            height = 5
          )
        }
      )
    )
  
  .aa
  
}


# brain heatmap -----------------------------------------------------------

se_group_de_enrichment %>% 
  dplyr::filter(seq == "Brain") %>% 
  dplyr::mutate(
    deg = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color != "grey") %>% 
          dplyr::pull(GeneName)
      }
    ),
    deg_up = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color == "red") %>% 
          dplyr::pull(GeneName)
      }
    ),
    deg_down = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color == "green") %>% 
          dplyr::pull(GeneName)
      }
    )
  ) ->
  se_group_de_enrichment_brain


{
  pdf(
    file = "analysis/2022-07-20/result/integration/heatmap/Brain_heatmap.pdf", 
    width = 8, 
    height = 10
  )
  ComplexHeatmap::draw(
    object = fn_heatmap_se(
      .se = se[se_group_de_enrichment_brain$deg %>% purrr::reduce(.f = union), se$seq=="Brain"],
      .type = "Brain"
    )
  )
  dev.off()
}

# skull heatmap -----------------------------------------------------------

se_group_de_enrichment %>% 
  dplyr::filter(seq == "Skull") %>% 
  dplyr::mutate(
    deg = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color != "grey") %>% 
          dplyr::pull(GeneName)
      }
    ),
    deg_up = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color == "red") %>% 
          dplyr::pull(GeneName)
      }
    ),
    deg_down = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color == "green") %>% 
          dplyr::pull(GeneName)
      }
    )
  ) ->
  se_group_de_enrichment_skull


{
  pdf(
    file = "analysis/2022-07-20/result/integration/heatmap/Skull_heatmap.pdf", 
    width = 8, 
    height = 10
  )
  ComplexHeatmap::draw(
    object = fn_heatmap_se(
      .se = se[se_group_de_enrichment_skull$deg %>% purrr::reduce(.f = union), se$seq=="Skull"],
      .type = "Skull"
    )
  )
  dev.off()
}

# brain deg up ------------------------------------------------------------
de_color_up <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(5, 6)]
de_color_down <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(1, 2)]

se_group_de_enrichment_brain %>% 
  dplyr::select(vs, deg_up) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", " vs. ", x = vs)) %>% 
  tidyr::unnest(cols = deg_up) %>% 
  dplyr::group_by(deg_up) %>% 
  dplyr::summarize(Group = list(vs)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = Group, .f = length)) %>% 
  dplyr::arrange(-n) %>% 
  ggplot(aes(x = Group)) +
  geom_bar(width = 0.6, fill= de_color_up[1]) +
  scale_x_upset() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = de_color_up[1],
    combmatrix.panel.line.size = 0,
    combmatrix.label.make_space = TRUE,
  ) +
  labs(
    y = "# of Genes",
    x = "",
    title = "Brain Upnregulated"
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.8)
  ) ->
  upset_brain_up;upset_brain_up

ggsave(
  plot = upset_brain_up,
  filename = "UPSET_brain_up.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/de",
  width = 10, 
  height = 4
)

# brain deg down ----------------------------------------------------------


se_group_de_enrichment_brain %>% 
  dplyr::select(vs, deg_down) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", " vs. ", x = vs)) %>% 
  tidyr::unnest(cols = deg_down) %>% 
  dplyr::group_by(deg_down) %>% 
  dplyr::summarize(Group = list(vs)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = Group, .f = length)) %>% 
  dplyr::arrange(-n) %>% 
  ggplot(aes(x = Group)) +
  geom_bar(width = 0.6, fill= de_color_down[1]) +
  scale_x_upset() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = de_color_down[1],
    combmatrix.panel.line.size = 0,
    combmatrix.label.make_space = TRUE,
  ) +
  labs(
    y = "# of Genes",
    x = "",
    title = "Brain Downregulated"
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.8)
  ) ->
  upset_brain_down;upset_brain_down

ggsave(
  plot = upset_brain_down,
  filename = "UPSET_brain_down.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/de",
  width = 10, 
  height = 4
)

# Skull deg up ------------------------------------------------------------

se_group_de_enrichment_skull %>% 
  dplyr::select(vs, deg_up) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", " vs. ", x = vs)) %>% 
  tidyr::unnest(cols = deg_up) %>% 
  dplyr::group_by(deg_up) %>% 
  dplyr::summarize(Group = list(vs)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = Group, .f = length)) %>% 
  dplyr::arrange(-n) %>% 
  ggplot(aes(x = Group)) +
  geom_bar(width = 0.6, fill= de_color_up[2]) +
  scale_x_upset() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = de_color_up[2],
    combmatrix.panel.line.size = 0,
    combmatrix.label.make_space = TRUE,
  ) +
  labs(
    y = "# of Genes",
    x = "",
    title = "Skull Upnregulated"
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.8)
  ) ->
  upset_skull_up;upset_skull_up

ggsave(
  plot = upset_skull_up,
  filename = "UPSET_skull_up.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/de",
  width = 10, 
  height = 4
)

# skull deg down ----------------------------------------------------------


se_group_de_enrichment_skull %>% 
  dplyr::select(vs, deg_down) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", " vs. ", x = vs)) %>% 
  tidyr::unnest(cols = deg_down) %>% 
  dplyr::group_by(deg_down) %>% 
  dplyr::summarize(Group = list(vs)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = Group, .f = length)) %>% 
  dplyr::arrange(-n) %>% 
  ggplot(aes(x = Group)) +
  geom_bar(width = 0.6, fill= de_color_down[2]) +
  scale_x_upset() +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = de_color_down[2],
    combmatrix.panel.line.size = 0,
    combmatrix.label.make_space = TRUE,
  ) +
  labs(
    y = "# of Genes",
    x = "",
    title = "Skull Downregulated"
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.8)
  ) ->
  upset_skull_down;upset_skull_down

ggsave(
  plot = upset_skull_down,
  filename = "UPSET_skull_down.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/de",
  width = 10, 
  height = 4
)


# Brain skull comparison and gobp -----------------------------------------

comp <- c("B2B_vs_BC#S2S_vs_SC", "B6B_vs_BC#S4S_vs_SC", "B6B_vs_B2B#S4S_vs_S2S", "S2B_vs_BC#B2S_vs_SC", "S4B_vs_BC#B6S_vs_SC", "S4B_vs_S2B#B6S_vs_B2S", "S2B_vs_B2B#B2S_vs_S2S", "S4B_vs_B6B#B6S_vs_S4S")


se_group_de_enrichment %>% 
  dplyr::select(vs, des_color) %>% 
  dplyr::mutate(comp = rep(comp, 2)) %>% 
  dplyr::group_by(comp) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() ->
  se_group_de_enrichment_des_nest

se_group_de_enrichment_des_nest %>% 
  dplyr::mutate(
    tileplot = purrr::map2(
      .x = comp,
      .y = data,
      .f = function(.x, .y) {
        # .x <- se_group_de_enrichment_des_nest$comp[[1]]
        # .y <- se_group_de_enrichment_des_nest$data[[1]]
        
        .vs <- strsplit(x = .x, split = "#")[[1]]
        .b <- .vs[[1]] %>% 
          gsub(
            pattern = "_vs_", 
            replacement = " vs. ", 
            x = .
          )
        .s <- .vs[[2]] %>% 
          gsub(
            pattern = "_vs_", 
            replacement = " vs. ", 
            x = .
          )
        
        .y %>% 
          dplyr::mutate(
            des_color = purrr::map(
              .x = des_color,
              .f = function(.x) {
                .x %>% 
                  dplyr::select(GeneName, color) %>% 
                  dplyr::mutate(color = plyr::revalue(
                    x = color,
                    replace = c("red" = "Up", "green" = "Down", "grey" = "No Sig.")
                  ))
              }
            )
          ) %>% 
          tidyr::spread(key = vs, value = des_color) %>% 
          dplyr::select(dplyr::all_of(.vs)) ->
          .yy
        
        .bd <- .yy[[1]][[1]]
        .sd <- .yy[[2]][[1]]
        
        .bd %>% 
          dplyr::inner_join(.sd, by = "GeneName") %>% 
          dplyr::select(GeneName, color.x, color.y) %>% 
          dplyr::filter(!(color.x == "No Sig." & color.y == "No Sig.")) ->
          .bsd
        
        CPCOLS <- c("#191970", "#F8F8FF", "#FF4040")
        
        .bsd %>% 
          dplyr::group_by(color.x, color.y) %>% 
          tidyr::nest() %>% 
          dplyr::ungroup() %>% 
          dplyr::mutate(n = purrr::map_int(.x = data, .f = nrow)) %>% 
          dplyr::ungroup() ->
          .bsdn
        
        .bsdn %>% 
          ggplot(aes(x = color.x, y = color.y)) +
          geom_tile(aes(fill = n, )) +
          geom_text(aes(label = n, ), color = "black") +
          scale_fill_gradient2(
            low = CPCOLS[1], 
            mid = CPCOLS[2], 
            high = CPCOLS[3],
            name = "Count"
          ) +
          theme(
            axis.ticks = element_blank(),
            axis.text = element_text(size = 12, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.5),
            axis.title = element_text(size = 14, color = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
          ) +
          coord_fixed() +
          labs(
            x = .b,
            y = .s
          ) ->
          .p
        
        ggsave(
          filename = glue::glue("{.x}.pdf"),
          plot = .p,
          path = "analysis/2022-07-20/result/integration/de/comparison",
          device = "pdf",
          width = 5, 
          height = 5
        )
        
        .bsdn %>% 
          dplyr::filter(n > 100) %>% 
          dplyr::mutate(
            gobp = purrr::pmap(
              .l = list(
                .bb = color.x,
                .ss = color.y,
                .dd = data,
                .nn = n
              ),
              .f = function(.bb, .ss, .dd, .nn) {
                # .bb <- .bsdn$color.x[[1]]
                # .ss <- .bsdn$color.y[[1]]
                # .dd <- .bsdn$data[[1]]
                # .nn <- .bsdn$n[[1]]
                .cc <- c("red" = "#AE1700", "green" = "#112a13")
                .color <- ifelse(any(c(.bb, .ss) == "Up"), "red", "green") 
                
                .go_bp <- clusterProfiler::enrichGO(
                  gene = .dd$GeneName,
                  universe = rownames(se),
                  keyType = "SYMBOL",
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                )
                .go_bp %>%
                  tibble::as_tibble()  %>% 
                  dplyr::mutate(
                    Description = stringr::str_wrap(
                      stringr::str_to_sentence(string = Description), width = 60
                    )
                  ) %>%
                  dplyr::mutate(adjp = -log10(p.adjust)) %>%
                  dplyr::select(ID, Description, adjp, Count) %>%
                  head(20) %>%
                  dplyr::arrange(adjp, Count) %>%
                  dplyr::mutate(Description = factor(Description, levels = Description)) ->
                  .go_bp_for_plot
                
                .go_bp_for_plot %>%
                  ggplot(aes(x = Description, y = adjp)) +
                  geom_col(fill =  .cc[.color], color = NA, width = 0.7) +
                  geom_text(aes(label = Count), hjust = 4, color = "white", size = 5) +
                  labs(
                    title = glue::glue("{.b} ({.bb}) & {.s} ({.ss}) n=({.nn})"),
                    y = "-log10(Adj. P value)"
                    ) +
                  scale_y_continuous(expand = c(0, 0.02)) +
                  coord_flip() +
                  theme(
                    panel.background = element_rect(fill = NA),
                    panel.grid = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"),
                    axis.title.y = element_blank(),
                    axis.text.y = element_text(color = "black", size = 13, hjust = 1),
                    axis.ticks.length.y = unit(3, units = "mm"),
                    axis.text.x = element_text(color = "black", size = 12)
                  ) ->
                  .go_bp_plot
                
                
                .filename_prefix <- glue::glue("{.x}___{.b} ({.bb}) & {.s} ({.ss}) n=({.nn})")
                
                writexl::write_xlsx(
                  x = as.data.frame(.go_bp), 
                  path = file.path(
                    "analysis/2022-07-20/result/integration/de/comparison",
                    glue::glue("{.filename_prefix}.xlsx")
                  )
                )
                ggsave(
                  plot = .go_bp_plot,
                  filename = glue::glue("{.filename_prefix}.pdf"),
                  device = "pdf",
                  path = "analysis/2022-07-20/result/integration/de/comparison",
                  width = 10, 
                  height = 6.5
                )
              }
            )
          )
        
        .p
      }
    )
  ) ->
  se_group_de_enrichment_des_nest_tile_plot


se_group_de_enrichment_des_nest_tile_plot$tileplot[4:8] %>% 
  purrr::reduce(`+`) +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A") ->
  se_group_de_enrichment_des_nest_tile_plot_patch


ggsave(
  filename = "Brain_Skull_comparison_n_gene.pdf",
  plot = se_group_de_enrichment_des_nest_tile_plot_patch,
  path = "analysis/2022-07-20/result/integration/de/comparison",
  device = "pdf",
  width = 20, 
  height = 10
)

# gsea --------------------------------------------------------------------



se_group_de_enrichment %>% 
  dplyr::mutate(
    H = purrr::map(
      .x = gsea, 
      .f = function(.x) {
        .x$Cat %>% 
          table()
        .x %>% 
          dplyr::filter(Cat == "H") %>% 
          dplyr::filter(p.adjust < 0.05, abs(NES) > 1.5) %>% 
          dplyr::select(Description, NES, p.adjust, core_enrichment)
      }
    )
  ) ->
  se_group_de_enrichment_gsea

se_group_de_enrichment_gsea %>% 
  dplyr::mutate(
    Hn = purrr::map(
      .x = H, 
      .f = function(.x) {
        .n = nrow(.x)
        tibble::tibble(
          des = list(.x$Description),
          n = .n
        )
      }
    )
  ) %>% 
  tidyr::unnest(Hn) ->
  se_group_de_enrichment_gsea_des


se_group_de_enrichment_gsea_des %>% 
  dplyr::select(vs, des) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", " vs. ", x = vs)) %>% 
  tidyr::unnest(cols = des) ->
  se_group_de_enrichment_gsea_des_hallmark_group

se_group_de_enrichment_gsea_des_hallmark_group %>% 
  dplyr::group_by(des) %>% 
  dplyr::summarize(Group = list(vs)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = Group, .f = length)) %>% 
  dplyr::arrange(-n) %>% 
  ggplot(aes(x = Group, fill = n)) +
  geom_bar(width = 0.6) +
  scale_x_upset() +
  scale_fill_continuous(
    limits = c(0, 15),
    breaks = c(0, 3, 6, 9, 12, 15),
    guide = guide_colorbar(
      title = "# of Group",
      title.position = "top",
      title.theme = element_text(hjust = 0.5),
      barheight = 0.7,
      direction = "horizontal"
    )
  ) +
  labs(
    y = "# of Hallmarks",
    x = ""
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.8)
  ) ->
  upset_group;upset_group

ggsave(
  plot = upset_group,
  filename = "UPSET_group.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/go/gsea_upset",
  width = 10, 
  height = 6
)

se_group_de_enrichment_gsea_des_hallmark_group %>% 
  dplyr::group_by(vs) %>% 
  dplyr::summarise(H = list(des)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = H, .f = length)) %>%
  dplyr::arrange(-n) %>% 
  ggplot(aes(x = H, fill = n)) +
  geom_bar(width = 0.6) +
  scale_x_upset() +
  scale_fill_continuous(
    limits = c(0, 31),
    breaks = c(0, 10, 20, 30),
    guide = guide_colorbar(
      title = "# of Hallmarks",
      title.position = "top",
      title.theme = element_text(hjust = 0.5),
      barheight = 0.7,
      direction = "horizontal"
    )
  ) +
  labs(
    y = "# of Groups",
    x = ""
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.75)
  ) ->
  upset_h;upset_h

ggsave(
  plot = upset_h,
  filename = "UPSET_hallmark.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/go/gsea_upset",
  width = 13, 
  height = 10
)


# GSEA plot ---------------------------------------------------------------


se_group_de_enrichment_gsea_des_hallmark_group %>% 
  dplyr::group_by(des) %>% 
  dplyr::summarize(Group = list(vs)) %>% 
  dplyr::mutate(n = purrr::map_int(.x = Group, .f = length)) %>% 
  dplyr::arrange(-n)

se_group_de_enrichment_gsea %>% 
  dplyr::filter(purrr::map_lgl(.x = H, .f = function(.x){nrow(.x) != 0})) %>% 
  dplyr::select(seq, vs, des_color, H) %>% 
  dplyr::mutate(
    gsea_plot = purrr::pmap(
      .l = list(
        .x = seq,
        .y = vs,
        .d = des_color, 
        .h = H
      ),
      .f = fn_plot_single_gsea
    )
  ) ->
  se_group_de_enrichment_gsea_plot




# GSEA bubble plot --------------------------------------------------------

se_group_de_enrichment %>% 
  dplyr::mutate(
    H = purrr::map(
      .x = gsea, 
      .f = function(.x) {
        .x %>% 
          dplyr::filter(Cat == "H") 
      }
    )
  ) %>% 
  # dplyr::filter(purrr::map_lgl(.x = H, .f = function(.x){nrow(.x) != 0})) %>% 
  dplyr::select(seq, vs, H) %>% 
  tidyr::unnest(H) %>% 
  dplyr::select(seq, vs, Description, NES, pvalue, p.adjust) %>% 
  dplyr::mutate(des = gsub(pattern = "HALLMARK_", replacement = "", x = Description)) %>% 
  dplyr::mutate(des = gsub(pattern = "_", replacement = " ", x = des)) %>% 
  dplyr::mutate(p.adj  = -log10(p.adjust)) %>% 
  dplyr::mutate(nes_up = ifelse(NES > 0, "up", "down")) %>% 
  dplyr::mutate(vs = gsub(pattern = "_vs_", replacement = " vs. ", x = vs)) ->
  se_group_de_enrichment_gsea_for_plot

se_group_de_enrichment_gsea_for_plot %>% 
  dplyr::mutate(
    pattern = purrr::map2_dbl(
      .x = NES,
      .y = p.adjust, 
      .f = function(.x, .y) {
        if ((.x > 1.3) && (.y < 0.05)) {
          return(1)
        } else if ((.x < -1.3) && (.y < 0.05)) {
          return(-1)
        } else {
          return(0)
        }
      }
    )
  ) %>% 
  dplyr::select(vs, des, pattern) %>% 
  tidyr::spread(key = vs, value = pattern) %>% 
  dplyr::mutate_if(
    .predicate = is.numeric, 
    .fun = dplyr::funs(ifelse(is.na(.), 0, .))
  ) ->
  se_group_de_enrichment_gsea_for_plot_pattern


se_group_de_enrichment_gsea_for_plot %>% 
  dplyr::select(seq, vs) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(color = ifelse(seq == "Brain", "red", "black")) ->
  seq_color

se_group_de_enrichment_gsea_for_plot_pattern %>% 
  dplyr::rowwise() %>%
  dplyr::do(
    symbol = .$des,
    rank =  unlist(.[-1], use.names = F) %>% abs %>% sum(),
    up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
    down = (unlist(.[-1], use.names = F) == -1) %>% sum()
  ) %>%
  dplyr::ungroup() %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(up_p = up / 15, down_p = down / 15, none = 1 - up_p - down_p) %>% 
  dplyr::arrange(-rank, -up) ->
  des_rank

se_group_de_enrichment_gsea_for_plot_pattern %>% 
  dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(.))) %>%
  tidyr::gather(key = vs, value = rank) %>%
  dplyr::left_join(seq_color, by = "vs") %>% 
  dplyr::arrange(color, rank) ->
  group_rank


se_group_de_enrichment_gsea_for_plot %>% 
  ggplot(aes(x = des, y = vs, color = NES, size = p.adj)) +
  geom_point() +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-3, 3, length.out = 7),
    labels = c("-3", "-2", "-1", "0", "1",  "2", "3"),
    name = "NES"
  ) +
  scale_size_continuous(
    name = "p.adjust",
    limit = c(-log10(0.01), 9),
    range = c(1, 6),
    breaks = c(-log10(0.01), 3, 6, 9),
    labels = c("0.01", latex2exp::TeX("$10^{-3}$"), latex2exp::TeX("$10^{-6}$"), latex2exp::TeX("$< 10^{-9}$"))
  ) +
  scale_x_discrete(limit = des_rank$symbol) +
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
  coord_cartesian(clip = 'off', xlim = c(1, 39)) +
  annotate(geom = "segment",  x = -2.7, xend = -2.7, y = 8, yend = 15) +
  annotate(geom = "text", x = -3.2, y = 11, label = "Brain", size = 6, angle = 90) +
  annotate(geom = "segment",  x = -2.7, xend = -2.7, y = 1, yend = 7) +
  annotate(geom = "text", x = -3.2, y = 4, label = "Skull", size = 6, angle = 90) ->
  gsea_bubble_plot

ggsave(
  plot = gsea_bubble_plot,
  filename = "HALLMARK_gsea_bubble_plot.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/go/gsea_upset",
  width = 15, 
  height = 8
)



# GSVA --------------------------------------------------------------------

colData(se) %>% as.data.frame() -> se_coldata
assay(se) -> se_count_matrix

msig_df_h %>% 
  dplyr::group_by(gs_name) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    gsva = purrr::map2(
      .x = gs_name,
      .y = data,
      .f = fn_gsva,
      .d = se_count_matrix
    )
  ) %>% 
  dplyr::select(-data) ->
  se_gsva


se_gsva %>% 
  dplyr::mutate(
    p = purrr::map2(
      .x = gs_name,
      .y = gsva,
      .f = fn_gsva_plot,
      .d = se_coldata
    )
  ) %>% 
  dplyr::select(-gsva) %>% 
  tidyr::unnest(cols = p) ->
  se_gsva_plot


# save gsva score ---------------------------------------------------------

se_gsva %>% 
  tidyr::unnest(cols = gsva) %>% 
  readr::write_rds(
    file.path(
      "analysis/2022-07-20/rda",
      "se_hallmark_GSVA.rds.gz"
    )
  )

list(
  `HALLMARK GSVA Score` = se_gsva %>% 
    tidyr::unnest(cols = gsva) ,
  `HALLMARK GSVA Score t test P vlaue` = se_gsva_plot %>% 
    dplyr::select(gs_name, pval_t_test) %>% 
    tidyr::unnest(cols = pval_t_test) %>% 
    tidyr::spread(key = gs_name, value = pval_t_test) 
) %>% 
  writexl::write_xlsx(
    path = "analysis/2022-07-20/result/integration/go/GSVA_hallmark.xlsx"
  )






# GSEA gene select --------------------------------------------------------

se_group_de_enrichment %>% 
  dplyr::mutate(
    H = purrr::map2(
      .x = gsea, 
      .y = des_color,
      .f = function(.x, .y) {
        .x %>% 
          dplyr::filter(Cat == "H") %>% 
          dplyr::filter(p.adjust < 0.05, abs(NES) > 1.5) %>% 
          dplyr::select(des = Description, core_enrichment) ->
          .xx
        
        .xx
        
        .y %>% 
          dplyr::filter(color != "grey") %>% 
          dplyr::select(GeneName, log2FC, FDR, color) ->
          .yy
        
        .xx %>% 
          dplyr::mutate(
            g = purrr::map(
              .x = core_enrichment,
              .f = function(.g) {
                .gg <- strsplit(x = .g, split = "/")[[1]]
                .yy %>% 
                  dplyr::filter(GeneName %in% .gg) %>% 
                  dplyr::arrange(-FDR, abs(log2FC))
              }
            )
          ) %>% 
          dplyr::select(des, g)
      }
    )
  ) %>% 
  dplyr::select(seq, vs, H) %>% 
  tidyr::unnest(cols = H) %>% 
  tidyr::unnest(cols = g) %>% 
  dplyr::arrange(-FDR) ->
  hallmark_deg

readr::write_rds(
  x = hallmark_deg,
  file = file.path(
    "analysis/2022-07-20/rda",
    "hallmark_deg.se.rds.gz"
  )
)

hallmark_deg %>% 
  dplyr::group_by(GeneName) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(-n) ->
  hallmark_deg_arrange
  
hallmark_deg %>% 
  dplyr::filter(seq == "Brain") %>% 
  dplyr::group_by(GeneName) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(-n) ->
  hallmark_deg_arrange_brain


hallmark_deg %>% 
  dplyr::filter(seq == "Skull") %>% 
  dplyr::group_by(GeneName) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(-n) ->
  hallmark_deg_arrange_skull

hallmark_deg %>% 
  dplyr::select(des, GeneName) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(des) %>% 
  dplyr::count() %>% 
  dplyr::arrange(-n) %>% 
  print(n = Inf)

hallmark_deg %>% 
  dplyr::select(des, GeneName) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(des == "HALLMARK_P53_PATHWAY") %>% 
  dplyr::pull(GeneName)->
  gene_set_p53

hallmark_deg %>% 
  dplyr::select(des, GeneName) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(des == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  dplyr::pull(GeneName) ->
  gene_set_gama


# gene_set_cuproptosis <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7B") %>% stringr::str_to_sentence()
# gene_set_ferroptosis <- msig_df_h_ferroptosis$gene_symbol
# gene_set_pyroptosis <- msig_df_h_pyroptosis$gene_symbol
# gene_set_apoptosis <- setdiff(setdiff(msig_df_h_apoptosis$gene_symbol, gene_set_ferroptosis), gene_set_pyroptosis)

gene_set_cuproptosis <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7B") %>% stringr::str_to_sentence()
gene_set_ferroptosis <- msig_df_h_ferroptosis$gene_symbol
gene_set_necroptosis <- msig_df_h_necroptosis$gene_symbol
gene_set_pyroptosis <- msig_df_h_pyroptosis$gene_symbol
gene_set_apoptosis <- setdiff(setdiff(msig_df_h_apoptosis$gene_symbol, gene_set_ferroptosis), gene_set_pyroptosis)


g_color <- c("black", RColorBrewer::brewer.pal(n = 7, name = "Set2"))
g_pathway <- c("others", "P53 Pathway", "INTERFERON GAMMA RESPONSE", "Cuproptosis Pathway", "Ferroptosis Pathway", "Necroptosis Pathway", "Pyroptosis Pathway", "Apoptosis Pathway")
names(g_color) <- g_pathway

hallmark_deg_arrange_brain %>% 
  dplyr::filter(n > 15) %>% 
  dplyr::select(g = GeneName) %>%  
  dplyr::mutate(t = ifelse(g %in% gene_set_p53, "P53 Pathway", "others")) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_gama, "INTERFERON GAMMA RESPONSE", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_cuproptosis, "Cuproptosis Pathway", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_ferroptosis, "Ferroptosis Pathway", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_necroptosis, "Necroptosis Pathway", t)) %>%
  dplyr::mutate(t = ifelse(g %in% gene_set_pyroptosis, "Pyroptosis Pathway", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_apoptosis, "Apoptosis Pathway", t)) %>% 
  dplyr::mutate(color = plyr::revalue(x = t, g_color)) ->
  hallmark_deg_arrange_brain_g

hallmark_deg_arrange_skull %>% 
  dplyr::filter(n > 15) %>% 
  dplyr::select(g = GeneName) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_p53, "P53 Pathway", "others")) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_gama, "INTERFERON GAMMA RESPONSE", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_cuproptosis, "Cuproptosis Pathway", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_ferroptosis, "Ferroptosis Pathway", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_necroptosis, "Necroptosis Pathway", t)) %>%
  dplyr::mutate(t = ifelse(g %in% gene_set_pyroptosis, "Pyroptosis Pathway", t)) %>% 
  dplyr::mutate(t = ifelse(g %in% gene_set_apoptosis, "Apoptosis Pathway", t)) %>% 
  dplyr::mutate(color = plyr::revalue(x = t, g_color)) ->
  hallmark_deg_arrange_skull_g

g_color %>% 
  tibble::enframe() %>% 
  dplyr::mutate(name = stringr::str_to_sentence(name)) %>% 
  dplyr::arrange(name) %>% 
  dplyr::mutate(name = factor(name, levels = name)) ->
  g_color_legend

g_color_legend %>% 
  ggplot(aes(x = 1, y = name, fill = name)) +
  geom_tile() +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = g_color_legend$value) +
  coord_fixed() +
  theme(
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    legend.position = "none",
  ) ->
  g_color_legend_plot

ggsave(
  filename = "heatmap_text_legend.pdf",
  plot = g_color_legend_plot,
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/heatmap",
  width = 3,
  height = 7
)



# Mark heatmap ------------------------------------------------------------

{
  pdf(file = "analysis/2022-07-20/result/integration/heatmap/Brain_heatmap_mark.pdf", width = 8, height = 10)
  ComplexHeatmap::draw(
    object = fn_heatmap_se_g(
      .se = se_brain_deg, 
      .g = hallmark_deg_arrange_brain_g
    )
  )
  dev.off()
  
  pdf(file = "analysis/2022-07-20/result/integration/heatmap/Skull_heatmap_mark.pdf", width = 8, height = 10)
  ComplexHeatmap::draw(
    object = fn_heatmap_se_g(
      .se = se_skull_deg, 
      .g = hallmark_deg_arrange_skull_g
      )
    )
  dev.off()
}

# Cell death heatmap --------------------------------------------------------------
gene_set_cuproptosis <- c("FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7B") %>% stringr::str_to_sentence()
gene_set_ferroptosis <- msig_df_h_ferroptosis$gene_symbol
gene_set_necroptosis <- msig_df_h_necroptosis$gene_symbol
gene_set_pyroptosis <- msig_df_h_pyroptosis$gene_symbol
gene_set_apoptosis <- setdiff(setdiff(msig_df_h_apoptosis$gene_symbol, gene_set_ferroptosis), gene_set_pyroptosis)

list(
  cuproptosis = gene_set_cuproptosis, 
  ferroptosis = gene_set_ferroptosis, 
  necroptosis = gene_set_necroptosis,
  pyroptosis = gene_set_pyroptosis, 
  apoptosis = gene_set_apoptosis
) %>% 
  tibble::enframe(name = "celldeath", value = "GeneName") %>%
  dplyr::mutate(celldeath = stringr::str_to_sentence(string = celldeath)) %>% 
  dplyr::mutate(color = g_color[c(4:8)]) %>% 
  tidyr::unnest(cols = GeneName) %>% 
  dplyr::filter(GeneName %in% rownames(se)) ->
  gene_set_celldeath

readr::write_rds(
  x = gene_set_celldeath,
  file = file.path(
    "analysis/2022-07-20/rda",
    "gene_set_celldeath.rds.gz"
  )
)


{
  pdf(
    file = "analysis/2022-07-20/result/integration/heatmap/celldeath_heatmap.pdf", 
    width = 18, 
    height = 6
  )
  ComplexHeatmap::draw(
    object = fn_heatmap_se_celldeath(
      .se = se[gene_set_celldeath$GeneName], 
      .g = gene_set_celldeath
    ),
      merge_legend = TRUE,
      heatmap_legend_side  = "bottom",
      annotation_legend_side  = "bottom"
  )
  dev.off()
  
}




# Cell death gene example -------------------------------------------------

gene_set_celldeath %>% 
  dplyr::mutate(
    p = purrr::map2(
      .x = GeneName,
      .y = celldeath,
      .f = function(.x, .y) {
        
        c("Brain", "Skull") %>% 
        purrr::map(
          .f = function(.z) {
            print(.z)
            .comparisons <- if(.z == "Brain") {
              list(
                c("BC", "B2B"),
                c("BC", "B6B"),
                c("BC", "S2B"),
                c("BC", "S4B")
              )
            } else {
              list(
                c("SC", "B2S"),
                c("SC", "B6S"),
                c("SC", "S2S"),
                c("SC", "S4S")
              )
            }
            
            assay(se[, se$seq == .z])[.x, ] %>% 
              tibble::enframe(name = "barcode", value = "expr") %>% 
              dplyr::inner_join(se_coldata, by = "barcode") %>% 
              dplyr::left_join(colors, by = "group") %>% 
              ggpubr::ggboxplot(
                x = "group", 
                y = "expr", 
                color = "black",
                fill = "group", 
                palette = Ripk1 %>% 
                  dplyr::select(group, color) %>% 
                  dplyr::distinct() %>% 
                  dplyr::pull(color)
              ) +
              ggpubr::stat_compare_means(
                comparisons = .comparisons,
                label = "p.signif",
                method = "t.test",
              ) +
              scale_y_log10()+
              theme(legend.position = "none") +
              labs(
                x = "Group",
                y = "Normalized count",
                title = "{.z} - {.y} - {.x}" %>% glue::glue()
              ) ->
              .p
            ggsave(
              filename = "{.z}_{.y}_{.x}.pdf" %>% glue::glue(),
              plot = .p,
              device = "pdf",
              path = "analysis/2022-07-20/result/integration/heatmap/celldeathboxplot",
              width = 5,
              height = 4
            )
          }
        )
        
      }
    )
  )

# Ripk1 # for brain

# Gls # for skull


# DEG HALLMARK ------------------------------------------------------------


# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/08-heatmap.rda")
load(file = "analysis/2022-07-20/rda/08-heatmap.rda")