
fn_volcano_plot <- function(.xx) {
  .xx %>% 
    dplyr::mutate(log2FC = ifelse(log2FC > 4, 4, log2FC)) %>% 
    dplyr::mutate(log2FC = ifelse(log2FC < -4, -4, log2FC)) ->
    .xd
  
  .xd %>% 
    ggplot(aes(x = log2FC, y = FDR, color = color)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values =  c("#6BAB62", "grey", "#9A84B2")) +
    geom_segment(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color = "black",
      linetype = 88,
      data = tibble::tibble(
        x1 = c(-Inf, 1, -1, 1),
        y1 = -log10(0.05),
        x2 = c(-1, Inf, -1, 1 ),
        y2 = c(-log10(0.05), -log10(0.05), Inf, Inf))
    ) +
    ggrepel::geom_text_repel(
      aes(label = GeneName),
      data = subset(.xd, color == "red") %>% 
        dplyr::arrange(-FDR, -abs(log2FC)) %>% 
        dplyr::slice(1:10),
      box.padding = 0.5,
      max.overlaps = Inf,
      size = 6
    ) +
    ggrepel::geom_text_repel(
      aes(label = GeneName),
      data = subset(.xd, color == "green") %>% 
        dplyr::arrange(-FDR, -abs(log2FC)) %>% 
        dplyr::slice(1:10),
      box.padding = 0.5,
      max.overlaps = Inf,
      size = 6
    ) +
    scale_x_continuous(
      limits = c(-4, 4),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      expand = c(0.01, 0)
    ) +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      axis.line.x.bottom = element_line(color = "black"),
      axis.line.y.left = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(color = "black", size = 18),
      
      legend.position = "none",
      
    ) +
    labs(
      x = "log2FC",
      y = "-log10(FDR)"
    )
}

fn_selectdata <- function(.x, .y, outpath) {
  # .y <- "UVB0_BC"
  .yy <- strsplit(x = .y, split = "_")[[1]]
  
  
  .x %>% 
    dplyr::select(GeneName, 5:7) %>% 
    dplyr::mutate(
      a = rowMeans(
        dplyr::across(2:4)
      )
    )
  
  .x %>% 
    dplyr::mutate(
      v = rowMeans(dplyr::across(5:7)),
      s = rowMeans(dplyr::across(8:10))
    ) %>% 
    dplyr::filter(!(v < 5 & s < 5)) %>% 
    # dplyr::filter(padj < 0.05, abs(log2FC) > 1) %>% 
    dplyr::mutate(FDR = -log10(padj)) %>% 
    dplyr::mutate(color = dplyr::case_when(
      log2FC > log2(2) & padj < 0.05 ~ "red",
      log2FC < log2(1/2) & padj < 0.05 ~ "green",
      TRUE ~ "grey"
    )) ->
    .xx
  
  .xx %>% 
    dplyr::filter(color != "grey") %>% 
    dplyr::group_by(color) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    tibble::deframe() ->
    .xxx
  
  # volcano
  .plotfilename <- glue::glue("volcano_plot_{.y}.pdf")
  .p <- fn_volcano_plot(.xx) +
    labs(
      title = glue::glue("Volcano plot for {.y}\n Up (n={.xxx[2]}), Down (n={.xxx[1]})")
    )
  
  ggsave(
    filename = .plotfilename,
    plot = .p,
    device = "pdf",
    path = outpath,
    width = 7,
    height = 5
  )
  
  .xx
}

fn_union_genes <- function(.d) {
  .d %>% 
    purrr::map(.f = function(.x) {
      .x %>% 
        dplyr::filter(color != "grey") %>% 
        dplyr::pull(GeneName)
    }) %>% 
    purrr::reduce(.f = union) 
}

fn_combine_data <- function(.d, .name = "UVB0_BC") {
  
  
  .d %>% 
    purrr::map(
      .f = function(.x) {
        .n <- strsplit(.name, split = "_")[[1]]
        
        .colnames <- colnames(.x)
        
        if(!any(grepl(pattern = .n[1], x = .colnames))) {
          return(
            .x %>% 
              dplyr::select(GeneName, 5:10)
          )
        } else if (
          any(grepl(pattern = .n[2], x = .colnames))
        ) {
          return(NULL)
        } else {
          .x %>% 
            dplyr::select(GeneName, dplyr::starts_with(.n[1]))
        }
      }
    ) ->
    .dd
  
  Filter(Negate(is.null), .dd) %>% 
    purrr::reduce(
      .f = dplyr::inner_join,
      by = "GeneName"
    )
  
}

fn_heatmap <- function(.x, .order) {
  
  .col <- colnames(.x)[-1]
  
  .x %>%
    dplyr::distinct(GeneName, .keep_all = TRUE) %>% 
    tibble::column_to_rownames("GeneName") %>%
    dplyr::select(
      which(grepl(pattern = .order[3], .col)),
      which(grepl(pattern = .order[2], .col)),
      which(grepl(pattern = .order[1], .col)),
    ) ->
    .xx
  
  .x_mat <- .xx %>% 
    as.matrix() %>%
    apply(1, scale) %>%
    t()
  
  colnames(.x_mat) <- colnames(.xx)
  
  # gsub(
  #   pattern = "_\\d|D", 
  #   replacement = "", 
  #   x = colnames(x = .x[,-1])
  # ) %>% 
  #   unique() ->
  #   .colname
  .colname <- .order
  
  
  cluster_col <- c("black", "#AE1700","#039B9E")
  names(cluster_col) <- .colname 
  
  library(ComplexHeatmap)
  hma_top = ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(
      Group = gsub(
        pattern="_\\d", 
        replacement = "", 
        x = colnames(.x_mat)
      )
    ),
    gap = unit(c(2,2), "mm"),
    col = list(Group = cluster_col),
    which = "column"
  )
  
  ComplexHeatmap::Heatmap(
    # data and color
    matrix = .x_mat,
    col =  circlize::colorRamp2(c(-1.1, 0, 1.1), c("blue", "white", "red"), space = "RGB"),
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
    column_order = c(7,8,9,4,5,6,1,2,3),
    
    # row labels
    row_labels = rownames(.x_mat),
    row_names_side = 'right',
    show_row_names = T,
    row_names_max_width = unit(6, 'cm'),
    row_names_gp = gpar(fontsize = 6),
    row_names_rot = 0,
    row_names_centered = FALSE,
    
    # column labels
    # column_labels = colnames(BD6_BD2_BC_deg_mat),
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
    row_km = 3,
    row_km_repeats = 10,
    row_split = NULL,
    # column_km = 3,
    # column_km_repeats = 10,
    column_split = rep(c("C", "B", "A"), each = 3),
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

fn_upgenes <- function(.x, .y) {
  # .x <- UVB_filter_upgenes$data[[1]]
  # .y <- UVB_filter_upgenes$upgenes[[1]]
  
  .x %>% 
    dplyr::filter(GeneName %in% .y$GeneName) %>% 
    dplyr::select(-c(2,3,4,5,6,7,8,9,10))
}


fn_gobp <- function(.x, .y, .z, outpath) {
  # .x <- UVB_filter_upgenes_upd$data[[1]]
  # .y <- UVB_filter_upgenes_upd$upgenes[[1]]
  # .z <- UVB_filter_upgenes_upd$up[[1]]
  
  .gobp <- clusterProfiler::enrichGO(
    gene = .y$GeneName,
    universe = .x$GeneName,
    keyType = "SYMBOL",
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
  )
  
  .gobp %>% 
    tibble::as_tibble()  %>% 
    dplyr::mutate(Description = stringr::str_to_sentence(string = Description)) %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, Count) %>%
    head(20) %>%
    dplyr::arrange(adjp, Count) %>%
    # tibble::rowid_to_column() %>%
    dplyr::mutate(Description = factor(Description, levels = Description)) ->
    .gobp_for_plot
  
  .gobp_for_plot %>% 
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = "#9A84B2", color = NA, width = 0.7) +
    geom_text(
      aes(label = Count), 
      hjust = 6, 
      color = "white", 
      size = 5
    ) +
    labs(y = "-log10(Adj. P value)") +
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
    .gobp_plot
  
  ggsave(
    plot = .gobp_plot,
    filename = glue::glue("{.z}-up-GO-BP.pdf"),
    device = "pdf",
    path = outpath,
    width = 7, 
    height = 6.5
  )
  writexl::write_xlsx(
    x = .gobp %>% 
      as.data.frame(), 
    path = file.path(outpath, glue::glue("{.z}-up-GO-BP.xlsx"))
  )
  
  .gs <- .x %>% 
    dplyr::arrange(-log2FC) %>%
    dplyr::select(GeneID, log2FC) %>%
    tidyr::drop_na() %>%
    tibble::deframe()
  
  .z_gsea <- clusterProfiler::GSEA(
    geneList = .gs,
    exponent = 1,
    pAdjustMethod = "BH",
    TERM2GENE = msig_df %>%
      dplyr::mutate(gs_name = glue::glue("{gs_cat}#{gs_name}")) %>%
      dplyr::select(gs_name, entrez_gene),
    verbose = FALSE
  )
  
  .z_gsea %>% 
    tibble::as_tibble() %>%
    tidyr::separate(Description, into = c("Cat", "Description"), sep = "#") ->
    .z_gsea_sep
  
  .z_gsea_sep %>%
    dplyr::filter(Cat == "C5") ->
    .z_gsea_sep_C5
  
  .z_gsea_sep_C5 %>% 
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, NES) %>%
    # dplyr::filter(grepl(pattern = "gobp", Description, ignore.case = TRUE)) %>%
    # dplyr::mutate(Description = gsub(pattern = "GOBP_", replacement = "", x = Description)) %>%
    dplyr::mutate(Description = gsub(pattern = "_", replacement = " ", x = Description)) %>%
    dplyr::mutate(Description = stringr::str_to_lower(Description)) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(color = ifelse(NES > 0, "P", "N")) %>%
    dplyr::slice(1:3, (dplyr::n()-3):dplyr::n()) %>%
    dplyr::arrange(NES, -adjp) %>%
    dplyr::mutate(Description = factor(Description, levels = Description)) ->
    .z_gsea_sep_for_plot
  
  .z_gsea_sep_for_plot %>% 
    ggplot(aes(x = Description, y = NES)) +
    geom_col(aes(fill = color)) +
    scale_fill_manual(values = c("#54AE59", "#9F82B5")) +
    scale_x_discrete(
      limit = .z_gsea_sep_for_plot$Description,
      labels =stringr::str_wrap(
        stringr::str_replace_all(
          string = .z_gsea_sep_for_plot$Description,
          pattern = "_",
          replacement = " "
        ),
        width = 50)
    )+
    labs(y = "Normalized Enrichment Score") +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black"),
      
      axis.title.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      legend.position = "none"
    )  ->
    .z_gsea_plot
  ggsave(
    plot = .z_gsea_plot,
    filename = glue::glue("{.z}-GSEA-MSigDB.pdf"),
    device = "pdf",
    path = outpath,
    width = 7, height = 6
  )
  
  .gobp %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(p.adj = -log10(p.adjust)) %>% 
    dplyr::arrange(-p.adj, -Count) %>% 
    head(3) %>% 
    .$geneID %>% 
    strsplit(split = "/") %>% 
    unlist() %>% 
    sort() %>% 
    unique() %>% 
    head(20)
}