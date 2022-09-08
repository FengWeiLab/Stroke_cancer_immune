# @DESCRIPTION: enrichment.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(DESeq2)

# Load data ---------------------------------------------------------------


se_group_de <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "se_group_de_volcano.rds.gz"
  )
) %>% 
  dplyr::select(bs, group, seq, se, vs, des_color)

msig_df <- msigdbr::msigdbr(species = "Mus musculus") 

msig_df %>%
  dplyr::mutate(gs_name = glue::glue("{gs_cat}#{gs_name}")) %>%
  dplyr::select(gs_name, gene_symbol) ->
  msig_df_s
  

S4_S2_SC_color <- c("SC" = "#112a13", "S2" = "#AE1700", "S4" = "#039B9E")

# function ----------------------------------------------------------------

fn_gobp <- function(.de, .color = "red") {
  .cc <- c("red" = "#AE1700", "green" = "#112a13")
  
  .d <- .de %>% dplyr::filter(color == .color)
  
  .go_bp <- clusterProfiler::enrichGO(
    gene = .d$GeneName,
    universe = .de$GeneName,
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
    geom_col(fill = .cc[.color], color = NA, width = 0.7) +
    geom_text(aes(label = Count), hjust = 4, color = "white", size = 5) +
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
    .go_bp_plot
  
  tibble::tibble(
    gobp = list(.go_bp),
    goplot = list(.go_bp_plot)
  )
  
}

fn_gsea <- function(.de) {
  .de %>% 
    dplyr::arrange(-log2FC) %>%
    dplyr::select(GeneName, log2FC) %>%
    tidyr::drop_na() %>%
    tibble::deframe() ->
    geneList
  
  
  .GSEA <- clusterProfiler::GSEA(
    geneList = geneList,
    pAdjustMethod = "BH",
    TERM2GENE = msig_df_s,
    verbose = FALSE
  )
  
  .GSEA %>%
    tibble::as_tibble() %>%
    tidyr::separate(Description, into = c("Cat", "Description"), sep = "#") ->
    .GSEA_sep
  
  .GSEA_sep$Cat %>% table()
    
  .GSEA_sep %>%
    dplyr::filter(Cat == "C5") ->
    .GSEA_sep_C5
  
  .GSEA_sep_C5 %>% 
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, NES) %>%
    dplyr::filter(grepl(pattern = "gobp", Description, ignore.case = TRUE)) %>%
    dplyr::mutate(Description = gsub(pattern = "GOBP_", replacement = "", x = Description)) %>%
    dplyr::mutate(Description = gsub(pattern = "_", replacement = " ", x = Description)) %>%
    dplyr::mutate(Description = stringr::str_wrap(
      stringr::str_to_sentence(string = Description), width = 60
    )) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(color = ifelse(NES > 0, "P", "N")) %>%
    dplyr::slice(1:10, (dplyr::n()-9):dplyr::n()) %>%
    dplyr::distinct() %>% 
    dplyr::arrange(NES, -adjp) %>%
    dplyr::mutate(Description = factor(Description, levels = Description)) ->
    .GSEA_sep_C5_for_plot
  
  .GSEA_sep_C5_for_plot %>%
    ggplot(aes(x = Description, y = NES)) +
    geom_col(aes(fill = color)) +
    scale_fill_manual(values = c("#54AE59", "#9F82B5")) +
    scale_x_discrete(
      limit = .GSEA_sep_C5_for_plot$Description,
      labels =stringr::str_wrap(
        stringr::str_replace_all(
          string = .GSEA_sep_C5_for_plot$Description,
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
    ) ->
    gseaplotc5
  
  .GSEA_sep %>%
    dplyr::filter(Cat == "H") ->
    .GSEA_sep_H
  
  
  .GSEA_sep_H %>% 
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, NES) %>%
    # dplyr::filter(grepl(pattern = "gobp", Description, ignore.case = TRUE)) %>%
    # dplyr::mutate(Description = gsub(pattern = "GOBP_", replacement = "", x = Description)) %>%
    dplyr::mutate(Description = gsub(pattern = "_", replacement = " ", x = Description)) %>%
    dplyr::mutate(Description = stringr::str_wrap(
      Description, width = 60
    )) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(color = ifelse(NES > 0, "P", "N")) ->
    .GSEA_sep_H_d
    
  .GSEA_sep_H_for_plot <- if(nrow(.GSEA_sep_H) > 20) {
    .GSEA_sep_H_d %>% 
      dplyr::slice(1:10, (dplyr::n()-9):dplyr::n()) %>%
      dplyr::distinct() %>% 
      dplyr::arrange(NES, -adjp) %>%
      dplyr::mutate(Description = factor(Description, levels = Description))
  } else {
    .GSEA_sep_H_d %>% 
      dplyr::arrange(NES, -adjp) %>%
      dplyr::mutate(Description = factor(Description, levels = Description))
  }
   
  
  .GSEA_sep_H_for_plot %>%
    ggplot(aes(x = Description, y = NES)) +
    geom_col(aes(fill = color)) +
    scale_fill_manual(values = c("#54AE59", "#9F82B5")) +
    scale_x_discrete(
      limit = .GSEA_sep_H_for_plot$Description,
      labels =stringr::str_wrap(
        stringr::str_replace_all(
          string = .GSEA_sep_H_for_plot$Description,
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
    ) ->
    gseaploth
  
  tibble::tibble(
    gsea = list(.GSEA_sep),
    gseaplotc5 = list(gseaplotc5),
    gseaploth = list(gseaploth)
  )
}

fn_go_enrich <- function(.vs, .de) {
  # .vs <- se_group_de$vs[[12]]
  # .de <- se_group_de$des_color[[12]]
  print(.vs)
  
  .up <- fn_gobp(.de, .color = "red")
  .down <- fn_gobp(.de, .color = "green")
  
  dplyr::bind_rows(
    .up, .down
  ) %>% 
    dplyr::mutate(reg = c("go_up", "go_down")) ->
    .gobp
  
  tryCatch(
    expr = {
      .gobp %>% 
        purrr::pmap(
          .f = function(gobp, goplot, reg) {
            .xlsx_filename <- glue::glue("GOBP-{.vs}-{reg}.xlsx")
            
            writexl::write_xlsx(
              x = as.data.frame(gobp), 
              path = file.path(
                "analysis/2022-07-20/result/integration/go",
                .xlsx_filename
              )
            )
            
            .plot_filename <- glue::glue("GOBP-{.vs}-{reg}.pdf")
            
            ggsave(
              plot = goplot,
              filename = .plot_filename,
              device = "pdf",
              path = "analysis/2022-07-20/result/integration/go",
              width = 10, 
              height = 6.5
            )
          }
        )
    },
    error = function(e) {NULL}
  )

  
  .gsea <- fn_gsea(.de = .de)
  
  tryCatch(
    expr = {
      .gsea %>% 
        purrr::pmap(
          .f = function(gsea, gseaplotc5, gseaploth) {
            .xlsx_filename <- glue::glue("GSEA-MSigDB-{.vs}.xlsx")
            
            writexl::write_xlsx(
              x =  gsea %>% as.data.frame() %>% dplyr::select(-ID, -Cat), 
              path = file.path(
                "analysis/2022-07-20/result/integration/go",
                .xlsx_filename
              )
            )
            
            .plotc5_filename <- glue::glue("GSEA-MSigDB-C5-{.vs}.pdf")
            
            ggsave(
              plot = gseaplotc5,
              filename = .plotc5_filename,
              device = "pdf",
              path = "analysis/2022-07-20/result/integration/go",
              width = 10, 
              height = 6.5
            )
            
            .ploth_filename <- glue::glue("GSEA-MSigDB-H-{.vs}.pdf")
            
            ggsave(
              plot = gseaploth,
              filename = .ploth_filename,
              device = "pdf",
              path = "analysis/2022-07-20/result/integration/go",
              width = 10, 
              height = 6.5
            )
          }
        )
    },
    error = function(e) {NULL}
  )

  
  .gobp %>% 
    dplyr::select(-goplot) %>% 
    tidyr::spread(key = reg, value = gobp) %>% 
    dplyr::bind_cols(.gsea %>% dplyr::select(gsea))
}

# enrichment --------------------------------------------------------------


se_group_de %>% 
  dplyr::mutate(
    enrichment = purrr::map2(
      .x = vs,
      .y = des_color,
      .f = fn_go_enrich
    )
  ) ->
  se_group_de_enrichment

readr::write_rds(
  x = se_group_de_enrichment,
  file = "analysis/2022-07-20/rda/07-enrichment.rds.gz"
)


# save image --------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/06-enrichment.rda")
load(file = "analysis/2022-07-20/rda/06-enrichment.rda")
