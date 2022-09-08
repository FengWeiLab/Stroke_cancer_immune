# @DESCRIPTION: de

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(DESeq2)

# Load data ---------------------------------------------------------------

se <- readr::read_rds(
  file = file.path(
    "analysis/2022-07-20/rda",
    "count_matrix.se.rds.gz"
  )
)


# seq brain
bb <- c("BC", "B2B", "B6B")

bcb2b <- c("BC", "B2B")
bcb6b <- c("BC", "B6B")
b2b6b <- c("B2B", "B6B")

sb <- c("BC", "S2B", "S4B")

bcs2b <- c("BC", "S2B")
bcs4b <- c("BC", "S4B")
s2s4b <- c("S2B", "S4B")

b2s2b <- c("B2B", "S2B")
b6s4b <- c("B6B", "S4B")


# seq skull

ss <- c("SC", "S2S", "S4S")

scs2s <- c("SC", "S2S")
scs4s <- c("SC", "S4S")
s2s4s <- c("S2S", "S4S")

bs <- c("SC", "B2S", "B6S")

scb2s <- c("SC", "B2S")
scb6s <- c("SC", "B6S")
b2b6s <- c("B2S", "B6S")

s2b2s <- c("S2S", "B2S")
s4b6s <- c("S4S", "B6S")

colData(se) %>% as.data.frame() -> se_coldata

# list(bb = bb, ss = ss, bs = bs, sb = sb) %>% 
list(
  # bb
  bcb2b = bcb2b,
  bcb6b = bcb6b,
  b2b6b = b2b6b,
  
  # sb
  bcs2b = bcs2b,
  bcs4b = bcs4b,
  s2s4b = s2s4b,
  
  b2s2b = b2s2b,
  b6s4b = b6s4b,
  
  # ss
  scs2s = scs2s,
  scs4s = scs4s,
  s2s4s = s2s4s,
  
  # bs
  scb2s = scb2s,
  scb6s = scb6s,
  b2b6s = b2b6s,
  
  s2b2s = s2b2s,
  s4b6s = s4b6s
) %>% 
  tibble::enframe(name = "bs", value = "group") %>% 
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "b$", x = bs), "Brain", "Skull")
  ) %>% 
  dplyr::mutate(
    se = purrr::map(
      .x = group,
      .f = function(.x) {
        se_coldata %>% 
          dplyr::filter(group %in% .x) ->
          .xx
        
        .xxx <- assay(se)[, .xx$barcode]
        .se <- SummarizedExperiment::SummarizedExperiment(
          assays = .xxx,
          colData = .xx
        )
      }
    )
  ) ->
  se_group

readr::write_rds(
  x = se_group,
  file = file.path(
    "analysis/2022-07-20/rda",
    "count_matrix.se.gruop.rds.gz"
  )
)

# Function ----------------------------------------------------------------

fn_se_de <- function(.x, .y) {
  # .x <- se_group$group[[1]]
  # .y <- se_group$se[[1]]
  
  .x1 <- .x[c(1, 2)]
  # .x2 <- .x[c(1, 3)]
  # .x3 <- .x[c(2, 3)]
  
  .y1 <- .y[,.y@colData[.y$group %in% .x1, ]$barcode]
  .y1$group <- factor(.y1$group, .x1)
  .y1_dds <- DESeqDataSet(.y1, design = ~ group)
  .y1_res <- results(DESeq(.y1_dds))
  
  # .y2 <- .y[,.y@colData[.y$group %in% .x2, ]$barcode]
  # .y2$group <- factor(.y2$group, .x2)
  # .y2_dds <- DESeqDataSet(.y2, design = ~ group)
  # .y2_res <- results(DESeq(.y2_dds))
  # 
  # .y3 <- .y[,.y@colData[.y$group %in% .x3, ]$barcode]
  # .y3$group <- factor(.y3$group, .x3)
  # .y3_dds <- DESeqDataSet(.y3, design = ~ group)
  # .y3_res <- results(DESeq(.y3_dds))
  
  .n1 <- paste0(rev(.x1), collapse = "_vs_")
  # .n2 <- paste0(rev(.x2), collapse = "_vs_")
  # .n3 <- paste0(rev(.x3), collapse = "_vs_")
  
  tibble::tibble(
    n1 = list(.y1_res),
    # n2 = list(.y2_res),
    # n3 = list(.y3_res)
  ) ->
    .d
  # colnames(.d) <- c(.n1, .n2, .n3)
  colnames(.d) <- c(.n1)
  
  .d
}

fn_volcano <- function(.x, .y) {
  # .x <- se_group_de$seq[[2]]
  # .y <- se_group_de$de[[2]]
  .prefix <- .x
  
  .y %>% 
    tidyr::gather(key = "vs", value = "des") ->
    .yy
  
  .yy %>% 
    dplyr::mutate(
      des_color = purrr::map(
        .x = des,
        .f = function(.des) {
          .des %>% 
            as.data.frame() %>%
            tibble::rownames_to_column(var = "GeneName") %>% 
            dplyr::filter(!is.na(!padj)) %>% 
            dplyr::filter(!is.na(log2FoldChange)) %>% 
            # dplyr::mutate(log2FC = ifelse(log2FoldChange > 4, 4, log2FoldChange)) %>% 
            # dplyr::mutate(log2FC = ifelse(log2FC < -4, -4, log2FC)) %>% 
            dplyr::mutate(log2FC = log2FoldChange) %>% 
            dplyr::mutate(FDR = -log10(padj)) %>% 
            dplyr::mutate(
              color = dplyr::case_when(
                log2FC > log2(2) & padj < 0.05 ~ "red",
                log2FC < log2(1/2) & padj < 0.05 ~ "green",
                TRUE ~ "grey"
              )
            )
        }
      )
    ) ->
    .yyy
  
  .yyy %>% 
    dplyr::mutate(
      vp = purrr::map2(
        .x = vs,
        .y = des_color,
        .f = function(.x, .y) {
          # .x <- .yyy$vs[[1]]
          # .y <- .yyy$des_color[[1]]
          .xd <- .y
          
          .xd %>% 
            dplyr::filter(color != "grey") %>% 
            dplyr::group_by(color) %>% 
            dplyr::count() %>% 
            dplyr::ungroup() %>% 
            tibble::deframe() ->
            .xxx
          
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
              # size = 6
            ) +
            ggrepel::geom_text_repel(
              aes(label = GeneName),
              data = subset(.xd, color == "green") %>% 
                dplyr::arrange(-FDR, -abs(log2FC)) %>% 
                dplyr::slice(1:10),
              box.padding = 0.5,
              max.overlaps = Inf,
              # size = 6
            ) +
            scale_x_continuous(
              # limits = c(-4, 4),
              expand = c(0.02, 0)
            ) +
            scale_y_continuous(
              expand = c(0.01, 0),
              limits = c(
                0,
                ceiling(
                  max(
                    .xd %>% 
                      dplyr::filter(!is.infinite(FDR)) %>% 
                      dplyr::pull(FDR)
                  ) / 10
                ) * 10
              )
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
              y = "-log10(FDR)",
              title = glue::glue("{.x}; Up (n={.xxx[2]}), Down (n={.xxx[1]})")
            ) ->
            .p
          .plotfilename <- glue::glue("volcano_plot_{.prefix}_{.x}.pdf")
          ggsave(
            filename = .plotfilename,
            plot = .p,
            device = "pdf",
            path = "analysis/2022-07-20/result/integration/de",
            width = 7,
            height = 5
          )
          .p
        }
      )
    )
}

# DE ----------------------------------------------------------------------

se_group %>% 
  dplyr::mutate(
    de = purrr::map2(
      .x = group,
      .y = se,
      .f = fn_se_de
    )
  ) ->
  se_group_de

readr::write_rds(
  x = se_group_de,
  file = file.path(
    "analysis/2022-07-20/rda",
    "se_group_de.rds.gz"
  )
)



unlist(se_group_de$de) %>% 
  purrr::map(.f = function(.x) {
    .x %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "GeneName") %>% 
      dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1) %>% 
      dplyr::arrange(padj, -abs(log2FoldChange))
  }) %>% 
  writexl::write_xlsx(
    path = file.path(
      "analysis/2022-07-20/result/integration/de",
      "01-basic_de.xlsx"
    )
  )


# volcano plot ------------------------------------------------------------

se_group_de %>% 
  dplyr::mutate(
    volcano_plot = purrr::map2(
      .x = seq,
      .y = de,
      .f = fn_volcano
    )
  ) %>% 
  tidyr::unnest(cols = volcano_plot) ->
  se_group_de_volcano

readr::write_rds(
  x = se_group_de_volcano,
  file = file.path(
    "analysis/2022-07-20/rda",
    "se_group_de_volcano.rds.gz"
  )
)


se_group_de_volcano %>% 
  dplyr::select(bs, seq, vp) %>% 
  dplyr::group_by(seq) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    patchwork = purrr::map2(
      .x = seq,
      .y = data,
      .f = function(.x, .y) {
        
        
        .plotfilename <- glue::glue("volcano_plot_{.x}.pdf")
        (.y$vp[[1]] | .y$vp[[2]] | .y$vp[[3]]) / (.y$vp[[4]] | .y$vp[[5]] | .y$vp[[6]]) / (.y$vp[[7]] | .y$vp[[8]] | plot_spacer() ) +
          plot_annotation(
            title = glue::glue("Differential expression {.x}"),  
            tag_levels = "A"
          ) ->
          .p
        
        ggsave(
          filename = .plotfilename,
          plot = .p,
          device = "pdf",
          path = "analysis/2022-07-20/result/integration/de",
          width = 17,
          height = 12
        )
        .p
      }
    )
  ) ->
  se_group_de_volcano_patchworks

# save image --------------------------------------------------------------


save.image(file = "analysis/2022-07-20/rda/06-de.rda")
