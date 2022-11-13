# Metainfo ----------------------------------------------------------------

# @DATE: Tue Sep  6 18:36:33 2022
# @DESCRIPTION: tcga 5op 50 genes

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# load data ---------------------------------------------------------------


tcga_deg <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/DifferentialExpressionTable.xlsx")

tcga_gsva <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/DifferentialGSVATable.xlsx")
tcga_gsva_immune <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/GsvaImmuTable.xlsx")
tcga_gsva_survival <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/GsvaSurvivalTable.xlsx")

tcga_gsea <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/GseaTable.xlsx")

tcga_expr_survival <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/ExpressionAndSurvivalTable.xlsx")

tcga_expr_immune <- readxl::read_excel(path = "analysis/2022-07-20/result/integration/gsca/ImmuneAndExprTable.xlsx")


# DEG ---------------------------------------------------------------------


tcga_deg %>% 
  dplyr::filter(fdr < 0.05) %>% 
  dplyr::mutate(log2fc = log2(fc)) %>% 
  dplyr::mutate(FDR = -log10(fdr)) %>% 
  dplyr::mutate(FDR = ifelse(FDR > 15, 15, FDR)) %>% 
  dplyr::mutate(log2fc = ifelse(log2fc > 3, 3, log2fc)) %>% 
  dplyr::mutate(log2fc = ifelse(log2fc < -3, -3, log2fc)) ->
  tcga_deg_bubble

tcga_deg_bubble %>% 
  dplyr::mutate(
    expr_pattern = purrr::map2_dbl(
      .x = fc,
      .y = fdr,
      .f = function(fc, fdr) {
        if ((fc > 3 / 2) && (fdr < 0.05)) {
          return(1)
        } else if ((fc < 2 / 3) && (fdr < 0.05)) {
          return(-1)
        } else {
          return(0)
        }
      }
    )
  ) %>% 
  dplyr::select(cancertype, symbol, expr_pattern) %>% 
  tidyr::spread(key = cancertype, value = expr_pattern) %>%
  dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) ->
  tcga_deg_bubble_pattern

tcga_deg_bubble_pattern %>% 
  dplyr::rowwise() %>%
  dplyr::do(
    symbol = .$symbol,
    rank =  unlist(.[-1], use.names = F) %>% sum(),
    up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
    down = (unlist(.[-1], use.names = F) == -1) %>% sum()
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::mutate(up_p = up / 14, down_p = down / 14, none = 1 - up_p - down_p) %>% 
  dplyr::arrange(-rank) ->
  gene_rank



tcga_deg_bubble_pattern %>% 
  dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(abs(.)))) %>%
  tidyr::gather(key = cancer_types, value = rank) %>%
  dplyr::arrange(dplyr::desc(rank)) ->
  cancer_rank

tcga_deg_bubble %>% 
  ggplot(aes(y = cancertype, x = symbol)) +
  geom_point(aes(size = FDR, col = log2fc)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-3, 3, length.out = 5),
    labels = c("<= -3", "-1.5", "0", "1.5", ">= 3"),
    name = "log2 FC"
  ) +
  scale_size_continuous(
    limit = c(-log10(0.05), 15),
    range = c(1, 6),
    breaks = c(-log10(0.05), 5, 10, 15),
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    name = "FDR"
  ) +
  scale_x_discrete(limit = gene_rank$symbol) +
  scale_y_discrete(limit = cancer_rank$cancer_types) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(color = "black", size = 12)
  ) ->
  tcga_deg_bubble_plot;tcga_deg_bubble_plot

ggsave(
  plot = tcga_deg_bubble_plot,
  filename = "tcga-expr-bubble-plot.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/tcga-result/",
  width = 11, 
  height = 4
)


gene_rank %>% 
  dplyr::select(symbol, up, down) %>% 
  tidyr::gather(key = type, value = n, - symbol) %>% 
  ggplot(aes(x = symbol, y = n, fill = type)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_x_discrete(limit = gene_rank$symbol) +
  scale_fill_manual(
    name = NULL,
    values = c("blue", "red"),
    label = c("Downregulation", "Upregulation")
  ) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  theme(
  panel.background = element_rect(colour = "black", fill = "white"),
  panel.grid = element_line(colour = "grey", linetype = "dashed"),
  panel.grid.major = element_line(
    colour = "grey",
    linetype = "dashed",
    size = 0.2
  ),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  legend.direction = "horizontal",
  legend.position = c(0.5, 0.8),
  legend.key.size = unit(2, 'mm')
) +
  labs(y = "# of cancer types") ->
  gene_rank_barplot;gene_rank_barplot
ggsave(
  plot = gene_rank_barplot,
  filename = "tcga-expr-bubble-plot-barplot.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/tcga-result/",
  width = 10, 
  height = 2
)

# GSVA --------------------------------------------------------------------

tcga_gsva %>% 
  dplyr::group_by(cancertype) %>% 
  tidyr::nest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(m = purrr::map(
    .x = data,
    .f = function(.x) {
      .m <- median(.x$gsva)
      .mark <- tryCatch(
        expr = {
          .t <- t.test(gsva ~ type, data = .x)
          if(.t$p.value > 0.05) {
            ""
          } else if(.t$p.value > 0.01) {
            "*"
          } else {
            "**"
          } 
        },
        error = function(e) {
          ""
        }
      )
      tibble::tibble(
        median = .m,
        mark = .mark
      )
    }
  )) %>% 
  tidyr::unnest(cols = m) %>% 
  dplyr::select(-data) %>% 
  dplyr::arrange(-median) %>% 
  dplyr::mutate(y = median + 1) ->
  tcga_gsva_rank

tcga_gsva %>% 
  dplyr::mutate(cancertype = factor(cancertype, levels = tcga_gsva_rank$cancertype)) %>% 
  ggplot(aes(y = cancertype, x = gsva, color = type)) +
  geom_boxplot(outlier.color = NA) +
  ggsci::scale_color_aaas(
    name = "Group",
    limit = c("normal", "tumor"),
    label = c("Noraml", "Tumor")
  ) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(color = "black", size = 12),
    legend.position = "top"
  ) +
  labs(
    y = "Cancer types",
    x = "GSVA Score"
  ) +
  annotate(
    geom = "text",
    y = tcga_gsva_rank$cancertype,
    x = tcga_gsva_rank$y,
    label = tcga_gsva_rank$mark
  ) ->
  gsva_boxplot;gsva_boxplot

ggsave(
  plot = gsva_boxplot,
  filename = "tcga-gsva-boxplot.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/tcga-result/",
  width = 3.5, 
  height = 8
)




# GSVA survival -----------------------------------------------------------


tcga_gsva_survival %>%
  dplyr::mutate(group=ifelse(coxp_categorical>0.05,">0.05","<=0.05")) %>%
  dplyr::rename(HR=hr_categorical) %>%
  dplyr::mutate(HR = ifelse(HR > 2.5, 2.5, HR)) %>% 
  dplyr::mutate(logp = -log10(coxp_categorical))-> 
  for_plot


fn_filter_pattern <- function(trend,value,trend1,trend2,p_cutoff) {
  if(!is.na(trend) && !is.na(value)){
    if ((trend == trend1) && (value < p_cutoff)) {
      return(1)
    } else if ((trend == trend2) && (value <p_cutoff)) {
      return(-1)
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}

fn_get_pattern <- function(
    .x,trend1,
    trend2,
    p_cutoff,
    selections) {
  .x %>%
    dplyr::mutate(pattern = purrr::map2_dbl(trend,value, fn_filter_pattern,trend1=trend1,trend2=trend2,p_cutoff=p_cutoff)) %>%
    dplyr::select(all_of(selections), pattern ) %>%
    unique() %>%
    tidyr::spread(key = cancertype, value = pattern) %>%
    dplyr::mutate_if(.predicate = is.numeric, .funs = function(.) {ifelse(is.na(.), 0, .)})
}


fn_get_cancer_types_rank_v2 <- function(.x) {
  .x %>%
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum((.)))) %>%
    tidyr::gather(key = cancertype, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
}

fetched_data_clean_pattern <- fn_get_pattern(
  .x = for_plot %>% 
    dplyr::rename(
      value=coxp_categorical,
      trend=higher_risk_of_death
      ) %>% 
    dplyr::filter(sur_type=="OS"),
  trend1="Higher GSVA",
  trend2="Lower GSVA",
  p_cutoff=0.05,
  selections =c("cancertype")
  )

cancer_rank <- fn_get_cancer_types_rank_v2(
  .x = fetched_data_clean_pattern
  )

bubble_plot <- function(data, cancer, gene, size, fill, fillmipoint=0, fillbreaks,colorgroup, ylab="Symbol",xlab="Cancer types",facet_exp,cancer_rank, gene_rank, sizename,colorvalue, colorbreaks, colorname, fillname, title) {
  CPCOLS <- c("red", "white", "blue")
  data %>%
    ggplot(aes_string(y = gene, x = cancer)) +
    geom_point(
      aes_string(
        size = size, 
        color = fill, 
        ), 
    ) +
    scale_y_discrete(limit = gene_rank) +
    scale_x_discrete(limit = cancer_rank) +
    labs(
      x = xlab,
      y = ylab,
      title = ""
    ) +
    scale_size_continuous(
      name = sizename,#  "-Log10(FDR)"
      breaks = c(0,1.3,2,3,4),
      labels = c("1","0.05","0.01","0.001","<=0.0001")
      #guide=FALSE
    ) +
    scale_color_gradient2(
      name = fillname, # "Methylation diff (T - N)",
      low = CPCOLS[3],
      mid = CPCOLS[2],
      high = CPCOLS[1],
      midpoint = fillmipoint,
      limits=c(min(fillbreaks),max(fillbreaks)),
      breaks=fillbreaks
    ) +
    guides(color=guide_colourbar(title.position="top",reverse=TRUE)) +
    # scale_color_manual(
    #   values = colorvalue, #c("black","grey"),
    #   breaks = colorbreaks, #c("FDR<0.05","FDR>0.05"),
    #   name=colorname
    # ) +
    # guides(color=guide_legend(title.position="top")) 
    theme(
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      axis.text.y = element_text(size = 12,colour = "black"),
      axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 12,colour = "black"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size=12),
      legend.title = element_text(size = 12),
      legend.key = element_rect(fill = "white", colour = "white"),
      legend.key.size = unit(0.5, "cm")
    ) -> p
  if(!is.na(facet_exp)){
    p +
      facet_grid(as.formula(facet_exp)) ->p
  } else {
    p -> p
  }
  return(p)
}



for_plot %>%
  dplyr::filter(!is.na(HR)) %>%
  dplyr::mutate(HR=ifelse(HR>=2.5, 2.5,HR)) %>%
  .$HR -> HR_value
min(HR_value) %>% floor() -> min
max(HR_value) %>% ceiling() -> max
fillbreaks <- sort(unique(c(1,min,max,seq(min,max,length.out = 3))))
title <- "Survival between high and low GSVA score group"
gsva_survival_bubble <- bubble_plot(data=for_plot%>%
              dplyr::filter(HR<10), 
            cancer="sur_type",
            gene="cancertype", 
            ylab="Cancer type", 
            xlab="", 
            facet_exp = NA,
            size="logp", 
            fill="HR", 
            fillmipoint =1,
            fillbreaks =fillbreaks,
            colorgroup="group",
            gene_rank=cancer_rank$cancertype, 
            cancer_rank=c("OS","PFS","DSS","DFI"),
            sizename= "Cox P", 
            fillname="Hazard ratio", 
            colorvalue=c("black","grey"), 
            colorbreaks=c("<=0.05",">0.05"),
            colorname="Cox P value",
            title=title)


ggsave(
  plot = gsva_survival_bubble,
  filename = "tcga-gsva-survival-bubble.pdf",
  device = "pdf",
  path = "analysis/2022-07-20/result/integration/tcga-result/",
  width = 4, 
  height = 8
)

# save --------------------------------------------------------------------

save.image(file = "analysis/2022-07-20/rda/10-tcga.rda")
