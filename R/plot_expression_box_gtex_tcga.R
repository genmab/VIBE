#' Plot a boxplot comparing GTEX and TCGA samples for selected gene with Percentage of samples above median, Kruskal Wallis test, Fold change and sample counts
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggsignif
#' @import ggpubr
#'
#'
#' @param df dataframe harmonized by VIBE's harmonize_df function, with the indication column being a concatenation of the study and primary_site columns of the XENA dataset.
#' @param gene name of genes for creating the plot
#' @param threshold_method the summary stat method to use, either "mean" or "median"
#' @param dataset character string with name of dataset, used in caption
#' @param plot_primary_site Character vector with primary_site values (name of tissue column in Xena dataset) you wish to include in the plot, however threshold will be calculated for target1 and target2 on all data in df.
#'
#' @return boxplot comparing the statistical differences of the sample distribution per gene in GTEX vs TCGA
#'
#'
plot_expression_box_gtex_tcga <- function(df, gene, threshold_method, dataset, plot_primary_site = NA) {

  .Deprecated(new = "plot_expression_box_split")

  # add tissue and dataset columns (the study column is not a mandatory column in harmonize_df function, so put it here to be safe)
  df_plot <- df %>%
    mutate(tissue  = str_remove_all(intr, "GTEX|TCGA") %>% str_trim(),
           dataset = str_extract(intr, "GTEX|TCGA")) %>%
    filter(symbol == gene)

  # extract sample numbers
  df_numbers <- get_numbers(df_plot) %>%
    mutate(tissue  = str_remove_all(intr, "GTEX|TCGA") %>% str_trim(),
           dataset = str_extract(intr, "GTEX|TCGA")) %>%
    arrange(tissue, dataset) %>%
    select(-intr) %>%
    pivot_wider(names_from = dataset, values_from = nr) %>%
    mutate(TCGA = ifelse(is.na(TCGA), 0, TCGA),
           GTEX = ifelse(is.na(GTEX), 0, GTEX),
           nr = paste0("GTEX n=", GTEX, "\n", "TCGA n=", TCGA),
           value = -0.2)

  # get fold change and kruskal-wallis p-value between TCGA and GTEX on tissue level

  df_kw_fc <- calc_fc_gtex_tcga(df = df, threshold_method = "median", return_primary_site = plot_primary_site) %>%
    dplyr::filter(symbol == gene)


  # extract % over threshold
  label_pos = ceiling(max(df_plot$value))

  df_perc <- get_perc_gene(df_plot, gene, threshold_method) %>%
    mutate(tissue  = str_remove_all(intr, "GTEX|TCGA") %>% str_trim(),
           dataset = str_extract(intr, "GTEX|TCGA")) %>%
    arrange(tissue, dataset) %>%
    select(-intr) %>%
    pivot_wider(names_from = dataset, values_from = 1) %>%
    mutate(TCGA = ifelse(is.na(TCGA), 0, TCGA),
           GTEX = ifelse(is.na(GTEX), 0, GTEX),
           perc = paste0("GTEX: ", GTEX, "%\n", "TCGA: ", TCGA, "%"),
           value = label_pos)

  # Combine percentage and statistics together
  df_kw_fc$fc_pval <- paste("FC:",round(df_kw_fc$fold_change,2),"\nPval:", df_kw_fc$p_value)
  df_perc_fc <- left_join(df_perc, df_kw_fc[,c("tissue","fc_pval")], by = "tissue")
  df_perc_fc$perc_stats <- paste(df_perc_fc$perc,"\n",df_perc_fc$fc_pval)


  # calculate threshold
  threshold <- get_threshold(df, gene = gene, method = threshold_method)[[1]]

  # other variables
  dataset_colors = c("TCGA" = palette_light[["green"]], "GTEX" = palette_light[["grey"]])
  unit2 = df$unit[1]


  # Filter data to plot only selected primary_site/tissues, while calculations are done on entire datset
  if(!any(is.na(plot_primary_site))) {
    df_plot <- dplyr::filter(df_plot, tissue %in% plot_primary_site)
    df_numbers <- dplyr::filter(df_numbers, tissue %in% plot_primary_site)
    df_perc_fc <- dplyr::filter(df_perc_fc, tissue %in% plot_primary_site)
  }

  # plot
  plot = ggplot(df_plot, aes(x = tissue, y = value)) +
    geom_boxplot(mapping = aes(fill = dataset)) +
    geom_text(data = df_numbers, aes(label = nr), size = 2) +
    geom_label(data = df_perc_fc, aes(label = perc_stats), size = 2) +
    geom_hline(yintercept = threshold,
               col = gmb_palette[["orange"]], linetype = "dashed") +
    expand_limits(y = label_pos+2) +
    labs(title = paste(gene, "expression in GTEX (normal) and TCGA (cancer) datasets"),
         caption = paste0("Dataset: ", dataset, "\n",
                          "The number of patient samples for each indication is depicted below the bars\n",
                          "Line: ", threshold_method, " gene expression over all samples is indicated with dashed line\n",
                          "FC: ", "Fold change are calculated per tissue using the formula TCGA/GTEX\n",
                          "Pval: ", " P vales are calculated based on the Kruskal Wallis test" )) +
    ylab(paste0("Expression log2(", unit2, ")")) +
    scale_fill_manual(values = dataset_colors) +
    ggpubr::theme_classic2() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  return(plot)

}


