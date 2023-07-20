#' Function to make boxplot of median/mean expression of a gene per indication
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#'
#'
#' @param df VIBE harmonized dataframe
#' @param gene Official gene symbol
#' @param method "median" or "mean"; method to calculate the summarized expression for the gene
#' @param order "T" or "F"; ordering x-axis of the plot in decreasing order
#' @param dataset character string with name of dataset, used in caption
#' @param plot_groups Character vector with values you wish to include in the plot, however threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#' @param grouping_var groups object based on this variable; intr if not provided
#'
#' @return ggplot barplot with expression median/mean for the gene
#'
#' @examples
#' # prepare data
#'data("VIBE_data")
#'df_harmonized = harmonize_df(df = VIBE_data,
#'                             col_patientid = "patient_no",
#'                           col_sampleid = "analysis",
#'                           col_indication = "tumor",
#'                           col_treatment = "treatment_flag",
#'                           col_symbol = "gene",
#'                           col_value = "log2_tpm",
#'                           unit = "tpm",
#'                           keep_extra_cols = c("database"))
#'
#' # example1: plot box for each unique intr value
#'plot_expression_box(df = df_harmonized, gene = "Tumor target", method = "median", order = "F",
#'                    dataset = "VIBE dummy data (not real data!)")
#' # example2: order boxes based on percentage above threshold
#'plot_expression_box(df = df_harmonized, gene = "Tumor target", method = "median", order = "T",
#'                    dataset = "VIBE dummy data (not real data!)")
#' # example3: group data by indication & plot selection of indications
#'plot_expression_box(df = df_harmonized, gene = "Tumor target", method = "median", order = "T",
#'dataset = "VIBE dummy data (not real data!)",
#'grouping_var = "indication", plot_groups = c("RCC", "NSCLC", "PDAC"))

plot_expression_box <- function (df, gene, method = c("median", "mean"), order = c("T", "F"), dataset =NA, plot_groups = NA, grouping_var = NA) {

  ifelse(!is.na(dataset), dataset <- dataset, dataset <- "Dataset")
  title = paste(gene, "expression")
  caption = paste0("Dataset: ", dataset, "  \n", "Each dot indicates a single patient sample\n",
                   "The number of patient samples for each indication is depicted below the bars\n",
                   "Percentage indicates the percentage of patients with ",
                   gene, " expression above the cut-off level\n", "Cut-off:  gene expression over all samples in dataset is indicated with dashed line")
  ifelse(!is.na(grouping_var), grouping_var <- grouping_var,
         grouping_var <- "intr")
  threshold = get_threshold(df, gene, method)[1, 1]
  nr = get_numbers(df, grouping_var = grouping_var)
  perc = get_perc_gene(df, gene, method, grouping_var) %>%
    rename(perc_high = 2)
  unit2 = df$unit[1]
  info = perc %>% full_join(nr, by = grouping_var) %>% mutate(perc_high = round(perc_high)) %>%
    arrange(-perc_high)
  df = df %>% dplyr::filter(symbol == gene)
  label_pos = ceiling(max(df$value))
  n_pos = ceiling(min(df$value))

  if (!is.na(plot_groups[1])) {
    df = dplyr::filter(df, !!sym(grouping_var) %in% plot_groups)
    info = dplyr::filter(info, !!sym(grouping_var) %in% plot_groups)
    n_pos = ceiling(min(df$value))
  }

  if (order == "T") {
    df[[grouping_var]] = factor(df[[grouping_var]], levels = info[[grouping_var]])
    plot = ggplot(df, aes(y = value, x = !!sym(grouping_var))) +
      geom_boxplot(fill = gmb_palette[["darkgreen"]], outlier.shape = NA) +
      geom_jitter(color = "black", size = 0.1) + xlab("") +
      ylab(paste0("Expression log2(", unit2, ")")) + geom_hline(yintercept = threshold,
                                                               col = gmb_palette[["orange"]]) + theme_classic2() + theme(axis.text.x = element_text(angle = 45,
                                                                                                                                           hjust = 1)) + labs(title = title, caption = caption) +
      expand_limits(y = c(n_pos, label_pos)) + geom_text(data = info,
                                               aes(label = paste0("n=", nr)), size = 2, y = n_pos) +
      geom_label(data = info, aes(label = paste0(perc_high,
                                                 "%")), size = 2, y = label_pos)
  }
  if (order == "F") {
    plot = ggplot(df, aes(y = value, x = !!sym(grouping_var))) +
      geom_boxplot(fill = gmb_palette[["darkgreen"]], outlier.shape = NA) +
      geom_jitter(color = "black", size = 0.1) + xlab("") +
      ylab(paste0("Expression log(", unit2, ")")) + geom_hline(yintercept = threshold,
                                                               col = gmb_palette[["orange"]]) + theme_classic2() + theme(axis.text.x = element_text(angle = 45,
                                                                                                                                           hjust = 1)) + labs(title = title, caption = caption) +
      expand_limits(y = c(n_pos, label_pos)) + geom_text(data = info,
                                               aes(label = paste0("n=", nr)), size = 2, y = n_pos) +
      geom_label(data = info, aes(label = paste0(perc_high,
                                                 "%")), size = 2, y = label_pos)
  }
  return(plot)
}



