#' Multi-panel plot representing scatterplot per quadrant between two different genes.
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import ggpubr

#'
#'
#' @param df dataframe with expression values for each indication (as obtained from Profiler)
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#' @param target3 Official gene symbol of gene 3. Will be used to color dots by. If left NA dots will be colored by quadrant.
#' @param QOI Which of the four quadrants are of interest for the question; Q1 (Top-left; high-low); Q2 (Top-right; high-high); Q3 (bottom-left; low-high) and Q4 (bottom right; low-low)
#' @param cor_method Correlation method used
#' @param threshold_method Threshold method to use; median or mean or Q3 (75th quantile)
#' @param dataset dataset to be mentioned for captioning
#' @param grouping_var groups object based on this variable; intr if not provided
#' @param plot_groups Character vector with values you wish to include in the plot, threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#' @param facet_ncol number of columns to be plotted
#'
#' @return ggplot2 scatterplot with information of the percentage of samples in the selected quadrant by selected correlation method; barplot with correlation values and the table with correlation values between two genes.
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
#'

#' # example1:
#'plot_expression_scatter(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'                        QOI = "Q2", cor_method = "spearman", threshold_method = "median",
#'                        dataset = "VIBE dummy data (not real data!)")
#' # example2: change grouping_var to indication
#'plot_expression_scatter(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'                        QOI = "Q2", cor_method = "spearman", threshold_method = "median",
#'                        dataset = "VIBE dummy data (not real data!)",
#'                        grouping_var = "indication")
#' # example3: limit plotted indications to selection,
#' # also color datapoints by expression of a third gene
#'plot_expression_scatter(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'                        QOI = "Q2", cor_method = "spearman", threshold_method = "median",
#'                        dataset = "VIBE dummy data (not real data!)",
#'                        grouping_var = "indication", plot_groups = c("BRCA", "NSCLC", "OV", "STAD"),
#'                        target3 = "Gene 9",
#'                        facet_ncol = 2)
#'
plot_expression_scatter <- function (df, target1, target2, target3 = NA, QOI = c("Q1", "Q2", "Q3", "Q4"),
                                     cor_method = c("spearman", "pearson"), threshold_method = c("median", "mean", "Q3"), dataset = NA,
                                     plot_groups = NA, grouping_var = NA, facet_ncol = 4) {

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")
  ifelse(!is.na(dataset), dataset <- dataset, dataset <- "Dataset")

  dir_axis <- direction_and_axis_quadrants(df, QOI, target1,target2)

  quad_summary <- add_quadrant_info(df, target1, target2,
                                    threshold_method, grouping_var)


  quad_perc <- calculate_percentage_per_quadrant(quad_summary,
                                                 grouping_var)
  xintercept_coord <- get_threshold(df, target1, method = threshold_method)[1,1]
  yintercept_coord <- get_threshold(df, target2, method = threshold_method)[1,1]

  caption = paste0("Dataset: ", dataset, "\n", "Cut-off: ",
                   threshold_method, " gene expression over all samples in dataset is indicated with dashed line\n",
                   "Percentage indicates the percentage of patients with ",
                   dir_axis[1], " ", target1, " expression and ", dir_axis[2],
                   " ", target2, " expression\n", "Correlation method: ",
                   cor_method)

  if (!is.na(plot_groups[1])) {
    df <- df %>%
      dplyr::filter(grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), !!sym(grouping_var)))
    quad_summary <- quad_summary %>% dplyr::filter(grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), !!sym(grouping_var)))

    quad_perc = dplyr::select(df, all_of(grouping_var)) %>% distinct() %>%
      right_join(quad_perc,by = grouping_var) %>%
      dplyr::filter(grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), !!sym(grouping_var)))}

  levels_plot <- as.data.frame(table(df[[grouping_var]]))$Var1
  quad_summary[[grouping_var]] <- factor(quad_summary[[grouping_var]],
                                         levels = levels_plot)
  quad_perc[[grouping_var]] <- factor(quad_perc[[grouping_var]],
                                      levels = levels_plot)



  if(is.na(target3)) {
    scatter_quadrant_base <- ggplot(quad_summary, aes(x = !!sym(target1),
                                                 y = !!sym(target2))) +
      geom_point(aes(color = as.factor(quadrant))) +
      scale_color_manual(values = as.vector(gmb_palette)) +
      theme(legend.position = "None") +
      ggtitle(paste0("Scatterplot between ", target1, " and ", target2))
  } else {
    scatter_quadrant_base <- ggplot(quad_summary, aes(x = !!sym(target1),
                                                 y = !!sym(target2))) +
      geom_point(aes(color = !!sym(target3))) +
      scale_color_gradient(low = gmb_palette[["lightblue"]], high = gmb_palette[["orange"]]) +
      ggtitle(paste("Scatterplot between", target1, "and", target2, "colored by", target3))
  }

  scatter_quadrant <- scatter_quadrant_base +  facet_wrap(grouping_var, ncol = facet_ncol) +
   ggpp::geom_quadrant_lines(xintercept = xintercept_coord,
                        yintercept = yintercept_coord) +
    labs(caption = caption) +
    theme(
      panel.background = element_blank(),
      strip.background = element_rect(linewidth = 0.5, color = "black", fill = "white"),
      strip.text.x = element_text(size = 12,face = "bold", color = "black"),
      plot.title = element_text(size = 16,face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12, color = "black"),
      plot.caption = element_text(size = 12)) +
    geom_text(data = quad_perc, aes(x = as.numeric(dir_axis[3]),
                                    y = as.numeric(dir_axis[4]), label = paste0(!!sym(QOI),"%"))) +
    stat_cor(method = cor_method, label.x = 0.2,label.y = 0.2)

  return(scatter_quadrant)
}
