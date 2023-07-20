#' Plot a boxplot comparing samples from multiple groups (unmatched) for selected gene with Percentage of samples above median, Kruskal Wallis test, Fold change and sample counts
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggsignif
#' @import ggpubr
#' @importFrom stats setNames
#'
#'
#' @param df dataframe harmonized by VIBE's harmonize_df function, with the indication column being a concatenation of the study and primary_site columns of the XENA dataset.
#' @param gene name of genes for creating the plot
#' @param threshold_method the summary stat method to use, either "mean" or "median"
#' @param dataset character string with name of dataset, used in caption
#' @param grouping_var column name which you want to group data by (intr if left to NA)
#' @param plot_groups Character vector with values you wish to include in the plot, threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#' @param splitting_var column name of column that should be used to color plot boxplots with.
#' @param reverse TRUE or FALSE, if FALSE fold chanae is calculated as b/a, if TRUE as a/b

#'
#' @return boxplot comparing the statistical differences of the sample distribution per gene in GTEX vs TCGA
#'
#' @examples
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

#' # example1: calculate fold change between two groups and show as boxplots
#'plot_expression_box_split(df = df_harmonized,
#'                          gene = "Tumor target",
#'                          threshold_method = "median",
#'                          dataset = "VIBE dummy data (not real data!)",
#'                          splitting_var = "database")
#'
#' # example2: group x-axis by indication instead of intr
#'plot_expression_box_split(df = df_harmonized,
#'                          gene = "Tumor target",
#'                          threshold_method = "median",
#'                          dataset = "VIBE dummy data (not real data!)",
#'                          splitting_var = "database",
#'                          grouping_var = "indication")
#'
#' # example3: limit number of groups plotted on x-axis and reverse numerator/denominator
#' # as the output is ggplot2 based, you can use additional ggplot2 functionality
#'plot_expression_box_split(df = df_harmonized,
#'                          gene = "Tumor target",
#'                          threshold_method = "median",
#'                          dataset = "VIBE dummy data (not real data!)",
#'                          splitting_var = "database",
#'                          grouping_var = "indication",
#'                          plot_groups = c("BLCA", "CRC", "TNBC"),
#'                          reverse = TRUE) +
#'  ggplot2::expand_limits(y = 0)
#'

plot_expression_box_split <- function(df, gene, threshold_method, dataset = "Dataset", grouping_var = "intr", plot_groups = NA, splitting_var, reverse = FALSE) {

  # check that gene is in df
  if(!gene %in% df$symbol) {
    stop("gene ", gene, " is not in the symbol column. Dit you spell it correctly?")
  }

  # check that there are only two groups in the splitting_var column
  unique_splitting_vars = sort(unique(df[[splitting_var]]))

  if(length(unique_splitting_vars) != 2) {
    stop("This function calculates the fold change between two groups. You have ", length(unique_splitting_vars), " groups: ",
         paste(unique_splitting_vars, collapse = ", "), ". Pleas ensure your splitting_var column only contains two unique values.")
  }


  # create new split_groups var to calculate statistics
  if(any(str_detect(unique(df[[grouping_var]]), pattern = "\\[sep\\]"))) {
    stop("Your ", grouping_var,  "column contains the pattern [sep]. As this pattern is used to merge columns within the function, these ",
         grouping_var,
         "column values stop the function from working. Please remove the [sep] pattern from these values:",
         paste(grep("\\[sep\\]", x = unique(df[[grouping_var]]), value = TRUE),  collapse = ", ")
    )
  } else if (any(str_detect(unique(df[[splitting_var]]), pattern = "\\[sep\\]"))) {
    stop("Your ", splitting_var,  " column contains the pattern [sep]. As this pattern is used to merge columns within the function, these ",
         splitting_var,
         " column values stop the function from working. Please remove the [sep] pattern from these values:",
         paste(grep("\\[sep\\]", x = unique(df[[splitting_var]]), value = TRUE),  collapse = ", ")
    )
  }

  df_plot <- mutate(df, split_groups = paste(!!sym(grouping_var), !!sym(splitting_var), sep = "[sep]")) %>%
    filter(symbol == gene)


  # extract sample numbers
  df_numbers <- get_numbers(df_plot, grouping_var = "split_groups") %>%
    separate(split_groups, into = c(grouping_var, splitting_var), sep = "\\[sep\\]") %>%
    pivot_wider(names_from = all_of(splitting_var), values_from = nr) %>%
    mutate(across(.cols = all_of(unique_splitting_vars), .fns = ~ifelse(is.na(.x), 0, .x)),
           value = min(df_plot$value) - 0.1*(max(df_plot$value - min(df_plot$value))))

  df_numbers$nr <- apply(df_numbers[, unique_splitting_vars, drop = FALSE], 1, function(row) {
    paste(paste(unique_splitting_vars, row, sep = " = "), collapse = "\n")
  })


    # get fold change and kruskal-wallis p-value between two splittin_vars
  df_kw_fc <- calc_fc(df = df_plot, threshold_method = threshold_method, grouping_var = grouping_var, return_groups = plot_groups, splitting_var = splitting_var, reverse = reverse)


  # extract % over threshold
  label_pos = ceiling(max(df_plot$value))

  df_perc <- get_perc_gene(df = df_plot, gene = gene, method = threshold_method, grouping_var = "split_groups") %>%
    separate(split_groups, into = c(grouping_var, splitting_var), sep = "\\[sep\\]") %>%
    pivot_wider(names_from = all_of(splitting_var), values_from = 3) %>%
  mutate(across(.cols = all_of(unique_splitting_vars), .fns = ~ifelse(is.na(.x), 0, .x)),
         value = label_pos)

  df_perc$perc <- apply(df_perc[, unique_splitting_vars, drop = FALSE], 1, function(row) {
    paste(paste(unique_splitting_vars, paste0(row, "%"), sep = " = "), collapse = "\n")

  })

  # Combine percentage and statistics together
  df_kw_fc$fc_pval <- paste("FC:",round(df_kw_fc$fold_change,2),"\nPval:", df_kw_fc$p_value)
  df_perc_fc <- left_join(df_perc, df_kw_fc[,c(grouping_var,"fc_pval")], by = grouping_var)
  df_perc_fc$perc_stats <- paste(df_perc_fc$perc,"\n",df_perc_fc$fc_pval)


  # calculate threshold
  threshold <- get_threshold(df, gene = gene, method = threshold_method)[[1]]

  # other variables
  dataset_colors = stats::setNames(c(palette_light[["grey"]], palette_light[["green"]]), nm = unique_splitting_vars)
  unit2 = df$unit[1]

  # Filter data to plot only selected primary_site/tissues, while calculations are done on entire datset
  if(!any(is.na(plot_groups))) {

    if(!all(plot_groups %in% unique(df[[grouping_var]]))) {
      warning("The following variables assigned to plot_groups are not present in column ", grouping_var,":",
        paste(plot_groups[!plot_groups %in% unique(df[[grouping_var]])], collapse = ", "),
        ". Only plot_groups variables that are present in ", grouping_var, " will be plotted."
        )
    }

    df_plot <- dplyr::filter(df_plot, !!sym(grouping_var) %in% plot_groups)
    df_numbers <- dplyr::filter(df_numbers, !!sym(grouping_var) %in% plot_groups) %>%
      mutate(value = min(df_plot$value) - 0.1*(max(df_plot$value - min(df_plot$value))))
    df_perc_fc <- dplyr::filter(df_perc_fc, !!sym(grouping_var) %in% plot_groups)
  }

  # plot
  plot = ggplot(df_plot, aes(x = !!sym(grouping_var), y = value)) +
    geom_boxplot(mapping = aes(fill = !!sym(splitting_var))) +
    geom_text(data = df_numbers, aes(label = nr), size = 2) +
    geom_label(data = df_perc_fc, aes(label = perc_stats), size = 2) +
    geom_hline(yintercept = threshold,
               col = gmb_palette[["orange"]], linetype = "dashed") +
    expand_limits(y = label_pos+2) +
    labs(title = paste(gene, "expression in", unique_splitting_vars[1], "and", unique_splitting_vars[2]),
         caption = paste0("Dataset: ", dataset, "\n",
                          "The number of patient samples for each indication is depicted below the bars\n",
                          "Line: ", threshold_method, " gene expression over all samples is indicated with dashed line\n",
                          "FC: ", "Fold change are calculated per tissue using the formula ",
                          ifelse(reverse, paste0(unique_splitting_vars[1], "/", unique_splitting_vars[2]),
                                 paste0(unique_splitting_vars[2], "/", unique_splitting_vars[1])), "\n",
                          "Pval: ", " P vales are calculated based on the Kruskal Wallis test" )) +
    ylab(paste0("Expression log2(", unit2, ")")) +
    scale_fill_manual(values = dataset_colors) +
    ggpubr::theme_classic2() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


  return(plot)

}


