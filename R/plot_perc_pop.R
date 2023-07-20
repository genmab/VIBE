#' Plot percentage of samples per quadrant or Quadrant of interest (QOI)
#' @export
#'
#' @import tidyverse
#'
#' @param df the dataframe
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#' @param threshold_method Threshold method to use; median or mean or Q3 (75th quantile)
#' @param QOI Quadrant of interest among 4 quadrants
#' @param order True or False to sort the a axis in decreasing order ofr the quadrant of interest (Will not work for all quadrants)
#' @param dataset Name of your dataset for captioning
#' @param grouping_var column name which you want to group data by (intr if left to NA)
#' @param plot_groups Character vector with values you wish to include in the plot, threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#'
#' @return the dataframe with percentage of samples per quadrant or Quadrant of interest (QOI)
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

#' # example1: percentage of samples in Q1 (low target1 and high target 2)
#'plot_perc_pop(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'              threshold_method = "median", QOI = "Q1",
#'              dataset = "VIBE dummy data", order = "F")
#' # example2: percentage of samples in Q1 (low target1 and high target 2) sorted
#'plot_perc_pop(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'              threshold_method = "median",
#'              QOI = "Q1", dataset = "VIBE dummy data", order = "T")
#' # example3: percentage of all quadrants per indication
#'plot_perc_pop(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'              threshold_method = "median",
#'              QOI = "all", dataset = "VIBE dummy data", order = "T",
#'              grouping_var = "indication", plot_groups = c("CRC", "PRAD", "UTEN"))


plot_perc_pop <- function (df, target1, target2,
                           threshold_method = c("median", "mean", "Q3"),
                           QOI = c("Q1", "Q2", "Q3", "Q4", "all"),
                           order = c("T", "F"), dataset=NA,
                           grouping_var = NA, plot_groups = NA) {

  # set grouping variable to intr if it is NA
  ifelse(!is.na(grouping_var), grouping_var <- grouping_var,
         grouping_var <- "intr")

  ifelse(!is.na(dataset), dataset <- dataset, dataset <- "Dataset")

  # calculate statistics
  perc_pop <- get_perc_pop(df, target1, target2, threshold_method,
                            QOI, grouping_var = grouping_var)
  perc_pop_long <- pivot_longer(perc_pop, cols = contains("Q"),
                                names_to = "variable")

  # filter data so only groups of interest are plotted
  if(!is.na(plot_groups[1])) {
    perc_pop_long <- perc_pop_long %>% dplyr::filter(grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), !!sym(grouping_var)))
    df <- df %>% dplyr::filter(grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), !!sym(grouping_var)))
  }

  # plot data
  if (QOI == "all") {
    title = paste("Percentage of patients in all quadrants of",
                  target1, "and", target2, "expression")
    caption = paste("Cut-off: ", threshold_method, " gene expression over all tumor types in dataset \n",
                    "Dataset: ", dataset)
    quad_perc_bar <- ggplot(data = perc_pop_long, aes(x = !!sym(grouping_var),
                                                      y = value, fill = variable)) + geom_bar(position = "dodge",
                                                                                              stat = "identity") + ggtitle("Percentage of samples in each quadrant per indication and treatment") +
      geom_text(aes(label = value), vjust = -0.1, position = position_dodge(0.9)) +
      scale_fill_manual(values = as.vector(gmb_palette)) + labs(fill = "Quadrants") +
      theme(legend.position = "bottom", legend.direction = "horizontal",
            panel.background = element_blank(),
            strip.background = element_rect(linewidth = 0.5,color = "black", fill = "white"),
            strip.text.x = element_text(size = 8, face = "bold", color = "black"),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 8, color = "black"),
            axis.text.x = element_text(size = 8, color = "black", angle = 60, hjust = 1)) +
      ylab("Percentage of samples per quadrant") + coord_flip() +
      labs(title = title, caption = caption)
  }
  else {
    direction_genes <- direction_and_axis_quadrants(df,
                                                    target1, target2, QOI = QOI)
    title = paste("Percentage of patients with", direction_genes[1],
                  target1, "and", direction_genes[2], target2, "expression")
    caption = paste("Cut-off: ", threshold_method, " gene expression over all tumor types in dataset \n",
                    "Dataset: ", dataset)
    quad_per_melt_select <- perc_pop_long[perc_pop_long$variable %in%
                                            QOI, ]
    if (order == "T") {
      quad_perc_bar <- ggplot(data = quad_per_melt_select,
                              aes(x = reorder(!!sym(grouping_var), -value), y = value, fill = variable)) +
        geom_bar(stat = "identity") + ggtitle("Percentage of samples in each quadrant per indication and treatment") +
        geom_text(aes(label = value), vjust = -0.1,
                  position = position_dodge(0.9)) + scale_fill_manual(values = gmb_palette[["lightgreen"]]) +
        theme(legend.position = "None", panel.background = element_blank(),
              strip.background = element_rect(linewidth = 0.5, color = "black", fill = "white"),
              strip.text.x = element_text(size = 8,face = "bold", color = "black"),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.text.y = element_text(size = 8,color = "black"),
              axis.text.x = element_text(size = 8, color = "black", angle = 60, hjust = 1)) +
        ylab(paste0("Percentage of samples in ", QOI)) +
        xlab("Indication and Treatment") + labs(title = title,
                                                caption = caption)
    }
    else {
      quad_perc_bar <- ggplot(data = quad_per_melt_select,
                              aes(x = !!sym(grouping_var), y = value, fill = variable)) +
        geom_bar(stat = "identity") + ggtitle("Percentage of samples in each quadrant per indication and treatment") +
        geom_text(aes(label = value), vjust = -0.1,
                  position = position_dodge(0.9)) + scale_fill_manual(values = gmb_palette[["lightgreen"]]) +
        theme(legend.position = "None", panel.background = element_blank(),
              strip.background = element_rect(linewidth = 0.5, color = "black", fill = "white"),
              strip.text.x = element_text(size = 8, face = "bold", color = "black"),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.text.y = element_text(size = 8,color = "black"),
              axis.text.x = element_text(size = 8, color = "black", angle = 60, hjust = 1)) +
        ylab(paste0("Percentage of samples in ", QOI)) +
        xlab("Indication and Treatment") + labs(title = title,
                                                caption = caption)
    }
  }
  return(quad_perc_bar)
}
