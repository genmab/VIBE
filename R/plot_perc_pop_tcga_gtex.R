#' Plot percentage of samples per quadrant or Quadrant of interest (QOI)
#' @export
#'
#' @param df the dataframe with TCGA and GTEX data.
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#' @param threshold_method Threshold method to use; median or mean or Q3 (75th quantile)
#' @param QOI Quadrant of interest among 4 quadrants (Q1, Q2, Q3 OR Q4)
#' @param order True or False to sort the a axis in decreasing order ofr the quadrant of interest (Will not work for all quadrants)
#' @param dataset Name of your dataset for captioning
#' @param grouping_var column name which you want to group data by (intr if left to NA)
#' @param plot_groups Character vector with values you wish to include in the plot, threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#' @param fill_by column name of column that should be used to color plot bars with. Standard is set to "study"
#'
#' @return the dataframe with percentage of samples per quadrant or Quadrant of interest (QOI)


plot_perc_pop_tcga_gtex <- function (df, target1, target2, threshold_method = c("median", "mean", "Q3"), QOI = c("Q1", "Q2", "Q3", "Q4"), order = c("T", "F"),
                                     dataset =NA, grouping_var = NA, plot_groups = NA, fill_by = "study") {

  .Deprecated("plot_perc_pop_split")

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
    perc_pop_long <- filter(perc_pop_long, !!sym(grouping_var) %in% plot_groups)
    df <- filter(df, !!sym(grouping_var) %in% plot_groups)
  }

  # colors
  if("TCGA" %in% toupper(df[[fill_by]]) & "GTEX" %in% toupper(df[[fill_by]])) {
    fill_palette = c("TCGA" = palette_light[["green"]], "GTEX" = palette_light[["grey"]])
    df[[fill_by]] <- toupper(df[[fill_by]])
  } else {
    fill_palette = as.vector(gmb_palette)
  }

  # plot data
    direction_genes <- direction_and_axis_quadrants(df,
                                                    target1, target2, QOI = QOI)
    title = paste("Percentage of patients with", direction_genes[1],
                  target1, "and", direction_genes[2], target2, "expression")
    caption = paste("Cut-off: ", threshold_method, " gene expression over all tumor types in dataset \n",
                    "Dataset: ", dataset)
    quad_per_melt_select <- perc_pop_long[perc_pop_long$variable %in%
                                            QOI, ] %>%
      left_join(distinct(select(df, any_of(c(grouping_var, fill_by)))),
                by = grouping_var
      )

    if (order == "T") {
      quad_perc_bar <- ggplot(data = quad_per_melt_select,
                              aes(x = reorder(!!sym(grouping_var), -value), y = value, fill = !!sym(fill_by))) +
        geom_bar(stat = "identity") + ggtitle("Percentage of samples in each quadrant per indication and treatment") +
        geom_text(aes(label = value), vjust = -0.1,
                  position = position_dodge(0.9)) +
        scale_fill_manual(values = fill_palette) +
        theme(panel.background = element_blank(),
              strip.background = element_rect(linewidth = 0.5,
                                              color = "black", fill = "white"), strip.text.x = element_text(size = 8,
                                                                                                            face = "bold", color = "black"), plot.title = element_text(size = 12,
                                                                                                                                                                       face = "bold", hjust = 0.5), axis.text.y = element_text(size = 8,
                                                                                                                                                                                                                               color = "black"), axis.text.x = element_text(size = 8,
                                                                                                                                                                                                                                                                            color = "black", angle = 60, hjust = 1)) +
        ylab(paste0("Percentage of samples in ", QOI)) +
        xlab("Indication and Treatment") + labs(title = title,
                                                caption = caption)
    }
    else {
      quad_perc_bar <- ggplot(data = quad_per_melt_select,
                              aes(x = !!sym(grouping_var), y = value, fill = !!sym(fill_by))) +
        geom_bar(stat = "identity") + ggtitle("Percentage of samples in each quadrant per indication and treatment") +
        geom_text(aes(label = value), vjust = -0.1,
                  position = position_dodge(0.9)) +
        scale_fill_manual(values = fill_palette) +
        theme(panel.background = element_blank(),
              strip.background = element_rect(linewidth = 0.5,
                                              color = "black", fill = "white"), strip.text.x = element_text(size = 8,
                                                                                                            face = "bold", color = "black"), plot.title = element_text(size = 12,
                                                                                                                                                                       face = "bold", hjust = 0.5), axis.text.y = element_text(size = 8,
                                                                                                                                                                                                                               color = "black"), axis.text.x = element_text(size = 8,
                                                                                                                                                                                                                                                                            color = "black", angle = 60, hjust = 1)) +
        ylab(paste0("Percentage of samples in ", QOI)) +
        xlab("Indication and Treatment") + labs(title = title,
                                                caption = caption)
    }

  return(quad_perc_bar)
}
