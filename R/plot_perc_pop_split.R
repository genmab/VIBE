#' Plot percentage of samples per quadrant or Quadrant of interest (QOI)
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import stringr
#' @import tidyr
#' @importFrom stats setNames
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
#' @param splitting_var column name of column that will be used to create subgroups, %of samples above threshold and bars will be done for each unique combination of grouping_var and splitting_var
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

#' # example1: without splitting up x-axis groups
#'plot_perc_pop_split(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'                    threshold_method = "median", QOI = "Q1", order = "F", 
#'                    dataset = "VIBE dummy data (not real data!)")
#' # example2: split based on database & group x-axis by indication
#'plot_perc_pop_split(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'                    threshold_method = "median", QOI = "Q1", order = "T", 
#'                    dataset = "VIBE dummy data (not real data!)",
#'                    grouping_var = "indication",
#'                    splitting_var = "database")
#'
#' # example3: up to 9 unique splitting_var groups can be plotted
#'
#'groups_9 <- paste("group", c(1:9))
#'set.seed(42)
#'df_extra_groups <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(df_harmonized, patientid),
#'                                 groups = sample(groups_9, size = 1)))
#'
#'plot_perc_pop_split(df = df_extra_groups, target1 = "Immune target", target2 = "Tumor target",
#'                    threshold_method = "median", QOI = "Q1", order = "T", 
#'                    dataset = "VIBE dummy data (not real data!)",
#'                    grouping_var = "indication",
#'                    plot_groups = c("PDAC", "CERV"),
#'                    splitting_var = "groups")
#'

plot_perc_pop_split <- function (df, target1, target2, threshold_method = c("median", "mean", "Q3"), QOI = c("Q1", "Q2", "Q3", "Q4"), order = c("T", "F"),
                                 dataset = "Dataset", grouping_var = NA, plot_groups = NA, splitting_var = NA) {

  # set grouping variable to intr if it is NA
  ifelse(!is.na(grouping_var), grouping_var <- grouping_var,
         grouping_var <- "intr")

  # create new split_groups var to calculate statistics (unless splitting_var == NA)
  if(is.na(splitting_var)) {

    df <- mutate(df, split_groups = !!sym(grouping_var))

  } else{
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
    df <- mutate(df, split_groups = paste(!!sym(grouping_var), !!sym(splitting_var), sep = "[sep]"))
  }


  # calculate statistics
  perc_pop <- get_perc_pop(df, target1, target2, threshold_method,
                           QOI, grouping_var = "split_groups")


  if(is.na(splitting_var)) {
    perc_pop_long <- pivot_longer(perc_pop, cols = contains("Q"),
                                  names_to = "variable")
    perc_pop_long[[grouping_var]] <- perc_pop_long$split_groups
  } else {
    perc_pop_long <- pivot_longer(perc_pop, cols = contains("Q"),
                                  names_to = "variable") %>%
      separate(split_groups, into = c(grouping_var, splitting_var), sep = "\\[sep\\]", remove = FALSE) %>%
      mutate(split_groups = str_replace(split_groups, pattern = "\\[sep\\]", " "))

  }

  # filter data so only groups of interest are plotted
  if(!any(is.na(plot_groups))) {

    perc_pop_long <- filter(perc_pop_long, !!sym(grouping_var) %in% plot_groups)
    df <- filter(df, !!sym(grouping_var) %in% plot_groups)
  }

  # colors
  if(all(is.na(splitting_var))) {
    fill_palette = palette_light[["green"]]
  } else {
    split_groups <- sort(unique(df[[splitting_var]]))

    if(length(split_groups)>9) {
      stop("There are more than 9 unique values in splitting_var. plot_perc_pop_split cannot handle more than 9 unique values.")
    } else {
      fill_palette = c(palette_light[["grey"]], palette_light[["green"]], as.vector(gmb_palette[names(gmb_palette) != "lightgreen"]))[1:length(split_groups)] %>%
        stats::setNames(nm = split_groups)
    }

  }

  # plot data
  direction_genes <- direction_and_axis_quadrants(df,
                                                  target1, target2, QOI = QOI)
  title = paste("Percentage of patients with", direction_genes[1],
                target1, "and", direction_genes[2], target2, "expression")
  caption = paste("Cut-off: ", threshold_method, " gene expression over all tumor types in dataset \n",
                  "Dataset: ", dataset)

  if (order == "T") {
    quad_perc_bar <- ggplot(data = perc_pop_long,
                            aes(x = reorder(split_groups, -value), y = value)) +
      # geom_bar(stat = "identity") +
      geom_text(aes(label = value), vjust = -0.1,
                position = position_dodge(0.9)) +
      scale_fill_manual(values = fill_palette) +
      theme(panel.background = element_blank(),
            strip.background = element_rect(linewidth = 0.5,
                                            color = "black", fill = "white"),
            strip.text.x = element_text(size = 8,face = "bold", color = "black"),
            plot.title = element_text(size = 12,face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 8,color = "black"),
            axis.text.x = element_text(size = 8,color = "black", angle = 60, hjust = 1)) +
      ylab(paste0("Percentage of samples in ", QOI)) +
      xlab("Indication and Treatment") + labs(title = title,
                                              caption = caption)
  }  else {
    quad_perc_bar <- ggplot(data = perc_pop_long,
                            aes(x = split_groups, y = value)) +
      # geom_bar(stat = "identity") +
      geom_text(aes(label = value), vjust = -0.1,
                position = position_dodge(0.9)) +
      scale_fill_manual(values = fill_palette) +
      theme(panel.background = element_blank(),
            strip.background = element_rect(linewidth = 0.5,
                                            color = "black", fill = "white"),
            strip.text.x = element_text(size = 8,face = "bold", color = "black"),
            plot.title = element_text(size = 12,face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 8,color = "black"),
            axis.text.x = element_text(size = 8,color = "black", angle = 60, hjust = 1)) +
      ylab(paste0("Percentage of samples in ", QOI)) +
      xlab("Indication and Treatment") + labs(title = title,
                                              caption = caption)
  }

  # fill color depends on whether splitting_var is used
  if(is.na(splitting_var)) {
    p_quad_perc_bar <- quad_perc_bar + geom_bar(stat = "identity", fill = fill_palette)
  } else {
    p_quad_perc_bar <- quad_perc_bar + geom_bar(stat = "identity", mapping = aes(fill = !!sym(splitting_var)))
  }

  return(p_quad_perc_bar)

}
