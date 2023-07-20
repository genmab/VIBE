#' Calculate the fold change and p-value for each symbol between two splitting_var groups per grouping_var group
#' @export
#'
#' @import dplyr
#'
#' @param df dataframe harmonized by VIBE's harmonize_df function
#' @param threshold_method the summary stat method to use, either "mean" or "median"
#' @param grouping_var column name which you want to group data by (intr if left to NA)
#' @param return_groups Character vector with grouping_var values you wish to include in the results table
#' @param splitting_var column name which splits the dataset into two groups between which the fold change is calculated
#' @param reverse TRUE or FALSE, if FALSE fold chanae is calculated as b/a, if TRUE as a/b
#'
#' @return table with fold change and Kruskal Wallis p-value per unique primary_site and gene combination
#'
#' @examples
#' # prepare data for use with VIBE functions
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
#'# example1: calculate foldchange between two databases, database2/database1
#'unique(df_harmonized$database)
#'calc_fc(df_harmonized, threshold_method = "median", splitting_var = "database")
#'
#'# example2: calculate foldchange between two databases, database1/database2
#'unique(df_harmonized$database)
#'calc_fc(df_harmonized, threshold_method = "median", splitting_var = "database",
#'        reverse = TRUE)
#'
#'# example3: group by indication instead of intr
#'unique(df_harmonized$database)
#'calc_fc(df_harmonized, threshold_method = "median", splitting_var = "database",
#'        grouping_var = "indication")


calc_fc <- function (df, threshold_method = c("mean", "median"), grouping_var = "intr", return_groups = NA, splitting_var, reverse = FALSE) {


  # check that there are only two groups in the splitting_var column
  unique_splitting_vars = sort(unique(df[[splitting_var]]))

  if(length(unique_splitting_vars) != 2) {
    stop("This function calculates the fold change between two groups. You have ", length(unique_splitting_vars), " groups: ",
         paste(unique_splitting_vars, collapse = ", "), ". Pleas ensure your splitting_var column only contains two unique values.")
  }

  # modify return_primary_site if it is set to NA
  if(all(is.na(return_groups))) {
    return_groups = as.vector(unique(df[[grouping_var]]))
    return_groups = return_groups[!is.na(return_groups)]
  }

  # filter out grouping_var values where there is data in only one out of two splitting_var groups
  group1_values = unique(filter(df, !!sym(splitting_var) == unique_splitting_vars[1])[[grouping_var]])
  group2_values = unique(filter(df, !!sym(splitting_var) == unique_splitting_vars[2])[[grouping_var]])

  list_drop <- c(setdiff(group1_values, group2_values),
                                setdiff(group2_values, group1_values))

  df_filtered <- dplyr::filter(df, !(!!sym(grouping_var) %in% list_drop), !!sym(grouping_var) %in% return_groups) %>%
    rename(splitting_var = !!sym(splitting_var),
           grouping_var = !!sym(grouping_var))

  # Kruskal wallis between two splitting_var groups
  kw <- df_filtered %>%
    group_by(grouping_var, symbol) %>%
    summarise(p_value = kruskal.test(value ~ splitting_var)$p.value) %>%
    ungroup() %>%
    mutate(p_value = ifelse(p_value<0.01, format(p_value, scientific = TRUE, digits = 3), round(p_value, digits = 2))
    )

  # Fold change between two splitting_var groups
  fc <- df_filtered %>%
    arrange(splitting_var) %>%
    group_by(grouping_var, symbol, splitting_var) %>%
    summarise(stat = if (threshold_method == "mean") {
      mean(value)
    }
    else if (threshold_method == "median") {
      median(value)
    }, y = max(value)) %>%
    group_by(grouping_var, symbol) %>%
    summarise(fold_change = ifelse(reverse, stat[1]/stat[2], stat[2]/stat[1])) %>% ungroup()

  return_na <- list_drop[list_drop %in% return_groups]

  df_result = left_join(kw, fc, by = c("grouping_var", "symbol")) %>%
    dplyr::bind_rows(data.frame(grouping_var = return_na)) %>%
    arrange(grouping_var)


  colnames(df_result)[colnames(df_result) == "grouping_var"] <- grouping_var


  return(df_result)

}

