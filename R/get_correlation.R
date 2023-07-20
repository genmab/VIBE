#' Calculate correlation between two genes
#' @export
#'
#' @import dplyr
#' @import tidyverse
#'
#'
#' @param df dataframe with expression values for each indication (as obtained from Profiler)
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#' @param cor_method Correlation method used
#' @param grouping_var groups object based on this variable; intr if not provided
#'
#' @return tibble with correlation value
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
#'# example 1: pearson correlation between target1 and target2 per intr group
#'get_correlation(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'cor_method = "pearson")
#'
#'# example 2: spearman correlation between target1 and target2 per indication group
#'get_correlation(df = df_harmonized, target1 = "Immune target", target2 = "Tumor target",
#'                cor_method = "spearman", grouping_var = "indication")
#'
get_correlation <- function(df, target1, target2, cor_method, grouping_var=NA){

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")

  df_wide<- create_wide_format(df)

  correlations_by_group  <- df_wide  %>%
    group_by(!!sym(grouping_var)) %>%
    summarize(correlation = cor(!!sym(target1), !!sym(target2), method = cor_method))

  return(correlations_by_group)
}
