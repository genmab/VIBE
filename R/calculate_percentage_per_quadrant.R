#' Calculate percentage of sample in each quadrant
#' @export
#'
#' @param quad_summary the wide format of the data including the quadrant information column as calculated by @add_quadrant_info function.
#' @param grouping_var groups object based on this variable; intr if not provided
#'
#' @return tibble with combined data_wide and quadrant information (as an additional column)
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
#' # example 1
#' df_quadrant <- add_quadrant_info(df = df_harmonized,
#'                                  target1 = "Immune target",
#'                                  target2 = "Tumor target",
#'                                  threshold_method = "median")
#'
#' calculate_percentage_per_quadrant(df_quadrant)
#'
#' # example 2
#' df_quadrant <- add_quadrant_info(df = df_harmonized,
#'                                  target1 = "Immune target",
#'                                  target2 = "Tumor target",
#'                                  grouping_var = "indication",
#'                                  threshold_method = "median")
#'
#' calculate_percentage_per_quadrant(df_quadrant, grouping_var = "indication")

calculate_percentage_per_quadrant <- function(quad_summary, grouping_var=NA){

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")

  ## Calculating the percentage of the samples per quadrant
  quad_counts <- as.data.frame.matrix(table(quad_summary[[grouping_var]], quad_summary$quadrant))
  colnames(quad_counts) <- c("Q1","Q2","Q3","Q4")
  quad_perc <- round(quad_counts/rowSums(quad_counts)*100,0)
  quad_perc[[grouping_var]] <- rownames(quad_perc)

  return(quad_perc)
}
