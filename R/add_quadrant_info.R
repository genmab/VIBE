#' Calculate the positioning of each sample in a quadrant based on the thresholds calculated by the median or 75th quantile
#' @export
#'
#' @import dplyr
#'
#'
#' @param df dataframe with expression values for each indication
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#' @param threshold_method Threshold method to use; "median" or "Q3" (75th quantile)
#' @param grouping_var groups object based on this variable; intr if not provided
#'
#' @return tibble with combined data_wide and quadrant information (as an additional column)
#' @examples
#' # prepare data for use with VIBE functions
#' data("VIBE_data")
#' df_harmonized = harmonize_df(df = VIBE_data,
#'                              col_patientid = "patient_no",
#'                           col_sampleid = "analysis",
#'                           col_indication = "tumor",
#'                           col_treatment = "treatment_flag",
#'                           col_symbol = "gene",
#'                           col_value = "log2_tpm",
#'                           unit = "tpm",
#'                           keep_extra_cols = c("database"))
#'
#' # example1: based on median threshold per intr group.
#' # Note: add_quadrant_info adds one column to the input dataframe called quadrant
#' df_quadrant_median <- add_quadrant_info(df = df_harmonized,
#'                                         target1 = "Immune target",
#'                                      target2 = "Tumor target",
#'                                      threshold_method = "median")
#' df_quadrant_median[c("patientid", "treatment", "Immune target", "Tumor target", "quadrant")]
#'
#'
#' # example2: 75th percentile threshold per intr group
#' df_quadrant_q3 <- add_quadrant_info(df = df_harmonized,
#'                                  target1 = "Immune target",
#'                                  target2 = "Tumor target",
#'                                  threshold_method = "Q3")
#' df_quadrant_q3[c("patientid", "treatment", "Immune target", "Tumor target", "quadrant")]
#'
#' # example3: median threshold calculated per indication group
#' df_quadrant <- add_quadrant_info(df = df_harmonized,
#'                               target1 = "Immune target",
#'                               target2 = "Tumor target",
#'                               threshold_method = "median",
#'                               grouping_var = "indication")
#' df_quadrant[c("patientid", "treatment", "Immune target", "Tumor target", "quadrant")]
#'
add_quadrant_info <- function (df, target1, target2, threshold_method, grouping_var = NA)
{

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var,
         grouping_var <- "intr")

  ifelse(!is.na(grouping_var),
         df_wide <- create_wide_format(df, additional_columns = grouping_var),
         df_wide <- create_wide_format(df))

  xintercept_coord <- get_threshold(df, target1, method = threshold_method)[1,
                                                                            1]
  yintercept_coord <- get_threshold(df, target2, method = threshold_method)[1,
                                                                            1]
  quad_summary <- df_wide %>% mutate(quadrant = sample_location_per_quadrant(x = df_wide[[target1]],
                                                                             y = df_wide[[target2]], xintercept = xintercept_coord,
                                                                             yintercept = yintercept_coord)) %>% group_by(!!sym(grouping_var),
                                                                                                                          quadrant) %>% mutate(across(contains("label"), max))
  return(quad_summary)
}

