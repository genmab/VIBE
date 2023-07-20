#' Returns a dataframe with percentage of samples per quadrant or Quadrant of interest (QOI)
#' Q1: Low target1; High target 2
#' Q2: High target1; High target 2
#' Q3: Low target1; Low target 2
#' Q4: High target1; Low target 2
#'
#'             #
#'     Q1      #       Q2
#'T            #
#'a            #
#'r            #
#'g ########################
#'e            #
#'t            #
#'2            #
#'     Q3      #       Q4
#'             #
#'         target 1
#'
#'
#'
#' @export
#'
#'
#' @param df the dataframe
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#' @param threshold_method Threshold method to use; median or mean or Q3 (75th quantile)
#' @param QOI Quadrant of interest among 4 quadrants
#' @param grouping_var groups object based on this variable; intr if not provided
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

#' # example1: percentage of samples with high expression of both targets (Q2)
#'get_perc_pop(df = df_harmonized,
#'             target1 = "Immune target",
#'             target2 = "Tumor target",
#'             threshold_method = "median", QOI = "Q2")
#'
#' # example2: percentage of samples for all quadrants of interest (QOI) per indication group
#'get_perc_pop(df = df_harmonized,
#'             target1 = "Immune target",
#'             target2 = "Tumor target",
#'             threshold_method = "median", QOI = "all",
#'             grouping_var = "indication")

get_perc_pop <- function(df, target1, target2,
                         threshold_method=c("median","mean","Q3"),
                         QOI=c("Q1","Q2","Q3","Q4","all"),
                         grouping_var=NA)
  {

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")

  quad_summary <- add_quadrant_info(df, target1, target2, threshold_method, grouping_var)
  percentage_per_quadrant <- calculate_percentage_per_quadrant(quad_summary, grouping_var)


  if (QOI=="all"){
    percentage_per_quadrant <- percentage_per_quadrant
  }
  else{
    percentage_per_quadrant <- percentage_per_quadrant[, c(QOI,grouping_var)]
  }
  return(percentage_per_quadrant)
}
