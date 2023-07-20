#' Function to calculate the threshold expression values for a specified gene
#' By calculating the median, mean or Q3 expression over all samples in the data set
#' @export
#'
#' @param df object with expression data frame as harmonized by @harmonized_df
#' @param gene Official gene symbol
#' @param method "median", "mean" or "Q3"; The method used for calculating the threshold value
#'
#' @return The threshold expression value for the specified gene
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

#' # example1: median threshold
#'get_threshold(df = df_harmonized, gene = "Immune target", method = "median")
#'
#' # example2: mean threshold
#'get_threshold(df = df_harmonized, gene = "Immune target", method = "mean")

get_threshold = function(df, gene, method = c("median","mean","Q3"))
{
  if(method == "median")
  {
    threshold = df %>%
      ungroup() %>%
      dplyr::filter(symbol==gene) %>%
      summarise(round(median(value, na.rm = TRUE),5)) %>%
      as.data.frame()
  }

  if(method == "mean")
  {
    threshold = df %>%
      ungroup() %>%
      dplyr::filter(symbol==gene) %>%
      summarise(round(mean(value, na.rm = TRUE),5)) %>%
      as.data.frame()
  }

  if(method == "Q3")
  {
    threshold = df %>%
      ungroup() %>%
      dplyr::filter(symbol==gene) %>%
      summarise(round(quantile(value, na.rm = TRUE)[4],5)) %>%
      as.data.frame()
  }
  return(threshold)
}



