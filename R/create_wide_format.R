#' Creates a wide format of the dataframe df
#' @export
#'
#' @import dplyr
#'
#' @param df dataframe with expression values for each indication (output of @harmonize_df )
#' @param additional_columns vector with column names to be kept as well (is used for grouping_var in some other vibe functions)
#'
#' @return dataframe with wide format of the data from profiler
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
#'# example 1
#'df_wide <- create_wide_format(df_harmonized)
#'dim(df_harmonized)
#'dim(df_wide)
#'dplyr::glimpse(df_wide)
#'
#'# example 2 keep additional column database
#'df_wide_extra <- create_wide_format(df_harmonized, additional_columns = "database")
#'dim(df_wide_extra)
#'dplyr::glimpse(df_wide_extra)

create_wide_format <- function (df, additional_columns="")
{
  df_wide <- df %>% ungroup() %>% dplyr::select(any_of(c("sampleid", "patientid", "treatment", "symbol", "value", "indication", "intr", "unit", additional_columns))) %>%
    pivot_wider(names_from = symbol, values_from = value)
  return(df_wide)
}


