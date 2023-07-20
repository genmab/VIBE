#' Function to calculate the number of samples per group in the data set
#' @export
#'
#' @param df object with expression data frame for each indication
#' @param grouping_var groups object based on this variable; intr if not provided
#' @return The number of samples per indication
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

#' # example1: sample numbers per intr group
#'get_numbers(df = df_harmonized)
#'
#' # example2: sample numbers per indiation group
#'get_numbers(df = df_harmonized, grouping_var = "indication")

get_numbers=function(df, grouping_var=NA){
  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")

  numbers_df <- df %>%
    dplyr::filter(symbol==df$symbol[1]) %>%
    group_by(across(c(!!sym(grouping_var),symbol))) %>%
    summarise(nr = n()) %>%
    as.data.frame() %>%
    dplyr::select(!!sym(grouping_var),nr)
  return(numbers_df)
}
