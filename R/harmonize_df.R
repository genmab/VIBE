#' harmonize dataframe for VIBE
#' @export
#'
#' @import dplyr
#' @import stringr
#'
#' @param df dataframe with expression values for each indication in long format
#' @param col_patientid column name for column that contains unique patient identifiers (optional column)
#' @param col_sampleid column name for colun that contains unique sample identifiers  (optional column)
#' @param col_indication column name for column that contains indications
#' @param col_treatment column name for column that contains indications (pre/post-treatment)  (optional column)
#' @param col_therapy column name for column that contains therapy (what type of treatment, e.g. chemo)  (optional column)
#' @param col_symbol column name for column that contains gene symbols
#' @param col_value column name for column that contains normalized count values
#' @param unit unit that is used for the normalized counts (e.g. tpm or cpm)
#' @param keep_extra_cols vector with column names if there are additional columns you wish to keep
#'
#' @return df dataframe in format that is appropriate for VIBE functions
#' @examples
#' df_harmonized <- harmonize_df(df = VIBE_data,
#' col_patientid = "patient_no",
#' col_sampleid = "analysis",
#' col_indication = "tumor",
#' col_treatment = "treatment_flag",
#' col_symbol = "gene",
#' col_value = "log2_tpm",
#' unit = "tpm")

harmonize_df <- function(df,
                         col_patientid = NA,
                         col_sampleid = NA,
                         col_indication,
                         col_treatment = NA,
                         col_therapy = NA,
                         col_symbol,
                         col_value,
                         unit,
                         keep_extra_cols = NA) {

  # keep following columns, these should always be there
  keep_cols = c("indication", "symbol", "value", "treatment")

  # keep&rename following optional columns if they are present
  if(!is.na(col_patientid)) {keep_cols = c(keep_cols, "patientid"); df <- rename(df, patientid = eval(col_patientid))}
  if(!is.na(col_sampleid)) {keep_cols = c(keep_cols, "sampleid"); df <- rename(df, sampleid = eval(col_sampleid))}
  if(!is.na(col_therapy)) {keep_cols = c(keep_cols, "therapy"); df <- rename(df, therapy = eval(col_therapy))}
  if(!is.na(col_treatment)) {df <- rename(df, treatment = eval(col_treatment))} else {  df$treatment <- NA  }

  # if there are additional columns to be kept, these will be added to the list here
  if(!is.na(keep_extra_cols[1])) {keep_cols = c(keep_cols, keep_extra_cols)}

  # restructure df
  df_harmonized <-  df %>%
    ungroup() %>%
    dplyr::rename(
      indication = eval(col_indication),
      symbol = eval(col_symbol),
      value = eval(col_value)
    )  %>%
    dplyr::select(all_of(keep_cols)) %>%
    dplyr::arrange(indication, desc(treatment)) %>%
    dplyr::mutate(treatment = tolower(treatment),
                  intr = paste(indication, str_extract(treatment, "pre|post")),
                  intr = ifelse(is.na(treatment), indication, intr),
                  intr = factor(intr, unique(intr)),
                  unit = unit
    )

  message("Returning dataframe with the following columns:", paste0(colnames(df_harmonized), collapse = ","))

  return(df_harmonized)
}
