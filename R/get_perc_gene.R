#' Function to calculate the percentage of samples above the threshold expression per group in the data set
#' @export
#'
#' @param df object with expression data frame for each indication
#' @param gene official gene symbol of the gene
#' @param method "median" or "mean"; The method used for calculating the threshold value.
#' @param grouping_var groups object based on this variable; intr if not provided

#'
#' @return The percentage of samples above the threshold
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

#' # example1: percentage of samples above median threshold per intr group
#'get_perc_gene(df = df_harmonized, gene = "Immune target",
#'              method = "median")
#'
#' # example2: percentage of samples above mean threshold per indication group
#'get_perc_gene(df = df_harmonized, gene = "Immune target",
#'              method = "mean", grouping_var = "indication")


get_perc_gene= function(df, gene, method=c("median","mean","Q3"),grouping_var=NA)

{
  threshold = get_threshold(df, gene,method)

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")

  col_name = paste0("pct_above_", method, "_",gene,"_expression")

  df_gene =      df %>%
    dplyr::filter(symbol==gene) %>%
    mutate(positive=(value >= threshold[1,])) %>%
    group_by_at(grouping_var) %>%
    summarise(!!col_name := round(100*mean(positive),2))
  return(df_gene)
}
