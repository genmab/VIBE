#' create table with sample numbers and median expression per intr group (indication + treatment)
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @param df dataframe in long format containing intr (indication + treatment), value and symbol columns
#' @param method the summary stat method to use, either "mean" or "median"
#' @param grouping_var groups object based on this variable; intr if not provided
#' @return sample number, median expression and median thresholds of all available genes per intr group
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

#' # example1: median expression of each intr group and all genes in df_harmonized
#'get_expr_table(df = df_harmonized, method = "median")
#'
#' # example2: mean expression of each indication group and all genes in df_harmonized
#'get_expr_table(df = df_harmonized, method = "mean", grouping_var = "indication")
#'


get_expr_table =function(df, method = c("mean", "median","Q3"), grouping_var = NA){

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")
  nr = get_numbers(df, grouping_var)

  stat = df %>%
    dplyr::select(!!sym(grouping_var), symbol, value) %>%
    group_by(!!sym(grouping_var), symbol) %>%
    {if(method == "mean")
    {summarise(., stat = round(mean(value, na.rm = TRUE),2))
    } else if (method == "median")
    {summarise(., stat = round(median(value, na.rm = TRUE),2))}
      else if (method == "Q3")
      {summarise(., stat = round(quantile(value, na.rm = TRUE)[4],2))}} %>%
    pivot_wider(names_from = symbol, values_from = stat) %>%
    arrange(!!sym(grouping_var))  %>%
    as.data.frame() %>%
    ungroup()


  # Generate threshold df to be added to the median expression table
  gene_panel <- unique(df$symbol)
  thresholds_genes <- data.frame()
  for(i in gene_panel){
    aa <- as.data.frame(get_threshold(df,i,method))
    aa$symbol=i
    thresholds_genes <- rbind(thresholds_genes,aa)
  }
  colnames(thresholds_genes) <- c(paste0(method,"_threshold"),"symbol")

  row_add <- data.frame(c(round(nrow(df)/length(gene_panel),0),"All samples"),c("nr",grouping_var))
  colnames(row_add) <- c(paste0(method,"_threshold"),"symbol")

  thresholds_genes <- rbind(row_add, thresholds_genes)

  thresholds_genes <- as.data.frame(t(thresholds_genes))
  colnames(thresholds_genes) <- thresholds_genes[2,]

  df_stat = stat %>% full_join(nr,by=grouping_var) %>% relocate(!!sym(grouping_var),nr)


  # Add threshold df to the median expression table
  df_expr_threshold <- rbind(df_stat, thresholds_genes[1,][names(df_stat)])

  return(df_expr_threshold)
}


