#' Calculate average expression of gene list combined for pathways
#' @export
#'
#' @import dplyr
#'
#'
#' @param df dataframe with expression values for each indication (as returned by harmonize_df function)
#' @param gene_panel_list Official gene symbol of genes compiled in list
#' @param grouping_cols groups object based on these columns; if default columns not provided
#'
#' @return dataframe with expression values for each gene/pathway and indication in long format
#'
#' @examples
#' # Creating lists ----------------------------------------------------------------------
#' pathway_1 <- c("Gene 14","Gene 30","Gene 20","Gene 22","Gene 29",
#' "Gene 23","Gene 1","Gene 32","Gene 33")
#' pathway_2 <- c("Gene 21","Gene 11","Immune target","Gene 18","Gene 8","Gene 16")
#' pathway_3 <- c("Gene 2" ,"Gene 3" ,"Gene 5" ,"Gene 6","Gene 7","Gene 10","Gene 15",
#' "Gene 17","Gene 18","Gene 24","Gene 25","Gene 28","Gene 31","Gene 35")
#' pathway_4 <- c('Gene 19',"Gene 12","Gene 26")
#' immune_target <- "Immune target"
#' tumor_target <- "Tumor target"
#' gene_panel_list <- list( pathway_1,pathway_2,pathway_3,pathway_4, immune_target, tumor_target)
#' names(gene_panel_list) <- c("pathway_1", "pathway_2","pathway_3","pathway_4",
#' "immune_target","tumor_target")
#'
#' # Calling function ----------------------------------------------------------------------
#'
#' df_harmonized <- harmonize_df(df = VIBE_data,
#' col_patientid = "patient_no",
#' col_sampleid = "analysis",
#' col_indication = "tumor",
#' col_treatment = "treatment_flag",
#' col_symbol = "gene",
#' col_value = "log2_tpm",
#' unit = "tpm")
#' df_avg_pathways <- harmonize_df_pathway(df_harmonized, gene_panel_list)


harmonize_df_pathway <- function(df, gene_panel_list, grouping_cols=NA){

  ifelse(!is.na(grouping_cols), grouping_cols <- grouping_cols,
         grouping_cols <- c("sampleid", "indication", "treatment", "intr"))

  myoutlist=list()

  for (i in 1:length(gene_panel_list)) {
    #print(primes_list[[i]])
    name=names(gene_panel_list)[[i]]

    df_sel <- df %>%
      filter(symbol %in% gene_panel_list[[i]]) %>%
      group_by_at(grouping_cols) %>%
      summarize(value = mean(value, na.rm=TRUE)) %>%
      as.data.frame()

    myoutlist[[i]] <- df_sel
  }

  names(myoutlist) <- names(gene_panel_list)

  avg_pathway = dplyr::bind_rows(myoutlist, .id = "symbol") %>% ungroup()
  avg_pathway$unit <- unique(df$unit)

  return(avg_pathway)
}
