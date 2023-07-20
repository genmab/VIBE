#' Calculate the fold change and p-value between GTEX and TCGA samples of the same primary site
#' @export
#'
#' @import dplyr
#'
#' @param df dataframe harmonized by VIBE's harmonize_df function, with the indication column being a concatenation of the study and primary_site columns of the XENA dataset.
#' @param threshold_method the summary stat method to use, either "mean" or "median"
#' @param return_primary_site Character vector with primary_site values (name of tissue column in Xena dataset) you wish to include in the results table
#'
#' @return table with fold change and Kruskal Wallis p-value per unique primary_site and gene combination
#'
#'


calc_fc_gtex_tcga <- function (df, threshold_method, return_primary_site = NA) {

  .Deprecated(new = "calc_fc")

  # modify return_primary_site if it is set to NA
  if(all(is.na(return_primary_site))) {
    return_primary_site = unique(df$primary_site)
    return_primary_site = return_primary_site[!is.na(return_primary_site)]
  }

  # Ensure that there is a study and tissue column, as this column is not required to be retained in the harmonize_df function
  df_table <- df %>% mutate(tissue = str_remove_all(intr, "GTEX|TCGA") %>% str_trim(),
                            study = str_extract(intr, "GTEX|TCGA"))

  # filter out primary_site values where there is no matched GTEX/TCGA data
  sites_gtex = unique(filter(df_table, study == "GTEX")$tissue)
  sites_tcga = unique(filter(df_table, study == "TCGA")$tissue)
  list_drop_primary_tissue <- c(setdiff(sites_gtex, sites_tcga),
                                setdiff(sites_tcga, sites_gtex))

  df_table_filtered <- dplyr::filter(df_table, !tissue %in% list_drop_primary_tissue, tissue %in% return_primary_site)

  # Kruskal wallis between GTEX primary site vs TCGA primary site
  kw <- df_table_filtered %>%
    group_by(tissue, symbol) %>%
    summarise(p_value = kruskal.test(value ~ study)$p.value) %>%
    ungroup() %>%
    mutate(p_value = ifelse(p_value<0.01, format(p_value, scientific = TRUE, digits = 3), round(p_value, digits = 2))
           )


  # Fold change between GTEX primary site vs TCGA primary site
  fc <- df_table_filtered %>%
    arrange(study) %>%
    group_by(tissue, symbol, study) %>%
    summarise(stat = if (threshold_method == "mean") {
      mean(value)
    }
    else if (threshold_method == "median") {
      median(value)
    }, y = max(value)) %>% group_by(tissue, symbol) %>%
    summarise(fold_change = stat[2]/stat[1]) %>% ungroup()

  return_na_primary_tissue <- list_drop_primary_tissue[list_drop_primary_tissue %in% return_primary_site]

  df_result = left_join(kw, fc, by = c("tissue", "symbol")) %>%
    dplyr::bind_rows(data.frame(tissue = return_na_primary_tissue)) %>%
    arrange( tissue)

  return(df_result)

}

