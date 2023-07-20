#' Get composite scores when compared between target when compared with other gene/pathway
#' @export
#'
#' @import tidyverse
#' @import dplyr
#'
#'
#' @param df dataframe with expression values for each indication
#' @param target1 Official gene symbol of gene 1
#' @param target2_list Official list of gene symbols
#' @param QOI Which of the four quadrants are of interest for the question; Q1 (Top-left; high-low); Q2 (Top-right; high-high); Q3 (bottom-left; low-high) and Q4 (bottom right; low-low)
#' @param threshold_method Threshold method to use; median or mean or Q3 (75th quantile)
#' @param grouping_var groups object based on this variable; intr if not provided
#' @param filter_threshold threshold for selecting rows relevant in terms of composite scores
#'
#' @return dataframe with composite scores
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
#'# Define target1 and the list of genes to compare
#' target1 <- "Immune target"
#' target2_list <- c("Gene 1","Gene 2","Gene 3","Gene 4","Gene 5","Gene 6","Gene 7","Gene 8","Gene 9","Gene 10")
#'
#'# example1: obtain composite scores for Q2 quadrant for each comparison for default grouping var "intr"
#' df =get_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     threshold_method = "median")
#'# example2: obtain composite scores for Q2 quadrant for each comparison for different grouping variable
#' df = get_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     grouping_var="indication",
#'                     threshold_method = "median")
#'
#'# example3: obtain composite scores filtered based on threshold of composite score 30 in any of the the grouping variables
#' df = get_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     threshold_method = "median",
#'                     filter_threshold = 30)
#'
get_composite_scores <- function(df, target1, target2_list,
                                 QOI = c("Q1", "Q2", "Q3", "Q4","all"),
                                 threshold_method = c("median", "mean", "Q3"),
                                 grouping_var=NA,
                                 filter_threshold = NA){

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")
  ifelse(!is.na(filter_threshold), filter_threshold <- filter_threshold, filter_threshold <- 0)

  myoutlist = list()

  for (i in target2_list) {
    df_perc_pop <- get_perc_pop(df, target1=target1,target2 = i,
                                threshold_method = threshold_method,
                                QOI=QOI)
    myoutlist[[i]] <- df_perc_pop

  }
  quartile_data <- as.data.frame(myoutlist)


  if (QOI == "all") {

    quartile_data <-  quartile_data %>%
      select(contains("Q1"), contains("Q2"),contains("Q3"), contains("Q4"))

    quartile_filter <- quartile_data  %>%
      filter(rowSums(. > filter_threshold) > 0)

  }

  else {

    quartile_data <- select(quartile_data,contains(QOI))
    colnames(quartile_data) <- gsub("\\.Q.*","",colnames(quartile_data))
    quartile_filter <- quartile_data  %>%
      filter(rowSums(. > filter_threshold) > 0)

  }

  return(quartile_filter)
}
