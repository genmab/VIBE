#' Heatmap representing the median expression of indications
#' @export
#'
#' @import ComplexHeatmap
#' @import circlize
#'
#' @param df dataframe in long format containing intr (indication + treatment), value and symbol columns
#' @param method the summary stat method to use, either "mean" or "median", median is standard
#' @param grouping_var groups object based on this variable; intr if not provided
#' @param cluster_switch if you would like to cluster your rows based on pathways etc. default is NULL
#' @param cluster_columns if the heatmap grouping variables are to be clustered (default is TRUE)
#' @param scale scale the expression of genes in range 0-1 (default is TRUE)
#' @param column_title provide a text for title of your heatmap (default is provided)
#' @param plot_groups Character vector with values you wish to include in the plot, threshold will be calculated
#' for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#' @param order_4_heatmap  Character vector of the order of genes in heatmap (helpful for reassigning rowsplit)
#'
#' @return heatmap with cells colored based on the percentage of samples darkblue the lowest and magenta the highest
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
#'# example1: Heatmap of expression data for each comparison for default grouping var "intr"
#'heatmap_sample_expression(df_harmonized, method = "median")
#'
#'# example2: Heatmap of expression data for each comparison for grouping var "indication"
#'heatmap_sample_expression(df_harmonized, method = "mean",  grouping_var="indication")
#'
#'# example3: Splitting the heatmap in groups
#'
#'rowsplit <-rep(c("Group1", "Group2","Group3", "Group4","IT","TT"), c(9,8,9,7,1,1))
#'heatmap_sample_expression(df_harmonized, method="median", cluster_switch=rowsplit, plot_groups = c("NSCLC","RCC"))
#'
#'# example4: Reordering the heatmap and redefining the groups based on user-defined pathways
#'
#'reorder_list <- c("Gene 14","Gene 30","Gene 20","Gene 22","Gene 29","Gene 23","Gene 1","Gene 32","Gene 33","Gene 21","Gene 11","Gene 18","Gene 8","Gene 16","Gene 2" ,"Gene 3" ,"Gene 5" ,"Gene 6","Gene 7","Gene 10","Gene 15","Gene 17","Gene 24","Gene 25","Gene 28","Gene 31","Gene 35",'Gene 19',"Gene 12","Gene 26","Gene 4","Gene 9","Immune target","Tumor target")
#'
#'rowsplit_new <-rep(c("Pathway1", "Pathway2","Pathway3", "Pathway4","G4","G9","IT","TT"), c(9,5,13,3,1,1,1,1))
#'heatmap_sample_expression(df_harmonized, method="median", cluster_switch=rowsplit_new, order_4_heatmap = reorder_list)
#'
#'
#'
heatmap_sample_expression <- function(df, method="median",grouping_var=NA,
                                      cluster_switch=NULL,
                                      scale=TRUE,
                                      cluster_columns=TRUE,
                                      column_title=NA,
                                      plot_groups = NA,
                                      order_4_heatmap = NA)
{
  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")


  gene_exp = get_expr_table(df, method, grouping_var)

  df_exp <-  gene_exp[!(row.names(gene_exp) %in% c("median_threshold","mean_threshold")),] #removes last row

  expr_mat <- df_exp
  rownames(expr_mat) <- expr_mat[[grouping_var]]

  expr_mat1 <- expr_mat[,-which(names(expr_mat) %in% c(grouping_var,"nr"))]
  expr_mat1 <- mutate_all(expr_mat1, function(x) as.numeric(as.character(x)))
  if (scale == TRUE){
  ## Rescale each column to range between 0 and 1
  scaled_expr_mat <- apply(expr_mat1, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
  # transpose matrix for plotting
  expr_4_plot <- t(scaled_expr_mat)

  ifelse(!is.na(column_title), column_title <- column_title, column_title <- "Heatmap representing scaled expression of the genes in selected group")
  }
  else if (scale == FALSE){
    expr_4_plot <- expr_mat1

    ifelse(!is.na(column_title), column_title <- column_title, column_title <- "Heatmap representing expression of the genes in selected group")
  }

   if (!is.na(plot_groups[1])) {
   expr_4_plot <- expr_4_plot[, grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), colnames(expr_4_plot))]
     }

  ifelse(!is.na(order_4_heatmap),
         expr_4_plot <- expr_4_plot[match(order_4_heatmap, row.names(expr_4_plot)), ],
         expr_4_plot <- expr_4_plot[order(row.names(expr_4_plot)), ])

  heatmap <- Heatmap(expr_4_plot,
          cluster_rows = FALSE,
          cluster_columns = cluster_columns,
          row_split = cluster_switch,
          show_column_names = TRUE,
          column_title = column_title,
          col=colorRamp2(c(min(expr_4_plot),max(expr_4_plot)/2,max(expr_4_plot)),c("#e6f5f4","#99d8d4","#004f4a")),
          column_names_side = "bottom",
          heatmap_legend_param = list(title = "expression"))

  return(heatmap)
}
