#' Heatmap representing the percentage of samples above median expression
#' @export
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#'
#' @param df dataframe in long format containing intr (indication + treatment), value and symbol columns
#' @param method the summary stat method to use, either "mean" or "median"
#' @param grouping_var groups object based on this variable; intr if not provided
#' @param cluster_switch if you would like to cluster your rows based on pathways etc. default is NULL
#' @param cluster_columns if the heatmap grouping variables are to be clustered (default is TRUE)
#' @param column_title provide a text for title of your heatmap (default is provided)
#' @param plot_groups Character vector which groups from column defined at grouping_var you wish to include in the plot. Median expression will be calculated from all groups in df.
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
#'heatmap_samples_above_median(df_harmonized, method = "median")
#'
#'# example2: Heatmap of expression data for each comparison for grouping var "indication"
#'heatmap_samples_above_median(df_harmonized, method = "median", grouping var "indication")
#'
#'# example3: Selecting the groups to be plotted
#'heatmap_samples_above_median(df_harmonized, method = "median", plot_groups = c("NSCLC","RCC))
#'
#'# example4: Splitting the heatmap in groups
#'rowsplit <-rep(c("Group1", "Group2","Group3", "Group4","IT","TT"), c(9,8,9,7,1,1))
#'heatmap_samples_above_median(df_harmonized, method = "median", cluster_switch = rowsplit, grouping_var = "indication")
#'
#'# example5: Reordering the heatmap and redefining the groups based on user-defined pathways
#'reorder_list <- c("Gene 14","Gene 30","Gene 20","Gene 22","Gene 29","Gene 23","Gene 1","Gene 32","Gene 33","Gene 21","Gene 11","Gene 18","Gene 8","Gene 16","Gene 2" ,"Gene 3" ,"Gene 5" ,"Gene 6","Gene 7","Gene 10","Gene 15","Gene 17","Gene 24","Gene 25","Gene 28","Gene 31","Gene 35",'Gene 19',"Gene 12","Gene 26","Gene 4","Gene 9","Immune target","Tumor target")
#'rowsplit_new <-rep(c("Pathway1", "Pathway2","Pathway3", "Pathway4","G4","G9","IT","TT"), c(9,5,13,3,1,1,1,1))
#'heatmap_samples_above_median(df_harmonized, method="median", cluster_switch=rowsplit_new, order_4_heatmap = reorder_list)
#'
#'
#'
heatmap_samples_above_median <- function(df, method=c("median","mean"),
                                         grouping_var=NA,
                                         cluster_switch=NULL,
                                         cluster_columns=TRUE,
                                         column_title=NA,
                                         plot_groups = NA,
                                         order_4_heatmap = NA){

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")


  gene_panel <- unique(df$symbol)

  df_percentage <- as.data.frame(get_perc_gene(df,gene_panel[1], method, grouping_var)) %>% select (!!sym(grouping_var))

  for (i in gene_panel){

    df_perc <- get_perc_gene(df,i, method, grouping_var)
    df_percentage <- cbind(df_percentage,df_perc[,2])
  }


  ifelse(!is.na(column_title),
         column_title <- column_title,
         column_title <- "Heatmap representing samples above median expression level of genes")


  if(is.na(plot_groups[1])) {
    df_4_plot=df_percentage
  } else {
    df_4_plot=df_percentage %>% dplyr::filter(grepl(paste0("(", paste(plot_groups, collapse = "|"), ")( pre| post)?"), !!sym(grouping_var)))
  }

  colnames(df_4_plot) <- gsub("pct_above_median_","",colnames(df_4_plot))
  colnames(df_4_plot) <- gsub("pct_above_mean_","",colnames(df_4_plot))
  colnames(df_4_plot) <- gsub("pct_above_Q3_","",colnames(df_4_plot))
  colnames(df_4_plot) <- gsub("_expression","",colnames(df_4_plot))

  mat_4_plot <- as.matrix(t(df_4_plot[,-1]))
  colnames(mat_4_plot) <- df_4_plot[,1]
  lt25 <- apply(mat_4_plot, 1, function(x) {between(x, 0, 24.99)}) %>% rowSums(na.rm = T)
  bw_25_50 <- apply(mat_4_plot, 1, function(x) {between(x, 25, 49.99)}) %>% rowSums(na.rm = T)
  bw_50_75 <- apply(mat_4_plot, 1, function(x) {between(x, 50, 74.99)}) %>% rowSums(na.rm = T)
  mt_75 <- apply(mat_4_plot, 1, function(x) {between(x, 75, 100)}) %>% rowSums(na.rm = T)

  #Annotations for the heatmap
  mat_4_barplot <- cbind(lt25, bw_25_50, bw_50_75, mt_75)
  rownames(mat_4_barplot) <- df_4_plot[,1]
  ann <- as.matrix(df_4_plot[,1])

  annotation = HeatmapAnnotation(
    pct_above_median = anno_barplot(mat_4_barplot, gp = gpar(fill = gmb_palette), bar_width = 0.5, height = unit(2, "cm")))


  ifelse(!is.na(order_4_heatmap),
         mat_4_plot <- mat_4_plot[match(order_4_heatmap, row.names(mat_4_plot)), ],
         mat_4_plot <- mat_4_plot[order(row.names(mat_4_plot)), ])

  ## Plotting the heatmap

  heatmap <- Heatmap(mat_4_plot,
                     cluster_rows = FALSE,
                     cluster_columns = cluster_columns,
                     show_column_names = TRUE,
                     row_order = sort(rownames(mat_4_plot)),
                     top_annotation = annotation,
                     column_title = column_title,
                     show_heatmap_legend = FALSE,
                     row_split = cluster_switch,
                     row_gap = unit(5, "mm"),
                     col=colorRamp2(c(0,25,50,75,100),c("darkblue","lightblue","white","#FFA070","#C71585")),
                     cell_fun = function(j, i, x, y, width, height, fill)
                     {
                       grid.text(sprintf("%.1f", mat_4_plot[i, j]), x, y, gp = gpar(fontsize = 10))
                     })
  return(heatmap)
}
