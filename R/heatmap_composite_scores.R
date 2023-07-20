#' Get composite scores when compared between target when compared with other gene/pathway
#' @export
#'
#' @import tidyverse
#' @import dplyr
#' @import tibble
#'
#'
#' @param df dataframe with expression values for each indication
#' @param target1 Official gene symbol of gene 1
#' @param target2_list Official list of gene symbols
#' @param QOI Which of the four quadrants are of interest for the question; Q1 (Top-left; high-low); Q2 (Top-right; high-high); Q3 (bottom-left; low-high) and Q4 (bottom right; low-low)
#' @param threshold_method Threshold method to use; median or mean or Q3 (75th quantile)
#' @param grouping_var groups object based on this variable; intr if not provided
#' @param filter_threshold threshold for selecting rows relevant in terms of composite scores
#' @param reverse if the plotted categories should be reversed. The horizontal and vertical axes will be flipped as well.
#' @param order_gene the gene column can be used for sorting the axes.
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
#'# example1: obtain heatmap of composite scores for Q2 quadrant for each comparison for default grouping var "intr"
#' heatmap_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     threshold_method = "median")
#'# example2: obtain heatmap of composite scores for Q2 quadrant for each comparison for different grouping variable
#' heatmap_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     grouping_var="indication",
#'                     threshold_method = "median")
#'
#'# example3: Heatmap filtered based on threshold of composite score 30 in any of the the grouping variables
#' heatmap_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     threshold_method = "median",
#'                     filter_threshold = 30)
#'
#'# example4: Heatmap filtered based on threshold of composite score 30 in any of the the grouping variables ordered by gene
#' heatmap_composite_scores(df_harmonized, target1, target2_list,
#'                     QOI = "Q2",
#'                     threshold_method = "median",
#'                     filter_threshold = 30,
#'                     order_gene = "Gene.9")
#'
#'

heatmap_composite_scores <- function (df, target1, target2_list,
                                      QOI = c("Q1", "Q2", "Q3", "Q4","all"),
                                      threshold_method = c("median", "mean", "Q3"),
                                      grouping_var=NA,
                                      filter_threshold = NA, reverse=T, order_gene = NA){

  ifelse(!is.na(grouping_var), grouping_var <- grouping_var, grouping_var <- "intr")
  ifelse(!is.na(filter_threshold), filter_threshold <- filter_threshold, filter_threshold <- 0)


  #order <- rownames(composite_scores)

  composite_scores = get_composite_scores(df, target1 = target1, target2_list = target2_list,threshold_method=threshold_method, QOI=QOI, filter_threshold = filter_threshold)

  ifelse(!is.na(order_gene), composite_scores <- composite_scores[order(-composite_scores[[order_gene]]), ] , composite_scores <- composite_scores)

  #rownames(composite_scores) <- factor(rownames(composite_scores), levels=(rownames(composite_scores))[order(composite_scores$CD27_CD40)])
  composite_scores_long <-   composite_scores %>%
    tibble::rownames_to_column(var = "category") %>%
    pivot_longer(-category, names_to = "symbol", values_to = "value") %>%
    group_by(symbol)

  if (reverse == T){

    composite_scores_tile = ggplot(composite_scores_long, aes(x=symbol, y=factor(category, levels = factor(rev(rownames(composite_scores)))))) + # x and y axes => Var1 and Var2
      geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
      geom_text(aes( label = round(value, 2))) + # write the values
      scale_fill_gradient2(low = "white",
                           mid = "#B2edea",
                           high = "#008b8b",
                           midpoint = 45) + # determine the colour
      theme(panel.grid.major.x=element_blank(), #no gridlines
            panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.background=element_rect(fill="white"), # background=white
            axis.text.x = element_text( hjust = 0,vjust=1,size = 12,face = "bold", angle=45, colour = "black"),
            plot.title = element_text(size=20,face="bold"),
            axis.text.y = element_text(size = 12,face = "bold", colour = "black")) +
      # ggtitle("Composite Scores") +
      theme(legend.title=element_text(face="bold", size=14)) +
      scale_x_discrete(name="", position = "top")  +
      scale_y_discrete(name="") +
      labs(fill="Composite scores")

  }
  else{

    composite_scores_tile= ggplot(composite_scores_long, aes(y=symbol, x=factor(category, levels = factor(rownames(composite_scores))))) + # x and y axes => Var1 and Var2
      geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
      geom_text(aes( label = round(value, 2))) + # write the values
      scale_fill_gradient2(low = "white",
                           mid = "#B2edea",
                           high = "#008b8b",
                           midpoint = 45) + # determine the colour
      theme(panel.grid.major.x=element_blank(), #no gridlines
            panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.background=element_rect(fill="white"), # background=white
            axis.text.y = element_text( hjust = 1,vjust=1,size = 12,face = "bold", angle=0, colour = "black"),
            plot.title = element_text(size=20,face="bold"),
            axis.text.x = element_text(size = 12,face = "bold", colour = "black", angle=45, hjust=1)) +
      # ggtitle("Composite Scores") +
      theme(legend.title=element_text(face="bold", size=14)) +
      scale_x_discrete(name="")  +
      scale_y_discrete(name="") +
      labs(fill="Composite scores")
  }

  return (composite_scores_tile)
}
