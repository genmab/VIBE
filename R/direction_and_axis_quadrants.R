#' Obtain the quadrant specific direction of selected genes and returns positioning of labels on the quadrant
#' @export
#'
#'
#' @param df Dataframe as harmonized by @harmonize_df
#' @param QOI Quadrant of interest among 4 quadrants
#' @param target1 Official gene symbol of gene 1
#' @param target2 Official gene symbol of gene 2
#'
#' @return a list comprising of the information of the direction of both the target (high/low) and the plotting location for each quadrant
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
#'# example 1: positioning of labels Q1 (samples low in target1, high in target2)
#'direction_and_axis_quadrants(df = df_harmonized, QOI = "Q1",
#'                             target1 = "Immune target", target2 = "Tumor target")

#'# example 2: positioning of labels Q2 (samples high in both targets)
#'direction_and_axis_quadrants(df = df_harmonized, QOI = "Q2",
#'                             target1 = "Immune target", target2 = "Tumor target")


#'# example 3: positioning of labels Q3 (samples high in target1, low in target 2)
#'direction_and_axis_quadrants(df = df_harmonized, QOI = "Q3",
#'                             target1 = "Immune target", target2 = "Tumor target")

#'# example 4: positioning of labels Q4 (samples low in both targets)
#'direction_and_axis_quadrants(df = df_harmonized, QOI = "Q4",
#'                             target1 = "Immune target", target2 = "Tumor target")
#'
direction_and_axis_quadrants <- function(df, QOI = c("Q1","Q2","Q3","Q4"), target1, target2){
  df_wide <- create_wide_format(df)
  if (QOI=="Q1"){
    dir_x ="low"
    dir_y="high"
    x=1
    y=max(df_wide[,target2])-0.5
  }
  else if(QOI=="Q2"){
    dir_x ="high"
    dir_y="high"
    x=max(df_wide[,target1])-0.5
    y=max(df_wide[,target2])-0.5
  }
  else if(QOI=="Q3"){
    dir_x ="low"
    dir_y="low"
    x=1
    y=1
  }
  else {
    dir_x ="low"
    dir_y="high"
    x=max(df_wide[,target1])-0.5
    y=1
  }
  return(c(dir_x,dir_y,x,y))
}
