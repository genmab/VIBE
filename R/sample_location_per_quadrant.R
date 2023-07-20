#' Obtain the quadrant specific direction of selected genes and returns positioning of labels on the quadrant
#' @export
#'
#'
#' @param x the expression values of target1
#' @param y the expression values of target2
#' @param xintercept threshold for target 1 as per the selected threshold method
#' @param yintercept threshold for target 2 as per the selected threshold method
#'
#' @note @add_quadrant_info acts as a wrapper around this function
#'
#' @return the quadrant information per sample based on the selected targets to be added to the quad_summary table

sample_location_per_quadrant <- function(x, y, xintercept, yintercept) {
  z <- ifelse(x >= xintercept & y >= yintercept,
              2L, ## Q2 high, high
              ifelse(x >= xintercept & y < yintercept,
                     4L, ## Q4 high, low
                     ifelse(x < xintercept & y < yintercept,
                            3L, ## Q3 low, low
                            1L ## Q1 low, high
                     )
              )
  )
  return (z)
}
