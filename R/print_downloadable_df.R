#' For HTML notebooks, print a downloadable table that also allows for searching
#' @export
#'
#' @import DT
#'
#' @param df a dataframe
#' @return a downloadable table that is nicely printed in a html notebook
#'
#' @examples
#' # example: print any dataframe to the viewer (or nicely printed in a notebook)
#' data("VIBE_data")
#' print_downloadable_df(VIBE_data[1:10,])

#creating a table that can be downloaded
print_downloadable_df = function(df) {
  datatable(
    df, extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = list('copy',list(
        extend = "collection",
        buttons=c('csv', 'excel', 'pdf'),
        text="Download",
        title=NULL)
      ),
      pageLength = 20,
      fixedHeader = TRUE
    )
  )
}
