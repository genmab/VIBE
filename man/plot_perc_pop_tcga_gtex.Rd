% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_perc_pop_tcga_gtex.R
\name{plot_perc_pop_tcga_gtex}
\alias{plot_perc_pop_tcga_gtex}
\title{Plot percentage of samples per quadrant or Quadrant of interest (QOI)}
\usage{
plot_perc_pop_tcga_gtex(
  df,
  target1,
  target2,
  threshold_method = c("median", "mean", "Q3"),
  QOI = c("Q1", "Q2", "Q3", "Q4"),
  order = c("T", "F"),
  dataset = NA,
  grouping_var = NA,
  plot_groups = NA,
  fill_by = "study"
)
}
\arguments{
\item{df}{the dataframe with TCGA and GTEX data.}

\item{target1}{Official gene symbol of gene 1}

\item{target2}{Official gene symbol of gene 2}

\item{threshold_method}{Threshold method to use; median or mean or Q3 (75th quantile)}

\item{QOI}{Quadrant of interest among 4 quadrants (Q1, Q2, Q3 OR Q4)}

\item{order}{True or False to sort the a axis in decreasing order ofr the quadrant of interest (Will not work for all quadrants)}

\item{dataset}{Name of your dataset for captioning}

\item{grouping_var}{column name which you want to group data by (intr if left to NA)}

\item{plot_groups}{Character vector with values you wish to include in the plot, threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.}

\item{fill_by}{column name of column that should be used to color plot bars with. Standard is set to "study"}
}
\value{
the dataframe with percentage of samples per quadrant or Quadrant of interest (QOI)
}
\description{
Plot percentage of samples per quadrant or Quadrant of interest (QOI)
}
