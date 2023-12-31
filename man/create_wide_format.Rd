% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_wide_format.R
\name{create_wide_format}
\alias{create_wide_format}
\title{Creates a wide format of the dataframe df}
\usage{
create_wide_format(df, additional_columns = "")
}
\arguments{
\item{df}{dataframe with expression values for each indication (output of @harmonize_df )}

\item{additional_columns}{vector with column names to be kept as well (is used for grouping_var in some other vibe functions)}
}
\value{
dataframe with wide format of the data from profiler
}
\description{
Creates a wide format of the dataframe df
}
\examples{
# prepare data
data("VIBE_data")
df_harmonized = harmonize_df(df = VIBE_data,
                            col_patientid = "patient_no",
                          col_sampleid = "analysis",
                          col_indication = "tumor",
                          col_treatment = "treatment_flag",
                          col_symbol = "gene",
                          col_value = "log2_tpm",
                          unit = "tpm",
                          keep_extra_cols = c("database"))


# example 1
df_wide <- create_wide_format(df_harmonized)
dim(df_harmonized)
dim(df_wide)
dplyr::glimpse(df_wide)

# example 2 keep additional column database
df_wide_extra <- create_wide_format(df_harmonized, additional_columns = "database")
dim(df_wide_extra)
dplyr::glimpse(df_wide_extra)
}
