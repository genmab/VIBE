#' @title VIBE data
#'
#' @description Randomly generated data was structured in a similar way to RNAseq datasets. This can be used to test out the VIBE functions.
#'
#' @description The data is in long format.
#'
#' @docType data
#'
#' @usage data(VIBE_data)
#'
#' @format A data frame with 56000 rows and 7 variables:
#' \describe{
#'   \item{tumor}{solid indications(BLCA, BRCA, CERV, CRC, HNSCC, NET, NSCLC, OV, PDAC, PRAD, RCC, SARC, SCLC, STAD, TNBC, UTEN)}
#'   \item{treatment_flag}{Whether the sample is collected from a pre- or post-treatment tumor.("pre", "post")}
#'   \item{gene}{Names for variables that act as a standin for gene symbols}
#'   \item{log2_tpm}{Randomly generated values that mimic log2(tpm+1) values (-4.781005 - 12.838787)}
#'   \item{analysis}{Randomly generated analysis codes, unique for each sample.}
#'   \item{patient_no}{Randomly generated patient codes consisting of 5 letters, four numbers and 1 letter. ("BHEKD7366G"-"ZTLUS4324O") }
#'   \item{database}{Patients are assigned to either database1 or database2. ("database1", "database2")}
#' }
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(VIBE_data)
"VIBE_data"
