#' IPD data from XX
#'
#' An `IPDfromKM` generated dataset based on the curves and tables of Figure 2
#' in [Thomas et al (2021) *Safety and Efficacy of the BNT162b2 mRNA Covid-19 Vaccine through 6 Months*](https://pubmed.ncbi.nlm.nih.gov/34525277/)
#'
#' @format ## `ipd_data`
#' A data frame with 44,939 and 3 columns:
#' \describe{
#'   \item{time}{day since first dose}
#'   \item{status}{value of 1 for occurrence of Covid-19 after first dose;
#'   0 censored}
#'   \item{treat}{value of 1 for treated; value of 0 for placebo}
#'   ...
#' }
#' @source <https://pubmed.ncbi.nlm.nih.gov/34525277/#&gid=article-figures&pid=figure-2-uid-1>
"ipd_data"
