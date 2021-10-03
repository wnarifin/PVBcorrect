#' PVBcorrect: A package to perform partial verification bias correction for estimates of accuracy measures in diagnostic accuracy studies
#'
#' The package contains a number of functions to perform partial verification bias
#' (PVB) correction for estimates of accuracy measures in diagnostic accuracy studies.
#' The available methods are: Begg and Greenes' method (as extended by Alonzo & Pepe, 2005),
#' Begg and Greenes' method 1 and 2 (with PPV and NPV as extended by deGroot et al, 2011),
#' multiple imputation method by logistic regression (Harel & Zhou, 2006),
#' and EM-based logistic regression method (Kosinski & Barnhart, 2003).
#'
#' @section General function:
#' \code{\link{view_table}}
#'
#' @section PVB correction main functions:
#' \code{\link{acc_cca}}, \code{\link{acc_ebg}}, \code{\link{acc_mi}}, \code{\link{acc_em}}
#'
#' @section PVB correction additional functions:
#' \code{\link{acc_dg1}}, \code{\link{acc_dg2}}
#'
#' @section Data set:
#' \code{\link{cad_pvb}}
#'
#' @references
#' \enumerate{
#'   \item{Alonzo, T. A., & Pepe, M. S. (2005). Assessing accuracy of a continuous screening test in the presence of verification bias. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 173–190.}
#'   \item{Begg, C. B., & Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. Biometrics, 207–215.}
#'   \item{de Groot, J. A. H., Janssen, K. J. M., Zwinderman, A. H., Bossuyt, P. M. M., Reitsma, J. B., & Moons, K. G. M. (2011). Correcting for partial verification bias: a comparison of methods. Annals of Epidemiology, 21(2), 139–148.}
#'   \item{Harel, O., & Zhou, X.-H. (2006). Multiple imputation for correcting verification bias. Statistics in Medicine, 25(22), 3769–3786.}
#'   \item{He, H., & McDermott, M. P. (2012). A robust method using propensity score stratification for correcting verification bias for binary tests. Biostatistics, 13(1), 32–47.}
#'   \item{Kosinski, A. S., & Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. Biometrics, 59(1), 163–171.}
#' }
#'
#' @docType package
#' @name PVBcorrect
NULL
#> NULL
