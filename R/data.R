#' SPECT Thallium test data set
#'
#' Single-photon-emission computed-tomography (SPECT) thallium is a non-invasive diagnostic test
#' used to diagnose coronary artery disease (CAD). SPECT thallium test was performed on 2688 patients.
#' CAD is diagnosed when stenosis exceeds 50\% of the artery, as evaluated by coronary
#' angiography (gold standard). Only 471 patients underwent the coronary angiography
#' for verification of the CAD status. The rest of the patients were unverified (82.5\%).
#' @format A data frame with 2688 rows and five variables:
#' \describe{
#'   \item{SPECT thallium test, \eqn{T}:}{Binary, 1 = Positive, 0 = Negative}
#'   \item{CAD, \eqn{D}:}{Binary, 1 = Yes, 0 = No}
#'   \item{Covariates, \eqn{X}:}
#'   \enumerate{
#'     \item{Gender, \eqn{X1}: Binary, 1 = Male, 0 = Female}
#'     \item{Stress mode, \eqn{X2}: Binary, 1 = Dipyridamole (Medication for stress test when the patient is unable to exercise), 0 = Exercise}
#'     \item{Age, \eqn{X3}: Binary, 1 = 60 years and above, 0 = Below 60 years}
#'   }
#' }
#' @source
#' \enumerate{
#'   \item{Cecil, M. P., Kosinski, A. S., Jones, M. T., Taylor, A., Alazraki, N. P., Pettigrew, R. I., & Weintraub, W. S. (1996). The importance of work-up (verification) bias correction in assessing the accuracy of SPECT thallium-201 testing for the diagnosis of coronary artery disease. Journal of Clinical Epidemiology, 49(7), 735–742.}
#'   \item{Kosinski, A. S., & Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. Biometrics, 59(1), 163–171.}
#' }
"cad_pvb"
