#' SPECT Thallium test data set
#'
#' Single-photon-emission computed-tomography (SPECT) thallium is a non-invasive diagnostic test used to diagnose coronary artery disease (CAD). SPECT thallium test was performed on 2688 patients. CAD is diagnosed when stenosis exceeds 50\% of the artery, as evaluated by coronary angiography (gold standard). Only 471 patients underwent the coronary angiography for verification of the CAD status. The rest of the patients were unverified (82.5\%).
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
#' Hepatic scintigraphy test data set
#'
#' The data set pertains to hepatic scintigraphy, a diagnostic imaging technique used for detecting liver cancer. The test was performed on 650 patients, where 344 patients were verified by liver pathological examination (gold standard test). The percentage of unverified patients is 47.1\%.
#' @format A data frame with 650 rows and three variables:
#' \describe{
#'    \item{Liver cancer, \eqn{disease}:}{Binary, 1 = Yes, 0 = No}
#'    \item{Hepatic scintigraphy, \eqn{test}:}{Binary, 1 = Positive, 0 = Negative}
#'    \item{Verified, \eqn{verified}:}{Binary, 1 = Yes, 0 = No}
#' }
#' @source
#' \enumerate{
#'   \item{Drum, D. E., & Christacopoulos, J. S. (1972). Hepatic scintigraphy in clinical decision making. Journal of Nuclear Medicine, 13(12), 908–915.}
#' }
"hepatic_pvb"
#' Diaphanography test data set
#'
#' Diaphanography test is a noninvasive method (diagnostic test) of breast examination by transillumination using visible or infrared light to detect the presence of breast cancer. The test was performed on 900 patients. Only 88 patients were verified by breast tissue biopsy for histological examination (gold standard test). The percentage of unverified patients is 90.2\%.
#' @format A data frame with 900 rows and three variables:
#' \describe{
#'    \item{Breast cancer, \eqn{disease}:}{Binary, 1 = Yes, 0 = No}
#'    \item{Diaphanography, \eqn{test}:}{Binary, 1 = Positive, 0 = Negative}
#'    \item{Verified, \eqn{verified}:}{Binary, 1 = Yes, 0 = No}
#' }
#' @source
#' \enumerate{
#'   \item{Marshall, V., Williams, D. C., & Smith, K. D. (1984). Diaphanography as a means of detecting breast cancer. Radiology, 150(2), 339–343.}
#' }
"diapha_pvb"
