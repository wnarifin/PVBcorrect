# Functions for PVB correction
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Author: Wan Nor Arifin
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# General Functions ====

#' Test vs Disease/Gold Standard cross-classification table
#'
#' @description View Test vs Disease/Gold Standard cross-classification table.
#' @param data A data frame, with at least "Test" and "Disease" variables.
#' @param test The "Test" variable name, i.e. the test result. The variable must be in binary; positive = 1, negative = 0 format.
#' @param disease The "Disease" variable name, i.e. the true disease status. The variable must be in binary; positive = 1, negative = 0 format.
#' @param show_unverified Optional. Set to \code{TRUE} to view observations with unverified disease status. The default is \code{FALSE}.
#' @param show_total Optional. Set to \code{TRUE} to view total by test result. The default is \code{FALSE}.
#' @return A cross-classification table.
#' @examples
#' str(cad_pvb)  # built-in data
#'
#' view_table(data = cad_pvb, test = "T", disease = "D")  # without unverified observations
#' view_table(data = cad_pvb, test = "T", disease = "D", show_total = TRUE)
#'   # also with total observations by test result
#'
#' view_table(data = cad_pvb, test = "T", disease = "D", show_unverified = TRUE)
#'   # with unverified observations
#' view_table(data = cad_pvb, test = "T", disease = "D", show_unverified = TRUE,
#'            show_total = TRUE)  # also with total observations by test result
#' @importFrom stats addmargins
#' @export
view_table = function(data, test, disease, show_unverified = FALSE, show_total = FALSE) {
  # assign variables
  test = data[, test]
  disease = data[, disease]

  # option: show verified/not
  if (show_unverified == TRUE) {
    tbl = table(Test = 2 - test, Disease = 2 - disease, useNA = "ifany")  # rearrange for epidemiology view
    if (ncol(tbl) < 3) {
      tbl = cbind(tbl, na = c(0, 0))
    }  # whenever there are no missing disease status
    dimnames(tbl) = list(Test = c("yes","no"),
                         Disease = c("yes","no","unverified"))
  } else if (show_unverified == FALSE) {
    tbl = table(Test = 2 - test, Disease = 2 - disease)  # rearrange for epidemiology view
    rownames(tbl) = colnames(tbl) = c("yes","no")
  }

  if (show_total == TRUE) {
    tbl = addmargins(tbl)[-3, ]
    colnames(tbl)[ncol(tbl)] = "Total"
  }

  return(tbl)
}

# Accuracy Measures Functions ====

#' Complete Case Analysis, CCA
#'
#' @description Perform Complete Case Analysis, CCA, used for complete data and multiple imputation, MI.
#' @inheritParams view_table
#' @param ci View confidence interval (CI). The default is \code{FALSE}.
#' @param ci_level Set the CI width. The default is 0.95 i.e. 95\% CI.
#' @param description Print the name of this analysis. The default is \code{TRUE}. This can be turned off for repeated analysis, for example in bootstrapped results.
#' @return A list object containing:
#' \describe{
#'   \item{acc_results}{The accuracy results.}
#' }
#' @examples
#' acc_cca(data = cad_pvb, test = "T", disease = "D")
#' acc_cca(data = cad_pvb, test = "T", disease = "D", ci = TRUE)
#' @importFrom stats qnorm
#' @export
acc_cca = function(data, test, disease,
                   ci = FALSE, ci_level = .95,
                   description = TRUE) {
  # assign variables
  test = data[, test]
  disease = data[, disease]

  # setup empty vector
  acc_values = acc_vars = acc_ses = rep(NA, 4)
  ci_values = matrix(rep(NA, 4*2), nrow = 4, ncol = 2)

  # setup table
  tbl = table(test, disease)  # imp! do not change to view_table here, acc_mi will break

  # get measures from table
  acc_values[1] = tbl[2, 2] / sum(tbl[, 2])  # Sn
  acc_values[2] = tbl[1, 1] / sum(tbl[, 1])  # Sp
  acc_values[3] = tbl[2, 2] / sum(tbl[2, ])  # PPV
  acc_values[4] = tbl[1, 1] / sum(tbl[1, ])  # NPV

  # only calculate when ci is requested
  if (ci == TRUE) {
    # get variances, var formula ref: Zhou (1994)
    acc_vars[1] = tbl[2, 2] * tbl[1, 2] / sum(tbl[, 2])^3  # Var Sn
    acc_vars[2] = tbl[1, 1] * tbl[2, 1] / sum(tbl[, 1])^3  # Var Sp
    acc_vars[3] = tbl[2, 2] * tbl[2, 1] / sum(tbl[2, ])^3  # Var PPV
    acc_vars[4] = tbl[1, 1] * tbl[1, 2] / sum(tbl[1, ])^3  # Var NPV

    # calculate SEs
    acc_ses = sqrt(acc_vars)

    # ci, default 95% ci
    for (i in 1:4) {
      ci_values[i, ] = acc_values[i] + c(-1, 1) * qnorm(1 - (1 - ci_level)/2) * acc_ses[i]
    }
  }

  # output
  df = data.frame(Est = acc_values,
                  row.names = c("Sn", "Sp", "PPV", "NPV"))
  # when ci is requested
  if (ci == TRUE) {
    df_ci = data.frame(SE = acc_ses, LowCI = ci_values[, 1], UppCI = ci_values[, 2])
    df = cbind(df, df_ci)
  }
  # when ci not requested
  else {
    df = df
  }
  # allows turning method description off for iterative procedure, e.g. MI
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nUncorrected for PVB: Complete Case Analysis\n\n")
  }
  acc_cca_list = list(acc_results = df)
  return(acc_cca_list)
}

#' PVB correction by extended Begg and Greenes' method
#'
#' @description Perform PVB correction by Begg and Greenes' method (as extended by Alonzo & Pepe, 2005).
#' @inheritParams acc_cca
#' @param covariate The name(s) of covariate(s), i.e. other variables associated with either test or disease status.
#'   Specify as name vector, e.g. c("X1", "X2") for two or more variables. The variables must be in formats acceptable to GLM.
#' @param saturated_model Set as \code{TRUE} to obtain the original Begg and Greenes' (1983) when all possible interactions are included.
#' @param show_fit Set to \code{TRUE} to view model fit summary for the logistic regression model.
#' @param show_boot Set to \code{TRUE} to show bootstrap iterations.
#' @param ci_type Set confidence interval (CI) type. Acceptable types are "norm", "basic", "perc", and "bca",
#'   for bootstrapped CI. See \code{\link[boot]{boot.ci}} for details.
#' @param R The number of bootstrap samples. Default \code{R = 999}.
#' @param seednum Set the seed number for the bootstrapped CI. The default is not set, so it depends on the user
#'   to set it outside or inside the function.
#' @param r_print_freq Print the current bootstrap sample number at each specified interval.
#'   Default \code{r_print_freq = 100}.
#' @return A list object containing:
#' \describe{
#'   \item{boot_data}{An object of class "boot" from \code{\link[boot]{boot}}.
#'   Contains Sensitivity, Specificity, PPV, and NPV}
#'   \item{boot_ci_data}{A list of objects of type "bootci" from \code{\link[boot]{boot.ci}}.
#'   Contains Sensitivity, Specificity, PPV, NPV.}
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Alonzo, T. A., & Pepe, M. S. (2005). Assessing accuracy of a continuous screening test in the presence of verification bias. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 173–190.}
#'   \item{Begg, C. B., & Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. Biometrics, 207–215.}
#'   \item{He, H., & McDermott, M. P. (2012). A robust method using propensity score stratification for correcting verification bias for binary tests. Biostatistics, 13(1), 32–47.}
#' }
#' @examples
#' # point estimates
#' acc_ebg(data = cad_pvb, test = "T", disease = "D")
#' acc_ebg(data = cad_pvb, test = "T", disease = "D", covariate = "X3")
#' acc_ebg(data = cad_pvb, test = "T", disease = "D", covariate = "X3", saturated_model = TRUE)
#'
#' # with bootstrapped confidence interval
#' acc_ebg(data = cad_pvb, test = "T", disease = "D", ci = TRUE, seednum = 12345)
#' acc_ebg(data = cad_pvb, test = "T", disease = "D", covariate = "X3", ci = TRUE, seednum = 12345)
#' acc_ebg(data = cad_pvb, test = "T", disease = "D", covariate = "X3", saturated_model = TRUE,
#'         ci = TRUE, seednum = 12345)
#' @export
acc_ebg = function(data, test, disease, covariate = NULL, saturated_model = FALSE,
                   ci = FALSE, ci_level = .95, ci_type = "basic",
                   R = 999, seednum = NULL,
                   show_fit = FALSE, show_boot = FALSE, r_print_freq = 100,
                   description = TRUE) {
  # w/out ci
  if (ci == FALSE) {
    acc_ebg_point(data, test, disease, covariate, saturated_model,
                  show_fit, show_boot = FALSE, description = description)
  }

  # w ci
  else {
    acc_ebg_boot_ci(data, test, disease, covariate, saturated_model,
                    show_fit, show_boot, description = description,
                    ci_level = ci_level, ci_type = ci_type,
                    R = R, seednum = seednum, r_print_freq = r_print_freq)
  }
}

#' PVB correction by Begg and Greenes' method with asymptotic normal CI
#'
#' @description PVB correction by Begg and Greenes' method with asymptotic normal CI. This is limited to no covariate.
#' @inheritParams acc_ebg
#' @param ci View confidence interval (CI). The default is \code{FALSE}.
#' @return A list object containing:
#' \describe{
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Begg, C. B., & Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. Biometrics, 207–215.}
#'   \item{Harel, O., & Zhou, X.-H. (2006). Multiple imputation for correcting verification bias. Statistics in Medicine, 25(22), 3769–3786.}
#'   \item{Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and specificity corrected for verification bias. Communications in Statistics-Theory and Methods, 22(11), 3177–3198.}
#'   \item{Zhou, X.-H. (1994). Effect of verification bias on positive and negative predictive values. Statistics in Medicine, 13(17), 1737–1745.}
#'   \item{Zhou, X.-H., Obuchowski, N. A., & McClish, D. K. (2011). Statistical Methods in Diagnostic Medicine (2nd ed.). John Wiley & Sons.}
#'
#' }
#' @examples
#' acc_bg(data = cad_pvb, test = "T", disease = "D")  # equivalent to result by acc_ebg()
#' acc_bg(data = cad_pvb, test = "T", disease = "D", ci = TRUE)
#'   # the CIs are slightly differerent from result by acc_ebg()
#' @importFrom stats addmargins
#' @importFrom stats qnorm
#' @export
acc_bg = function(data, test, disease,
                  ci = FALSE, ci_level = .95,
                  description = TRUE) {

  # setup empty vector
  acc_values = acc_vars = acc_ses = rep(NA, 4)
  ci_values = matrix(rep(NA, 4*2), nrow = 4, ncol = 2)

  # setup table
  tbl = view_table(data, test, disease, show_unverified = TRUE)
  tbl = addmargins(tbl)

  # get measures from table
  s1 = tbl[1, 1]; s0 = tbl[2, 1]
  r1 = tbl[1, 2]; r0 = tbl[2, 2]
  n1 = tbl[1, 4]; n0 = tbl[2, 4]; n = tbl[3, 4]

  acc_values[1] = n1*s1/(s1+r1) / (n1*s1/(s1+r1) + n0*s0/(s0+r0))  # Sn
  acc_values[2] = n0*r0/(s0+r0) / (n1*r1/(s1+r1) + n0*r0/(s0+r0))  # Sp
  acc_values[3] = s1 / (s1 + r1)  # PPV
  acc_values[4] = r0 / (s0 + r0)  # NPV
  acc_values

  # only calculate when ci is requested
  if (ci == TRUE) {
    # get variances, var formula ref: Zhou (1994)
    acc_vars[1] = (acc_values[1]*(1-acc_values[1]))^2 * (n/(n0*n1) + r1/(s1*(s1+r1)) + r0/(s0*(s0+r0)))  # Var Sn
    acc_vars[2] = (acc_values[2]*(1-acc_values[2]))^2 * (n/(n0*n1) + s1/(r1*(s1+r1)) + s0/(r0*(s0+r0)))  # Var Sp
    acc_vars[3] = s1*r1 / (s1+r1)^3  # Var PPV
    acc_vars[4] = s0*r0 / (s0+r0)^3  # Var NPV

    # calculate SEs
    acc_ses = sqrt(acc_vars)

    # ci, default 95% ci
    for (i in 1:4) {
      ci_values[i, ] = acc_values[i] + c(-1, 1) * qnorm(1 - (1 - ci_level)/2) * acc_ses[i]
    }
  }

  # output
  df = data.frame(Est = acc_values,
                  row.names = c("Sn", "Sp", "PPV", "NPV"))
  # when ci is requested
  if (ci == TRUE) {
    df_ci = data.frame(SE = acc_ses, LowCI = ci_values[, 1], UppCI = ci_values[, 2])
    df = cbind(df, df_ci)
  }
  # when ci not requested
  else {
    df = df
  }
  # allows turning method description off for iterative procedure, e.g. MI
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Begg and Greenes' Method\n\n")
  }
  acc_bg_list = list(acc_results = df)
  return(acc_bg_list)
}

#' PVB correction by Begg and Greenes' method 1 (deGroot et al, no covariate)
#'
#' @description Perform PVB correction by Begg and Greenes' method 1 as described in deGroot et al (2011),
#'   in which it also includes PPV and NPV calculation.
#' @inheritParams acc_ebg
#' @return A data frame object containing the accuracy results.
#' @references
#' \enumerate{
#'   \item{de Groot, J. A. H., Janssen, K. J. M., Zwinderman, A. H., Bossuyt, P. M. M., Reitsma, J. B., & Moons, K. G. M. (2011). Correcting for partial verification bias: a comparison of methods. Annals of Epidemiology, 21(2), 139–148.}
#' }
#' @examples
#' acc_dg1(data = cad_pvb, test = "T", disease = "D")  # equivalent to result by acc_ebg()
#' @export
acc_dg1 = function(data, test, disease,
                   description = TRUE) {
  # setup table counts
  tbl = view_table(data, test, disease, show_unverified = TRUE)
  a = tbl[1, 1]; b = tbl[1, 2]; c = tbl[2, 1]; d = tbl[2, 2]
  t01 = tbl[1, 3]; t00 = tbl[2, 3]
  a_ = a / (a + b) * t01
  b_ = b / (a + b) * t01
  c_ = c / (c + d) * t00
  d_ = d / (c + d) * t00

  # calculate acc
  acc_dg1 = rep(NA, 4)
  acc_dg1[1] = (a + a_) / (a + a_ + c + c_)
  acc_dg1[2] = (d + d_) / (d + d_ + b + b_)
  acc_dg1[3] = (a + a_) / (a + a_ + b + b_)
  acc_dg1[4] = (d + d_) / (d + d_ + c + c_)

  #output
  acc_dg1 = data.frame(Est = acc_dg1,
                      row.names = c("Sn", "Sp", "PPV", "NPV"))
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Begg and Greenes' Method 1 (deGroot, 2011)\n\n")
  }
  return(acc_dg1)
}

#' PVB correction by Begg and Greenes' method 2 (deGroot et al, one covariate)
#'
#' @description Perform PVB correction by Begg and Greenes' method 2 as described in deGroot et al (2011),
#'   in which it also includes PPV and NPV calculation. This is limited to only one covariate.
#' @inheritParams acc_ebg
#' @return A data frame object containing the accuracy results.
#' @references
#' \enumerate{
#'   \item{de Groot, J. A. H., Janssen, K. J. M., Zwinderman, A. H., Bossuyt, P. M. M., Reitsma, J. B., & Moons, K. G. M. (2011). Correcting for partial verification bias: a comparison of methods. Annals of Epidemiology, 21(2), 139–148.}
#' }
#' @examples
#' acc_dg2(data = cad_pvb, test = "T", disease = "D", covariate = "X3")
#'   # equivalent to acc_ebg(), saturated_model
#' @export
acc_dg2 = function(data, test, disease, covariate,
                   description = TRUE) {
  # setup table counts
  tbl2 = by(data, data[, test], function(x) view_table(data = x, test = covariate, disease, show_unverified = TRUE))
  a = tbl2$`1`[1, 1]; b = tbl2$`1`[1, 2]
  c = tbl2$`1`[2, 1]; d = tbl2$`1`[2, 2]
  e = tbl2$`0`[1, 1]; f = tbl2$`0`[1, 2]
  g = tbl2$`0`[2, 1]; h = tbl2$`0`[2, 2]
  t011 = tbl2$`1`[1, 3]
  t010 = tbl2$`1`[2, 3]
  t001 = tbl2$`0`[1, 3]
  t000 = tbl2$`0`[2, 3]
  a_ = a / (a + b) * t011
  b_ = b / (a + b) * t011
  c_ = c / (c + d) * t010
  d_ = d / (c + d) * t010
  e_ = e / (e + f) * t001
  f_ = f / (e + f) * t001
  g_ = g / (g + h) * t000
  h_ = h / (g + h) * t000

  # calculate acc
  acc_dg2 = rep(NA, 4)
  acc_dg2[1] = (a + a_ + c + c_) / (a + a_ + e + e_ + c + c_ + g + g_)
  acc_dg2[2] = (f + f_ + h + h_) / (b + b_ + d + d_ + f + f_ + h + h_)
  acc_dg2[3] = (a + a_ + c + c_) / (a + a_ + c + c_ + b + b_ + d + d_)
  acc_dg2[4] = (f + f_ + h + h_) / (e + e_ + g + g_ + f + f_ + h + h_)

  # output
  acc_dg2 = data.frame(Est = acc_dg2,
                       row.names = c("Sn", "Sp", "PPV", "NPV"))
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Begg and Greenes' Method 2 (deGroot, 2011)\n\n")
  }
  return(acc_dg2)
}

# IPB method
# Arifin & Yusof, 2022; Nahorniak et al., 2015
# IPB sampling for PVB correction is based on general IPB by Nahorniak et al, 2015
#' PVB correction by inverse probability bootstrap sampling (IPB)
#'
#' @description Perform PVB correction by inverse probability bootstrap sampling.
#' @inheritParams acc_ebg
#' @param b The number of bootstrap samples, b.
#' @param option 1 = IPW weight, 2 = W_h weight, described in Arifin (2023), modified weight of Krautenbacher (2017).
#'   The default is \code{option = 1}. For small weights, \code{option = 2} is more stable (Arifin, 2023).
#' @param ci_type Set confidence interval (CI) type. Acceptable types are "norm", "basic", "perc", and "bca",
#'   for bootstrapped CI.
#' @param return_data Return data for the bootstrapped samples.
#' @param return_detail Return accuracy measures for each of the bootstrapped samples.
#' @return A list object containing:
#' \describe{
#'   \item{data_each_sample}{Raw data for each bootstrap sample, available with \code{return_data = TRUE}}
#'   \item{acc_each_sample}{Accuracy results for each bootstrap sample, available with \code{return_detail = TRUE}}
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Arifin, W. N., & Yusof, U. K. (2022). Partial Verification Bias Correction Using Inverse Probability Bootstrap Sampling for Binary Diagnostic Tests. Diagnostics, 12(11), 2839.}
#'   \item{Arifin, W. N. (2023). Partial verification bias correction in diagnostic accuracy studies using propensity score-based methods (PhD thesis, Universiti Sains Malaysia). https://erepo.usm.my/handle/123456789/19184}
#'   \item{Krautenbacher, N., Theis, F. J., & Fuchs, C. (2017). Correcting Classifiers for Sample Selection Bias in Two-Phase Case-Control Studies. Computational and Mathematical Methods in Medicine, 2017, 1–18. https://doi.org/10.1155/2017/7847531}
#'   \item{Nahorniak, M., Larsen, D. P., Volk, C., & Jordan, C. E. (2015). Using inverse probability bootstrap sampling to eliminate sample induced bias in model based analysis of unequal probability samples. PLoS One, 10(6), e0131765.}
#' }
#' @examples
#' # no covariate
#' acc_ipb(data = cad_pvb, test = "T", disease = "D", b = 1000, seednum = 12345)
#'
#' # with three covariates
#' acc_ipb(data = cad_pvb, test = "T", disease = "D", covariate = c("X1","X2","X3"),
#'         b = 1000, seednum = 12345)
#'
#' # with confidence interval
#' acc_ipb(data = cad_pvb, test = "T", disease = "D", ci = TRUE, b = 1000, seednum = 12345)
#' @export
acc_ipb = function(data, test, disease, covariate = NULL, saturated_model = FALSE, option = 2,
                   ci = FALSE, ci_level = .95, ci_type = "norm",
                   b = 1000, seednum = NULL,
                   return_data = FALSE, return_detail = FALSE,
                   description = TRUE) {
  # data = original data
  # add verification status
  data$verified = 1  # verified: 1 = yes, 0 = no
  data[is.na(data[, disease]), "verified"] = 0

  # setup empty vector
  acc_values = rep(NA, 4)  # Sn Sp only

  # gen ps/pi
  # P(V = 1 | T = t)
  # if no covariate
  if (is.null(covariate)) {
    fmula = reformulate(test, "verified")
    fit = glm(fmula, data = data, family = "binomial")
  }
  # if covariate names supplied as vector, not NULL
  else {
    if (saturated_model == FALSE) {
      fmula = reformulate(c(test, covariate), "verified")
    } else {
      fmula = reformulate(paste(c(test, covariate), collapse = " * "), "verified")
    }
    fit = glm(fmula, data = data, family = "binomial")
  }
  data$ps = predict(fit, type = "response")

  # transform ps to weight
  if (option == 1) {
    # IPW weight
    data$ps_w = ifelse(data[, "verified"] == 1, 1/data$ps, 1/(1-data$ps))
  }
  if (option == 2) {
    # W_h weight, Krautenbacher 2017
    data$ps_w = ifelse(data[, "verified"] == 1, max(unique(data$ps)) / data$ps,
                       max(unique((1-data$ps))) / (1-data$ps))
  }

  # data1 = complete case
  data1 = na.omit(data)
  n_data = nrow(data1)
  data1$p_ipb = data1$ps_w / sum(data1$ps_w)

  # resample
  data_resample_list = NULL
  set.seed(seednum)
  i = 1
  while (i < b + 1) {
    # resample
    sampled = sample(1:n_data, n_data, TRUE, prob = data1$p_ipb)  # bootstrap w weight
    data_resample = data1[sampled, ]
    # exclude invalid samples: get the dimension sum, should be 4 for 2x2 epid table
    sum_tbl = sum(dim(table(data_resample[, test], data_resample[, disease])))
    # exclude invalid sample with sum < 4
    if (sum_tbl < 4) {
      # delete invalid sample
      data_resample_list[i] = NULL
      i = i
    } else {
      # save valid sample in list
      data_resample_list[i] = list(data_resample)
      i = i + 1
    }
  }

  acc_resample_list = t(sapply(data_resample_list, acc,
                               test = test, disease = disease))

  acc_values = apply(acc_resample_list, 2, mean)
  acc_ses = apply(acc_resample_list, 2, sd)

  # check ci_type, allowed "norm", "perc"
  ci_type_allowed = c("norm", "perc")
  # only 1 argument allowed
  if (length(ci_type) > 1) {
    stop("Invalid 'ci_type' argument.")
  }
  # if argument == 1, but type invalid
  else {
    if(!any(ci_type == ci_type_allowed)) {
      stop("Invalid 'ci_type' argument.")
    }
  }

  if (ci == TRUE) {
    if (ci_type == "norm") {
      z = qnorm((1 - (1 - ci_level)/2))
      acc_ll = acc_values - z * acc_ses
      acc_ul = acc_values + z * acc_ses
      acc_ci = cbind(acc_values, acc_ses, acc_ll, acc_ul)
      dimnames(acc_ci) = list(c("Sn", "Sp", "PPV", "NPV"), c("Est", "SE", "LowCI", "UppCI"))
    } else { # perc
      sn_percentile = quantile(acc_resample_list[, 1], probs = c(0.025, 0.975))
      sp_percentile = quantile(acc_resample_list[, 2], probs = c(0.025, 0.975))
      ppv_percentile = quantile(acc_resample_list[, 3], probs = c(0.025, 0.975))
      npv_percentile = quantile(acc_resample_list[, 4], probs = c(0.025, 0.975))
      acc_ci = cbind(acc_values, acc_ses, rbind(sn_percentile, sp_percentile,
                                                ppv_percentile, npv_percentile))
      dimnames(acc_ci) = list(c("Sn", "Sp", "PPV", "NPV"), c("Est", "SE", "LowCI", "UppCI"))
    }
    acc_ipb_list = list(acc_results = acc_ci)
  } else {
    acc_point = data.frame(Est = acc_values,
                           row.names = c("Sn", "Sp", "PPV", "NPV"))
    acc_ipb_list = list(acc_results = acc_point)
  }

  # output
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Inverse probability bootstrap sampling Method\n\n")
  }
  if (return_data == TRUE) {
    acc_ipb_list = append(acc_ipb_list, list(data_each_sample = data_resample_list))
  }
  if (return_detail == TRUE) {
    acc_ipb_list = append(acc_ipb_list, list(acc_each_sample = acc_resample_list))
  }
  return(acc_ipb_list)
}

# SIPW method
# Arifin, 2023, 2025
# SIPW is an extension of IPB method for PVB correction
#' PVB correction by scaled inverse probability weighted resampling (SIPW)
#'
#' @description Perform PVB correction by scaled inverse probability weighted resampling.
#' @inheritParams acc_ebg
#' @param b The number of repeated samples, b.
#' @param option 1 = IPW weight, 2 = W_h weight, described in Arifin (2023), modified weight of Krautenbacher (2017).
#'   The default is \code{option = 2}, which is more stable for small weights (Arifin, 2023).
#' @param return_data Return data for the bootstrapped samples.
#' @param return_detail Return accuracy measures for each of the bootstrapped samples.
#' @return A list object containing:
#' \describe{
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Arifin, W. N., & Yusof, U. K. (2025). Partial Verification Bias Correction Using Scaled Inverse Probability Resampling for Binary Diagnostic Tests. medRxiv. https://doi.org/10.1101/2025.03.09.25323631}
#'   \item{Arifin, W. N., & Yusof, U. K. (2022). Partial Verification Bias Correction Using Inverse Probability Bootstrap Sampling for Binary Diagnostic Tests. Diagnostics, 12(11), 2839.}
#'   \item{Arifin, W. N. (2023). Partial verification bias correction in diagnostic accuracy studies using propensity score-based methods (PhD thesis, Universiti Sains Malaysia). https://erepo.usm.my/handle/123456789/19184}
#'   \item{Krautenbacher, N., Theis, F. J., & Fuchs, C. (2017). Correcting Classifiers for Sample Selection Bias in Two-Phase Case-Control Studies. Computational and Mathematical Methods in Medicine, 2017, 1–18. https://doi.org/10.1155/2017/7847531}
#'   \item{Nahorniak, M., Larsen, D. P., Volk, C., & Jordan, C. E. (2015). Using inverse probability bootstrap sampling to eliminate sample induced bias in model based analysis of unequal probability samples. PLoS One, 10(6), e0131765.}
#' }
#' @examples
#' # point estimates
#' acc_sipw(data = cad_pvb, test = "T", disease = "D", b = 100, seednum = 12345)
#' acc_sipw(data = cad_pvb, test = "T", disease = "D", covariate = c("X1","X2","X3"),
#'          b = 100, seednum = 12345)
#'
#' # with bootstrapped confidence interval
#' acc_sipw(data = cad_pvb, test = "T", disease = "D", ci = TRUE, ci_type = "perc",
#'          b = 100, R = 9, seednum = 12345)
#' @export
acc_sipw = function(data, test, disease, covariate = NULL, saturated_model = FALSE, option = 2,
                    ci = FALSE, ci_level = .95, ci_type = "basic",
                    b = 1000, R = 999, seednum = NULL,
                    return_data = FALSE, return_detail = FALSE,
                    #show_boot = FALSE, r_print_freq = 100, # disabled for now
                    description = TRUE) {
  # w/out ci
  if (ci == FALSE) {
    acc_sipw_point(data, test, disease, covariate, saturated_model,
                   option = option,
                   #show_boot = FALSE, # disabled for now
                   b = b, seednum = seednum, return_data = return_data,
                   description = description)
  }

  # w ci
  else {
    acc_sipw_boot_ci(data, test, disease, covariate, saturated_model,
                     option = option,
                     #show_boot = show_boot, r_print_freq = r_print_freq, # disabled for now
                     ci_level = ci_level, ci_type = ci_type,
                     R = R, b = b, seednum = seednum,
                     description = description)
  }
}

# SIPW-B method
# Arifin, 2023, 2025
# SIPW-B is an extension of IPB method for PVB correction with balanced grouping
#' PVB correction by scaled inverse probability weighted balanced resampling (SIPW-B). SIPW-B only gives results
#' for Sensitivity and Specificity, for PPV and NPV please use SIPW instead.
#'
#' @description Perform PVB correction by scaled inverse probability weighted balanced resampling.
#' @inheritParams acc_ebg
#' @param b The number of repeated samples, b.
#' @param option 1 = IPW weight, 2 = W_h weight, described in Arifin (2023), modified weight of Krautenbacher (2017).
#'   The default is \code{option = 2}, which is more stable for small weights (Arifin, 2023).
#' @param rel_size ratio control:case, D=0:D=1. The default is 1.
#' @param return_data Return data for the bootstrapped samples.
#' @param return_detail Return accuracy measures for each of the bootstrapped samples.
#' @return A list object containing:
#' \describe{
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Arifin, W. N., & Yusof, U. K. (2025). Partial Verification Bias Correction Using Scaled Inverse Probability Resampling for Binary Diagnostic Tests. medRxiv. https://doi.org/10.1101/2025.03.09.25323631}
#'   \item{Arifin, W. N., & Yusof, U. K. (2022). Partial Verification Bias Correction Using Inverse Probability Bootstrap Sampling for Binary Diagnostic Tests. Diagnostics, 12(11), 2839.}
#'   \item{Arifin, W. N. (2023). Partial verification bias correction in diagnostic accuracy studies using propensity score-based methods (PhD thesis, Universiti Sains Malaysia). https://erepo.usm.my/handle/123456789/19184}
#'   \item{Krautenbacher, N., Theis, F. J., & Fuchs, C. (2017). Correcting Classifiers for Sample Selection Bias in Two-Phase Case-Control Studies. Computational and Mathematical Methods in Medicine, 2017, 1–18. https://doi.org/10.1155/2017/7847531}
#'   \item{Nahorniak, M., Larsen, D. P., Volk, C., & Jordan, C. E. (2015). Using inverse probability bootstrap sampling to eliminate sample induced bias in model based analysis of unequal probability samples. PLoS One, 10(6), e0131765.}
#' }
#' @examples
#' # point estimates
#' acc_sipwb(data = cad_pvb, test = "T", disease = "D", b = 100, seednum = 12345)
#' acc_sipwb(data = cad_pvb, test = "T", disease = "D", covariate = c("X1","X2","X3"),
#'          b = 100, seednum = 12345)
#'
#' # with bootstrapped confidence interval
#' acc_sipwb(data = cad_pvb, test = "T", disease = "D", ci = TRUE, ci_type = "perc",
#'          b = 100, R = 9, seednum = 12345)
#' @export
acc_sipwb = function(data, test, disease, covariate = NULL, saturated_model = FALSE, option = 2,
                     rel_size = 1, # ratio control:case, D=0:D=1
                     ci = FALSE, ci_level = .95, ci_type = "basic",
                     b = 1000, R = 999, seednum = NULL,
                     return_data = FALSE, return_detail = FALSE,
                     #show_boot = FALSE, r_print_freq = 100, # disabled for now
                     description = TRUE) {
  # w/out ci
  if (ci == FALSE) {
    acc_sipwb_point(data, test, disease, covariate, saturated_model,
                    option = option,
                    rel_size = rel_size, # ratio control:case, D=0:D=1
                    #show_boot = FALSE, # disabled for now
                    b = b, seednum = seednum, return_data = return_data,
                    description = description)
  }

  # w ci
  else {
    acc_sipwb_boot_ci(data, test, disease, covariate, saturated_model,
                      option = option,
                      rel_size = rel_size, # ratio control:case, D=0:D=1
                      #show_boot = show_boot, r_print_freq = r_print_freq, # disabled for now
                      ci_level = ci_level, ci_type = ci_type,
                      R = R, b = b, seednum = seednum,
                      description = description)
  }
}

# MI Method
# Harel & Zhou, 2006
#' PVB correction by multiple imputation
#'
#' @description Perform PVB correction by multiple imputation.
#' @inheritParams acc_ebg
#' @param m The number of imputation, m.
#' @param method Imputation method. The default is "logreg". Other allowed methods are
#'   "logreg.boot", "pmm", "midastouch", "sample", "cart", "rf".
#'   See \code{\link[mice]{mice}} for details of these methods.
#' @param mi_print Print multiple imputation history on console.
#'   This is \code{\link[mice]{mice}} \code{print} option. The default is \code{FALSE}.
#' @return A list object containing:
#' \describe{
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Harel, O., & Zhou, X.-H. (2006). Multiple imputation for correcting verification bias. Statistics in Medicine, 25(22), 3769–3786.}
#' }
#' @examples
#' # no covariate
#' acc_mi(data = cad_pvb, test = "T", disease = "D", ci = TRUE, seednum = 12345, m = 10)
#'
#' # with other imputation method. e.g. random forest "rf"
#' acc_mi(data = cad_pvb, test = "T", disease = "D", ci = TRUE, seednum = 12345, m = 10,
#'        method = "rf")
#'
#' # with three covariates
#' acc_mi(data = cad_pvb, test = "T", disease = "D", covariate = c("X1", "X2", "X3"),
#'        ci = TRUE, seednum = 12345, m = 10)
#' @importFrom stats qt
#' @export
acc_mi = function(data, test, disease, covariate = NULL,
                  ci = FALSE, ci_level = .95,
                  m = 100, seednum = NA,
                  method = "logreg",
                  mi_print = FALSE,
                  description = TRUE) {
  # check method
  method_allowed = c("logreg", "logreg.boot", "pmm", "midastouch", "sample", "cart", "rf")
  # only 1 argument allowed
  if (length(method) > 1) {
    stop("Invalid 'method' argument.")
  }
  # if argument == 1, but type invalid
  else {
    if(!any(method == method_allowed)) {
      stop("Invalid 'method' argument.")
    }
  }

  # setup subset of data
  data_mi = data[, c(test, disease, covariate)]
  data_mi[, disease] = as.factor(data_mi[, disease])
  data_mi_method = mice::mice(data_mi, m = m, method = method, seed = seednum, print = mi_print)
  data_mi_method_data = mice::complete(data_mi_method, "all")

  # w/out ci
  if (ci == FALSE) {
    acc_mi_data = lapply(data_mi_method_data, acc_cca, test = test, disease = disease, description = FALSE)
    acc_mi_data = as.data.frame(acc_mi_data)

    # output
    acc_mi_out = list(acc_results = data.frame(Est = rowMeans(acc_mi_data)))
    if (description == TRUE) {
      cat("Estimates of accuracy measures\nCorrected for PVB: Multiple Imputation Method\n\n")
    }
    return(acc_mi_out)
  }

  # w ci
  else {
    # Compile est & var from each imputed dataset
    acc_mi_data = lapply(data_mi_method_data, acc_cca, test = test, disease = disease, ci = TRUE, description = FALSE)
    acc_mi_data_est = as.matrix(as.data.frame(lapply(acc_mi_data, function(i) subset(i$acc_results, select = Est))))
    acc_mi_data_var = as.matrix(as.data.frame(lapply(acc_mi_data, function(i) subset(i$acc_results, select = SE)))^2)  # turn to Variance
    # when acc_cca returns data.frame
    # acc_mi_data_est = as.matrix(as.data.frame(lapply(acc_mi_data, subset, select = Est)))
    # acc_mi_data_var = as.matrix(as.data.frame(lapply(acc_mi_data, subset, select = SE))^2)  # turn to Variance

    # Pool in mice
    acc_mi_pool = list()
    for (i in 1:4) {
      acc_mi_pool[[i]] = mice::pool.scalar(Q = acc_mi_data_est[i, ], U = acc_mi_data_var[i, ])
    }

    # Obtain 'ingredients' to construct ci according to Rubin's rule
    acc_mi_out_qbar = unlist(lapply(acc_mi_pool, function(list) list$qbar))  # Est
    acc_mi_out_ubar = unlist(lapply(acc_mi_pool, function(list) list$ubar))  # w/in-imp Var
    acc_mi_out_b = unlist(lapply(acc_mi_pool, function(list) list$b))  # between-imp Var
    acc_mi_out_t = unlist(lapply(acc_mi_pool, function(list) list$t))  # total Var

    v_sn = (m-1)*(1+(acc_mi_out_ubar/((1+1/m)*acc_mi_out_b)))^2  # v_old / df
    t_sn = qt((1-(1-ci_level)/2), v_sn)  # t-stat for ci
    acc_mi_out_se = sqrt(acc_mi_out_t)  # SE
    acc_mi_out_ci_ll = acc_mi_out_qbar + -1 * t_sn * acc_mi_out_se  # LL ci
    acc_mi_out_ci_ul = acc_mi_out_qbar + 1 * t_sn * acc_mi_out_se  # UL ci

    # output
    acc_mi_out = list(acc_results = data.frame(Est = acc_mi_out_qbar, SE = acc_mi_out_se,
                                               LowCI = acc_mi_out_ci_ll, UppCI = acc_mi_out_ci_ul,
                                               row.names = c("Sn", "Sp", "PPV", "NPV")))
    if (description == TRUE) {
      cat("Estimates of accuracy measures\nCorrected for PVB: Multiple Imputation Method\n\n")
    }
    return(acc_mi_out)
  }
}

#' PVB correction by EM-based logistic regression method
#'
#' @description Perform PVB correction by EM-based logistic regression method.
#' @inheritParams acc_ebg
#' @param mnar The default is assuming missing not at random (MNAR) missing data mechanism, \code{MNAR = TRUE}.
#'   Set this to \code{FALSE} to obtain results assuming missing at random (MAR) missing data mechanism.
#'   This will be equivalent to using \code{\link{acc_ebg}}.
#' @param show_t Print the current EM iteration number t. The default is \code{TRUE}.
#' @param t_max The maximum iteration number for EM. Default \code{t_max = 500}. It is recommended to increase the
#'   number when covariates are included.
#' @param cutoff The cutoff value for the minimum change between iteration.
#'   This defines the convergence of the EM procedure. Default \code{cutoff = 0.0001}. This can be set to a larger value
#'   to test the procedure.
#' @param t_print_freq Print the current EM iteration number t at each specified interval.
#'   Default \code{t_print_freq = 100}.
#' @param return_t Return the final EM iteration number t.
#'   This can be used for the purpose of checking the EM convergence. The default is \code{FALSE}, but is set to
#'   TRUE when \code{ci = TRUE}.
#' @return A list object containing:
#' \describe{
#'   \item{boot_data}{An object of class "boot" from \code{\link[boot]{boot}}.
#'   Contains Sensitivity, Specificity, PPV, NPV and t (i.e. EM iteration taken for convergence).
#'   Use acc_em_object$boot_data$t[,5] to check the t.}
#'   \item{boot_ci_data}{A list of objects of type "bootci" from from \code{\link[boot]{boot.ci}}.
#'   Contains Sensitivity, Specificity, PPV, and NPV.}
#'   \item{acc_results}{The accuracy results.}
#' }
#' @references
#' \enumerate{
#'   \item{Kosinski, A. S., & Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. Biometrics, 59(1), 163–171.}
#' }
#' @example long_example/acc_em_ex.R
#' @export
acc_em = function(data, test, disease, covariate = NULL, mnar = TRUE,
                  ci = FALSE, ci_level = .95, ci_type = "basic",
                  R = 999, seednum = NULL,
                  show_t = TRUE, t_max = 500, cutoff = 0.0001,
                  t_print_freq = 100, return_t = FALSE, r_print_freq = 100,
                  description = TRUE) {
  # w/out ci
  if (ci == FALSE) {
  acc_em_point(data, test, disease, covariate, mnar,
               show_boot = FALSE, description,
               show_t,
               t_max, cutoff,
               t_print_freq = t_print_freq, return_t = return_t,
               r_print_freq = r_print_freq)
  }
  # w ci
  else {
    acc_em_boot_ci(data, test, disease, covariate, mnar,
                   show_boot = TRUE, description,
                   show_t, t_max, cutoff,
                   ci_level = ci_level, ci_type = ci_type,
                   R = R, seednum = seednum,
                   t_print_freq = t_print_freq,
                   r_print_freq = r_print_freq)
  }
}
