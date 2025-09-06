# Functions for internal use by other functions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Author: Wan Nor Arifin
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# For lengthy & complicated functions,
# to be called by PVBcorrect_functions.R
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Basic internal functions ====
# use by some functions
snsp = function(data, test, disease) {
  snsp_values = c(0, 0)
  tbl = table(data[, test], data[, disease])
  snsp_values[1] = tbl[2,2]/sum(tbl[,2])  # Sn
  snsp_values[2] = tbl[1,1]/sum(tbl[,1])  # Sp
  return(snsp_values)
}

acc = function(data, test, disease) {
  acc_values = c(0, 0, 0, 0)
  tbl = table(data[, test], data[, disease])
  acc_values[1] = tbl[2,2]/sum(tbl[,2])  # Sn
  acc_values[2] = tbl[1,1]/sum(tbl[,1])  # Sp
  acc_values[3] = tbl[2,2]/sum(tbl[2,])  # PPV
  acc_values[4] = tbl[1,1]/sum(tbl[1,])  # PPV
  return(acc_values)
}

# EBG ====
# Begg & Greenes, 1983; Alonzo, 2005; He & McDermott, 2012

# EBG using GLM, point estimate, internal use
#' @importFrom stats reformulate glm predict
acc_ebg_point = function(data, test, disease, covariate = NULL, saturated_model = FALSE,
                         show_fit = FALSE, show_boot = FALSE, description = TRUE,
                         data_all = NULL,
                         r_print_freq = 100) {
  # setup empty vector
  acc_values = rep(NA, 4)

  # data = complete case only for modeling logistic model
  # glm will automatically drop na cases
  # if no covariate
  if (is.null(covariate)) {
    fmula = reformulate(test, disease)
    fit = glm(fmula, data = data, family = "binomial")
  }
  # if covariate names supplied as vector, not NULL
  else {
    if (saturated_model == FALSE) {
      fmula = reformulate(c(test, covariate), disease)
    } else {
      fmula = reformulate(paste(c(test, covariate), collapse = " * "), disease)
    }
    fit = glm(fmula, data = data, family = "binomial")
  }

  # data_all = full data, preds for each data points, used with boot
  if (is.null(data_all)) {
    data_all = data
  }

  # get measures from fitted model
  preds = predict(fit, data_all, type = "response")  # Predicted for complete & incomplete data
  # Sn, P(T=1|D=1)
  acc_values[1] = sum(data_all[, test] * preds) / sum(preds)
  # Sp, P(T=0|D=0)
  acc_values[2] = sum((1 - data_all[, test]) * (1 - preds)) / sum(1 - preds)
  # PPV, P(D=1|T=1)
  acc_values[3] = sum(preds * acc_values[1]) / sum(preds * acc_values[1] + (1 - preds) * (1 - acc_values[2]))
  # NPV, P(D=0|T=0)
  acc_values[4] = sum((1 - preds) * acc_values[2]) / sum((1 - preds) * acc_values[2] + preds * (1- acc_values[1]))

  # output
  if (show_boot == TRUE) {  # don't remove, this is for *_boot_ci use
    show_fit = FALSE  # force set as FALSE
    if (exists("counter")) {  # counter is set in *_boot_ci
      if(counter %% r_print_freq == 0) {cat("=== Boot Iteration =", counter, "===\n")}
      counter <<- counter + 1  # updated in .GlobalEnv, else R will clear it!
    }
  }
  if (show_fit == TRUE) {
    # special option with formatted, show model fit summary
    cat("\nModel Fit Summary:\n")
    print(summary(fit))
  }
  acc_ebg_out = data.frame(Est = acc_values,
                           row.names = c("Sn", "Sp", "PPV", "NPV"))
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Extended Begg and Greenes' Method\n\n")
  }
  acc_ebg_list = list(acc_results = acc_ebg_out)
  return(acc_ebg_list)
}

# Bootstrap fun for EBG, internal use
acc_ebg_boot_fun = function(data_verified, indices, test, disease, covariate = NULL, saturated_model = FALSE,
                            show_boot = FALSE, data_all, r_print_freq = 100) {
  data = data_verified[indices, ]
  out = acc_ebg_point(data = data, test, disease, covariate, saturated_model,
                      show_fit = FALSE, show_boot = show_boot, description = FALSE,
                      data_all = data_all, r_print_freq = r_print_freq)
  out = as.matrix(out$acc_results)
  return(out)
}

# EBG using GLM, bootstrapped CI, internal use
#' @importFrom stats sd
acc_ebg_boot_ci = function(data, test, disease, covariate = NULL, saturated_model = FALSE,
                           show_fit = FALSE, show_boot = FALSE, description = TRUE,
                           ci_level = .95, ci_type = "basic",
                           R = 999, seednum = NULL,
                           r_print_freq = 100) {
  # check ci_type, based on 'boot' options
  ci_type_allowed = c("norm", "basic", "perc", "bca")
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

  # split data by verified & unverified
  data_verified = data[!is.na(data[, disease]), ]  # select disease != NA
  data_all = data

  # check for presence of seednum
  if (!is.null(seednum)) {
    set.seed(seednum)
  }  # else boot will select its own seednum

  # run ebg & get ci by bootstrap
  counter <<- 0  # set in .GlobalEnv, else R will clear it!
  acc_ebg_boot_data = boot::boot(data = data_verified, statistic = acc_ebg_boot_fun, R = R,
                           test = test, disease = disease, covariate = covariate, saturated_model = saturated_model,
                           show_boot = show_boot, data_all = data_all, r_print_freq = r_print_freq)
  if(show_boot == TRUE) {cat("[ Total Boot Iteration =", R, "]\n\n")}
  rm(counter, envir = .GlobalEnv)
  acc_ebg_boot_est = acc_ebg_boot_data$t0
  acc_ebg_boot_se = apply(acc_ebg_boot_data$t, 2, sd)
  acc_ebg_boot_ci_data = lapply(1:4, function(i) boot::boot.ci(acc_ebg_boot_data, conf = ci_level, type = ci_type, index = i))
  acc_ebg_boot_ci = lapply(acc_ebg_boot_ci_data, function(list) list[[4]][(length(list[[4]])-1):length(list[[4]])])
  acc_ebg_boot_out = t(sapply(1:4, function(i) c(acc_ebg_boot_est[i, ], acc_ebg_boot_se[i], acc_ebg_boot_ci[[i]])))
  dimnames(acc_ebg_boot_out) = list(c("Sn", "Sp", "PPV", "NPV"), c("Est", "SE", "LowCI", "UppCI"))

  # output
  if (show_fit == TRUE) {
    acc_ebg_point(data, test, disease, covariate, saturated_model, show_fit = TRUE, description = FALSE)
  }
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Extended Begg and Greenes' Method\n\n")
  }
  #print(acc_ebg_boot_out)
  acc_ebg_boot_list = list(boot_data = acc_ebg_boot_data,
                           boot_ci_data = acc_ebg_boot_ci_data,
                           acc_results = acc_ebg_boot_out)  # allows inspection of boot data
  return(acc_ebg_boot_list)
}

# IPW ====
# Alonzo & Pepe, 2005

# IPW internal functions
acc_ipw_point = function(data, test, disease, covariate = NULL, saturated_model = FALSE,
                         show_fit = FALSE, show_boot = FALSE, description = TRUE,
                         data_all = NULL,
                         r_print_freq = 100) {
  # data = original data
  # add verification status
  data$verified = 1  # verified: 1 = yes, 0 = no
  data[is.na(data[, disease]), "verified"] = 0

  # setup empty vector
  acc_values = rep(NA, 4)

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

  # data_all = full data, preds for each data points, used with boot
  if (is.null(data_all)) {
    data_all = data
  }

  # recode NA to -1
  data_ps = data
  data_ps[is.na(data_ps[, disease]), disease] = -1  # so as data$verified * data$disease == 0 for NA

  # get measures from fitted model
  # Sn, P(T=1|D=1)
  acc_values[1] = sum((data_ps[, test] * data_ps[, "verified"] * data_ps[, disease]) / data_ps$ps) /
    sum((data_ps[, "verified"] * data_ps[, disease]) / data_ps$ps)
  # Sp, P(T=0|D=0)
  acc_values[2] = sum(((1 - data_ps[, test]) * data_ps[, "verified"] * (1 - data_ps[, disease])) / data_ps$ps) /
    sum((data_ps[, "verified"] * (1 - data_ps[, disease])) / data_ps$ps)
  # For PPV & NPV, follow formulas from Arifin & Yusoff (2022), Stats in Med
  # Need to get P(D=1), preds
  # Rely on EBG for this
  # from EBG ---
  # if no covariate
  if (is.null(covariate)) {
    fmula = reformulate(test, disease)
    fit_ = glm(fmula, data = data, family = "binomial")  # for preds
  }
  # if covariate names supplied as vector, not NULL
  else {
    if (saturated_model == FALSE) {
      fmula = reformulate(c(test, covariate), disease)
    } else {
      fmula = reformulate(paste(c(test, covariate), collapse = " * "), disease)
    }
    fit_ = glm(fmula, data = data, family = "binomial")  # for preds
  }
  # data_all = full data, preds for each data points, used with boot
  if (is.null(data_all)) {
    data_all = data
  }
  preds = predict(fit_, data_all, type = "response")  # Predicted for complete & incomplete data
  # end from EBG
  # Use preds from EBG, and Sn & Sp from IPW
  # PPV, P(D=1|T=1)
  acc_values[3] = sum(preds * acc_values[1]) / sum(preds * acc_values[1] + (1 - preds) * (1 - acc_values[2]))
  # NPV, P(D=0|T=0)
  acc_values[4] = sum((1 - preds) * acc_values[2]) / sum((1 - preds) * acc_values[2] + preds * (1- acc_values[1]))


  # output
  if (show_boot == TRUE) {  # don't remove, this is for *_boot_ci use
    show_fit = FALSE  # force set as FALSE
    if (exists("counter")) {  # counter is set in *_boot_ci
      if(counter %% r_print_freq == 0) {cat("=== Boot Iteration =", counter, "===\n")}
      counter <<- counter + 1  # updated in .GlobalEnv, else R will clear it!
    }
  }
  if (show_fit == TRUE) {
    # special option with formatted, show model fit summary
    cat("\nModel Fit Summary:\n")
    print(summary(fit))
  }
  acc_ipw_out = data.frame(Est = acc_values,
                           row.names = c("Sn", "Sp", "PPV", "NPV"))
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Inverse Probability Weighting Estimator Method\n\n")
  }
  acc_ipw_list = list(acc_results = acc_ipw_out)
  return(acc_ipw_list)
}

acc_ipw_boot_fun = function(data, indices, test, disease, covariate = NULL, saturated_model = FALSE,
                            show_boot = FALSE, data_all, r_print_freq = 100) {
  data = data[indices, ] # resample from all (verified + unverified)
  # bcs ipw involves T & V
  out = acc_ipw_point(data = data, test, disease, covariate, saturated_model,
                      show_fit = FALSE, show_boot = show_boot, description = FALSE,
                      data_all = data_all, r_print_freq = r_print_freq)
  out = as.matrix(out$acc_results)
  return(out)
}

acc_ipw_boot_ci = function(data, test, disease, covariate = NULL, saturated_model = FALSE,
                           show_fit = FALSE, show_boot = FALSE, description = TRUE,
                           ci_level = .95, ci_type = "basic",
                           R = 999, seednum = NULL,
                           r_print_freq = 100) {
  # check ci_type, based on 'boot' options
  ci_type_allowed = c("norm", "basic", "perc", "bca")
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

  # split data by verified & unverified
  # data_verified = data[!is.na(data[, disease]), ]  # select disease != NA
  data_all = data

  # check for presence of seednum
  if (!is.null(seednum)) {
    set.seed(seednum)
  }  # else boot will select its own seednum

  # run ipw & get ci by bootstrap
  counter <<- 0  # set in .GlobalEnv, else R will clear it!
  acc_ipw_boot_data = boot::boot(data = data, statistic = acc_ipw_boot_fun, R = R,
                                 test = test, disease = disease, covariate = covariate, saturated_model = saturated_model,
                                 show_boot = show_boot, data_all = data_all, r_print_freq = r_print_freq)
  if(show_boot == TRUE) {cat("[ Total Boot Iteration =", R, "]\n\n")}
  rm(counter, envir = .GlobalEnv)
  acc_ipw_boot_est = acc_ipw_boot_data$t0
  acc_ipw_boot_se = apply(acc_ipw_boot_data$t, 2, sd)
  acc_ipw_boot_ci_data = lapply(1:4, function(i) boot::boot.ci(acc_ipw_boot_data, conf = ci_level, type = ci_type, index = i))
  acc_ipw_boot_ci = lapply(acc_ipw_boot_ci_data, function(list) list[[4]][(length(list[[4]])-1):length(list[[4]])])
  acc_ipw_boot_out = t(sapply(1:4, function(i) c(acc_ipw_boot_est[i, ], acc_ipw_boot_se[i], acc_ipw_boot_ci[[i]])))
  dimnames(acc_ipw_boot_out) = list(c("Sn", "Sp", "PPV", "NPV"), c("Est", "SE", "LowCI", "UppCI"))

  # output
  if (show_fit == TRUE) {
    acc_ipw_point(data, test, disease, verified, covariate, saturated_model, show_fit = TRUE, description = FALSE)
  }
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Inverse Probability Weighting Estimator Method\n\n")
  }
  acc_ipw_boot_list = list(boot_data = acc_ipw_boot_data,
                           boot_ci_data = acc_ipw_boot_ci_data,
                           acc_results = acc_ipw_boot_out)  # allows inspection of boot data
  return(acc_ipw_boot_list)
}

# SIPW ====
acc_sipw_point = function(data, test, disease, covariate = NULL, option = 2, saturated_model = FALSE,
                          b = 1000, seednum = NULL, return_data = FALSE, return_detail = FALSE,
                          show_boot = FALSE, description = TRUE,
                          r_print_freq = 100) {
  # data = original data
  N_data = nrow(data)  # original sample size
  # add verification status
  data$verified = 1  # verified: 1 = yes, 0 = no
  data[is.na(data[, disease]), "verified"] = 0

  # setup empty vector
  acc_values = rep(NA, 4)

  # gen ps
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    fmula = reformulate(test, "verified")
  }
  # if covariate names supplied as vector, not NULL
  else {
    if (saturated_model == FALSE) {
      fmula = reformulate(c(test, covariate), "verified")
    } else {
      fmula = reformulate(paste(c(test, covariate), collapse = " * "), "verified")
    }
  }
  fit = glm(fmula, data = data, family = "binomial")
  # save ps to original data
  data$ps = predict(fit, type = "response")
  # transform ps to weight
  if (option == 1) {
    # IPW weight
    data$ps_w = ifelse(data$verified == 1, 1/data$ps, 1/(1-data$ps))
  }
  if (option == 2) {
    # W_h weight, Krautenbacher 2017
    data$ps_w = ifelse(data$verified == 1, max(unique(data$ps)) / data$ps,
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
    sampled = sample(1:n_data, N_data, TRUE, prob = data1$p_ipb)  # bootstrap w weight
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

  # output
  if (show_boot == TRUE) {  # don't remove, this is for *_boot_ci use
    if (exists("counter")) {  # counter is set in *_boot_ci
      if(counter %% r_print_freq == 0) {cat("=== Boot Iteration =", counter, "===\n")}
      counter <<- counter + 1  # updated in .GlobalEnv, else R will clear it!
    }
  }

  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Scaled Inverse Probability Weighted Resampling Method\n\n")
  }

  # return data and results
  if (return_data == TRUE) {
    return(list(data = data_resample_list, acc = acc_values))
  }

  # return only results
  acc_sipw_out = data.frame(Est = acc_values,
                            row.names = c("Sn", "Sp", "PPV", "NPV"))
  acc_sipw_list = list(acc_results = acc_sipw_out)
  return(acc_sipw_list)
}



acc_sipw_boot_fun = function(data, indices, test, disease, covariate = NULL,
                             saturated_model = FALSE, option = 2, b = 1000,
                             seednum = NULL, # seednum important for sipw! replicable result
                             show_boot = FALSE, r_print_freq = 100) {
  data = data[indices, ] # resample from all (verified + unverified)
  # bcs ipw involves T & V, not D with NA
  out = acc_sipw_point(data = data, test, disease, covariate,
                       saturated_model = saturated_model, option = option, b = b, seednum = seednum,
                       description = FALSE, return_data = FALSE,
                       show_boot = show_boot, r_print_freq = r_print_freq)
  out = as.matrix(out$acc_results)
  return(out)
}

acc_sipw_boot_ci = function(data, test, disease, covariate = NULL,
                            saturated_model = FALSE, option = 2,
                            show_boot = FALSE, description = TRUE,
                            ci_level = .95, ci_type = "basic",
                            R = 999, b = 1000, seednum = NULL,
                            r_print_freq = 100) {
  # check ci_type, based on 'boot' options
  ci_type_allowed = c("norm", "basic", "perc", "bca")
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

  # check for presence of seednum
  if (!is.null(seednum)) {
    set.seed(seednum)
  }  # else boot will select its own seednum

  # run sipw & get ci by bootstrap
  counter <<- 0
  acc_sipw_boot_data = boot::boot(data = data, statistic = acc_sipw_boot_fun, R = R,
                                  test = test, disease = disease, covariate = covariate,
                                  saturated_model = saturated_model, option = option, b = b, seednum = seednum,
                                  show_boot = show_boot,
                                  r_print_freq = r_print_freq)
  if(show_boot == TRUE) {cat("[ Total Boot Iteration =", R, "]\n\n")}
  rm(counter, envir = .GlobalEnv)
  acc_sipw_boot_est = acc_sipw_boot_data$t0
  acc_sipw_boot_se = apply(acc_sipw_boot_data$t, 2, sd)
  acc_sipw_boot_ci_data = lapply(1:4, function(i) boot::boot.ci(acc_sipw_boot_data, conf = ci_level, type = ci_type, index = i))
  acc_sipw_boot_ci = lapply(acc_sipw_boot_ci_data, function(list) list[[4]][(length(list[[4]])-1):length(list[[4]])])
  acc_sipw_boot_out = t(sapply(1:4, function(i) c(acc_sipw_boot_est[i, ], acc_sipw_boot_se[i], acc_sipw_boot_ci[[i]])))
  dimnames(acc_sipw_boot_out) = list(c("Sn", "Sp", "PPV", "NPV"), c("Est", "SE", "LowCI", "UppCI"))

  # output
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Scaled Inverse Probability Weighted Resampling Method\n\n")
  }
  #print(acc_ebg_boot_out)
  acc_sipw_boot_list = list(boot_data = acc_sipw_boot_data,
                            boot_ci_data = acc_sipw_boot_ci_data,
                            acc_results = acc_sipw_boot_out)  # allows inspection of boot data
  return(acc_sipw_boot_list)
}

# SIPW-B ====
acc_sipwb_point = function(data, test, disease, covariate = NULL, option = 2, saturated_model = FALSE,
                           rel_size = 1, # ratio control:case, D=0:D=1
                           b = 1000, seednum = NULL, return_data = FALSE, return_detail = FALSE,
                           show_boot = FALSE, description = TRUE,
                           r_print_freq = 100) {
  # data = original data
  # rel_size = relative size of d0 to d1 in case-ctrl, default 1 i.e. 1:1
  N_data = nrow(data)  # original sample size
  # add verification status
  data$verified = 1  # verified: 1 = yes, 0 = no
  data[is.na(data[, disease]), "verified"] = 0

  # setup empty vector
  acc_values = rep(NA, 2)  # Sn Sp only

  # gen ps
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    fmula = reformulate(test, "verified")
  }
  # if covariate names supplied as vector, not NULL
  else {
    if (saturated_model == FALSE) {
      fmula = reformulate(c(test, covariate), "verified")
    } else {
      fmula = reformulate(paste(c(test, covariate), collapse = " * "), "verified")
    }
  }
  fit = glm(fmula, data = data, family = "binomial")
  # save ps to original data
  data$ps = predict(fit, type = "response")
  # transform ps to weight
  if (option == 1) {
    # IPW weight
    data$ps_w = ifelse(data$verified == 1, 1/data$ps, 1/(1-data$ps))
  }
  if (option == 2) {
    # W_h weight, Krautenbacher 2017
    data$ps_w = ifelse(data$verified == 1, max(unique(data$ps)) / data$ps,
                       max(unique((1-data$ps))) / (1-data$ps))
  }

  # data1 = complete case
  data1 = na.omit(data)
  n_data = nrow(data1)
  n_d = table(data1[, disease])
  multip = n_d[1] / n_d[2]  # d0:d1 ratio to inflate d1
  data1$ps_w_d1 = data1$ps_w
  data1$ps_w_d1[data1[, disease] == 1] = data1$ps_w_d1[data1[, disease] == 1] * multip
  k = sum(data1$ps_w_d1[data1[, disease] == 0])/  # correction factor
    sum(data1$ps_w_d1[data1[, disease] == 1])     # to allow almost equal d1:d0 ratio
  data1$ps_w_d1[data1[, disease] == 1] = data1$ps_w_d1[data1[, disease] == 1] * (k / rel_size)
  data1$p_ipb_d1 = data1$ps_w_d1 / sum(data1$ps_w_d1)

  # resample
  data_resample_list = NULL
  set.seed(seednum)
  i = 1
  while (i < b + 1) {
    # resample
    sampled = sample(1:n_data, N_data, TRUE, prob = data1$p_ipb_d1)  # bootstrap w weight
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

  snsp_resample_list = t(sapply(data_resample_list, snsp,
                                test = test, disease = disease))

  acc_values = apply(snsp_resample_list, 2, mean)

  # output
  if (show_boot == TRUE) {  # don't remove, this is for *_boot_ci use
    if (exists("counter")) {  # counter is set in *_boot_ci
      if(counter %% r_print_freq == 0) {cat("=== Boot Iteration =", counter, "===\n")}
      counter <<- counter + 1  # updated in .GlobalEnv, else R will clear it!
    }
  }

  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Scaled Inverse Probability Weighted Balanced Resampling Method\n\n")
  }

  # return data and results
  if (return_data == TRUE) {
    return(list(data = data_resample_list, acc = acc_values))
  }

  # return only results
  acc_sipwb_out = data.frame(Est = acc_values,
                             row.names = c("Sn", "Sp"))
  acc_sipwb_list = list(acc_results = acc_sipwb_out)
  return(acc_sipwb_list)
}

acc_sipwb_boot_fun = function(data, indices, test, disease, covariate = NULL,
                              saturated_model = FALSE, option = 2, b = 1000,
                              seednum = NULL, # seednum important for sipw! replicable result
                              rel_size = 1,
                              show_boot = FALSE, r_print_freq = 100) {
  data = data[indices, ] # resample from all (verified + unverified)
  # bcs ipw involves T & V, not D with NA
  out = acc_sipwb_point(data = data, test, disease, covariate,
                        saturated_model = saturated_model, option = option, b = b, seednum = seednum,
                        rel_size = rel_size,
                        description = FALSE, return_data = FALSE,
                        show_boot = show_boot, r_print_freq = r_print_freq)
  out = as.matrix(out$acc_results)
  return(out)
}

acc_sipwb_boot_ci = function(data, test, disease, covariate = NULL,
                             saturated_model = FALSE, option = 2,
                             rel_size = 1,
                             show_boot = FALSE, description = TRUE,
                             ci_level = .95, ci_type = "basic",
                             R = 999, b = 1000, seednum = NULL,
                             r_print_freq = 100) {
  # check ci_type, based on 'boot' options
  ci_type_allowed = c("norm", "basic", "perc", "bca")
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

  # check for presence of seednum
  if (!is.null(seednum)) {
    set.seed(seednum)
  }  # else boot will select its own seednum

  # run sipwb & get ci by bootstrap
  counter <<- 0
  acc_sipwb_boot_data = boot::boot(data = data, statistic = acc_sipwb_boot_fun, R = R,
                                   test = test, disease = disease, covariate = covariate,
                                   saturated_model = saturated_model, option = option, rel_size = rel_size,
                                   b = b, seednum = seednum,
                                   show_boot = show_boot,
                                   r_print_freq = r_print_freq)
  if(show_boot == TRUE) {cat("[ Total Boot Iteration =", R, "]\n\n")}
  rm(counter, envir = .GlobalEnv)
  acc_sipwb_boot_est = acc_sipwb_boot_data$t0
  acc_sipwb_boot_se = apply(acc_sipwb_boot_data$t, 2, sd)
  acc_sipwb_boot_ci_data = lapply(1:2, function(i) boot::boot.ci(acc_sipwb_boot_data, conf = ci_level, type = ci_type, index = i))
  acc_sipwb_boot_ci = lapply(acc_sipwb_boot_ci_data, function(list) list[[4]][(length(list[[4]])-1):length(list[[4]])])
  acc_sipwb_boot_out = t(sapply(1:2, function(i) c(acc_sipwb_boot_est[i, ], acc_sipwb_boot_se[i], acc_sipwb_boot_ci[[i]])))
  dimnames(acc_sipwb_boot_out) = list(c("Sn", "Sp"), c("Est", "SE", "LowCI", "UppCI"))

  # output
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: Scaled Inverse Probability Weighted Balanced Resampling Method\n\n")
  }
  #print(acc_ebg_boot_out)
  acc_sipwb_boot_list = list(boot_data = acc_sipwb_boot_data,
                             boot_ci_data = acc_sipwb_boot_ci_data,
                             acc_results = acc_sipwb_boot_out)  # allows inspection of boot data
  return(acc_sipwb_boot_list)
}

# EM ====
# Kosinski & Barnhart, 2003

# EM Algorithm Function
# Components: a = disease, b = diagnostic, c = missing data mechanism
#' @importFrom stats glm predict
em_fun = function(data_pseudo,
                  show_t = TRUE,
                  t_max = 500, cutoff = 0.0001,
                  a, b, c,
                  index_1, index_2, index_3,
                  t_print_freq = 100) {
  # Initialize values
  coef = vector("list", 0)  # empty vector list of coef

  # EM Iteration
  t = 1
  # weight_k must be declared outside this function, i.e. in .GlobalEnv
  while (t < t_max + 1) {
    # a -- P(D|X)
    model_a = glm(a, data = data_pseudo, family = "binomial", weights = weight_k)
    su_model_a = summary(model_a)
    coef_a = su_model_a$coefficients
    fitted_pa = model_a$fitted.values
    # b -- P(T|D,X)
    model_b = glm(b, data = data_pseudo, family = "binomial", weights = weight_k)
    su_model_b = summary(model_b)
    coef_b = su_model_b$coefficients
    fitted_pb = model_b$fitted.values
    # c -- P(V|T,X,D)
    model_c = glm(c, data = data_pseudo, family = "binomial", weights = weight_k)
    su_model_c = summary(model_c)
    coef_c = su_model_c$coefficients
    fitted_pc = model_c$fitted.values
    # all
    coef = append(coef, list(list(coef_a, coef_b, coef_c)))
    # weight
    fitted_ps = list(fitted_pa, fitted_pb, fitted_pc)
    ys = list(model_a$y, model_b$y, model_c$y)
    p0 = (fitted_ps[[1]]^ys[[1]]) * ((1 - fitted_ps[[1]])^(1 - ys[[1]]))  # P(D|X)
    p1 = (fitted_ps[[2]]^ys[[2]]) * ((1 - fitted_ps[[2]])^(1 - ys[[2]]))  # P(T|D,X)
    p2 = (fitted_ps[[3]]^ys[[3]]) * ((1 - fitted_ps[[3]])^(1 - ys[[3]]))  # P(V|T,X,D)
    pk = p0 * p1 * p2  # P(V,T,D)
    weight_k[index_2] <<- pk[index_2] / (pk[index_2] + pk[index_3])  # =/<- not working here
    weight_k[index_3] <<- 1 - weight_k[index_2]  # =/<- not working here
    # check if change in coef < 0.001 in abs value
    if (t < 2) {
      diffs_t = cutoff + 1
    } else {  # diffs in coef
      nrow_a = nrow(coef[[t]][[1]])
      nrow_b = nrow(coef[[t]][[2]])
      nrow_c = nrow(coef[[t]][[3]])
      diffs_t_a = coef[[t]][[1]][1:nrow_a] - coef[[t-1]][[1]][1:nrow_a]  # comp a
      diffs_t_b = coef[[t]][[2]][1:nrow_b] - coef[[t-1]][[2]][1:nrow_b]  # comp b
      diffs_t_c = coef[[t]][[3]][1:nrow_c] - coef[[t-1]][[3]][1:nrow_c]  # comp c
      diffs_t = c(diffs_t_a, diffs_t_b, diffs_t_c)
    }
    # early stopping rule
    if (all(abs(diffs_t) < cutoff)) {
      cat(paste("[ EM converged at t =", t, "when all changes <", cutoff, "]\n\n"))
      break
    }
    # if early stopping rule is not met, run until t max
    else {
      if(show_t == TRUE) {
        # printing frequency
        if(t %% t_print_freq == 0) {cat(paste("=== Current EM iteration t =", t, " ===\n"))}
      }
      t = t + 1
    }
  }

  # Output
  out = list(model_a = model_a, model_b = model_b, model_c = model_c,
             diffs_t = diffs_t, t = t, coef = coef)
  # this output model details and coef
  return(out)
}

# EM-based method, point estimate, internal use
#' @importFrom stats reformulate predict
acc_em_point = function(data, test, disease, covariate = NULL, mnar = TRUE,
                        show_t = TRUE, # it's better show_t = TRUE bcs EM takes long time to finish
                        show_boot = FALSE, description = TRUE,
                        t_max = 500, cutoff = 0.0001,
                        t_print_freq = 100, return_t = FALSE,
                        return_em_out = FALSE, # for debug & detailed check
                        r_print_freq = 100) {
  # setup data for EM
  # add verification status
  data$verified = 1  # verified: 1 = yes, 0 = no
  data[is.na(data[, disease]), "verified"] = 0
  data_verified = data[data$verified == 1, ]
  data_unverified = data[data$verified == 0, ]  # 1U, unverified U observation
  data_unverified_stacked = rbind(data_unverified, data[data$verified == 0, ])  # 2U, replicate U observation
  # create 0, 1 for 2U rows
  data_unverified_stacked[1:nrow(data[data$verified == 0, ]), disease] = 0        # 1st U rows
  data_unverified_stacked[(nrow(data[data$verified == 0, ])+1):nrow(data_unverified_stacked), disease] = 1  # 2nd U rows
  # pseudo data
  data_pseudo = rbind(data_verified, data_unverified_stacked)
  # index k
  index_1 = which(data_pseudo$verified == 1)  # verified
  index_2 = which(data_pseudo$verified == 0 & data_pseudo[, disease] == 0)  # unverified U
  index_3 = which(data_pseudo$verified == 0 & data_pseudo[, disease] == 1)  # unverified 2U

  # setup formula
  # if no covariate
  if (is.null(covariate)) {
    a = reformulate("1", disease)
    b = reformulate(disease, test)
    c = reformulate(c(test, disease), "verified")
    if (mnar == FALSE) {
      c = reformulate(c(test), "verified")
    }
  }
  # if covariate names supplied as vector, not NULL
  else {
    a = reformulate(covariate, disease)
    b = reformulate(c(disease, covariate), test)
    c = reformulate(c(test, covariate, disease), "verified")
    if (mnar == FALSE) {
      c = reformulate(c(test, covariate), "verified")
    }
  }

  # the EM run
  weight_k <<- rep(1, nrow(data_pseudo))  # to .GlobalEnv, else R will clear it!
  em_out = em_fun(data_pseudo,
                  show_t = show_t,
                  t_max = t_max, cutoff = cutoff,
                  a = a, b = b, c = c,
                  index_1, index_2, index_3,
                  t_print_freq = t_print_freq)
  rm(weight_k, envir = .GlobalEnv)  # rm from environment

  # calculate measures
  # setup empty vector
  acc_values = rep(NA, 4)
  # basic, no covariate, faster
  if (is.null(covariate)) {
    data_1 = data.frame(1); colnames(data_1) = disease
    data_0 = data.frame(0); colnames(data_0) = disease
    # P(T=1|D=1)
    acc_values[1] = predict(em_out$model_b, data_1, type="response")
    # P(T=0|D=0)
    acc_values[2] = 1 - predict(em_out$model_b, data_0, type="response")
    # PPV, P(D=1|T=1)
    acc_values[3] = (predict(em_out$model_a, data.frame(1), type="response") * acc_values[1]) /
      (predict(em_out$model_a, data.frame(1), type="response") * acc_values[1] + (1 - predict(em_out$model_a, data.frame(1), type="response")) * (1 - acc_values[2]))
    # NPV, P(D=0|T=0)
    acc_values[4] = ((1 - predict(em_out$model_a, data.frame(1), type="response")) * acc_values[2]) /
      ((1 - predict(em_out$model_a, data.frame(1), type="response")) * acc_values[2] + predict(em_out$model_a, data.frame(1), type="response") * (1- acc_values[1]))
  }
  # marginal estimates, w covariate
  else {
    # prepare data for marginal estimate
    data_1 = data_0 = data_pseudo[c(index_1, index_2),]
    data_1[, disease] = 1; data_0[, disease] = 0
    # Sn, P(T=1|D=1,X)
    acc_values[1] = sum(predict(em_out$model_b, data_1, type="response") * predict(em_out$model_a, data_1, type="response")) /
      sum(predict(em_out$model_a, data_1, type="response"))
    # Sp, P(T=0|D=0,X)
    acc_values[2] = sum((1 - predict(em_out$model_b, data_0, type="response")) * (1 - predict(em_out$model_a, data_0, type="response"))) /
      sum(1 - predict(em_out$model_a, data_0, type="response"))
    # PPV, P(D=1|T=1,X)
    acc_values[3] = sum(predict(em_out$model_a, data_1, type="response") * acc_values[1]) /
      sum(predict(em_out$model_a, data_1, type="response") * acc_values[1] + (1 - predict(em_out$model_a, data_1, type="response")) * (1 - acc_values[2]))
    # actually it's ok to use either data_1 / data_0 for model_a bcs it does not involve D as predictor
    # NPV, P(D=0|T=0,X)
    acc_values[4] = sum((1 - predict(em_out$model_a, data_0, type="response")) * acc_values[2]) /
      sum((1 - predict(em_out$model_a, data_0, type="response")) * acc_values[2] + predict(em_out$model_a, data_0, type="response") * (1- acc_values[1]))
  }

  if (show_boot == TRUE) {
    cat(paste("Finished Boot Iteration =", counter, "\n"))
    cat("=========================\n")
    cat("\n")
    counter <<- counter + 1  # <<- updated in .GlobalEnv,
    # so that counter won't reset with each boot iteration
  }

  # normal output
  if (return_t == FALSE) {
    acc_em_out = list(acc_results = data.frame(Est = acc_values,
                            row.names = c("Sn", "Sp", "PPV", "NPV")))
    # acc_em_out = data.frame(Est = acc_values,
    #                         row.names = c("Sn", "Sp"))
  }
  # if t is requested for check of convergence in multiple run
  else {
    acc_em_out = list(acc_results = data.frame(Est = acc_values,
                                       row.names = c("Sn", "Sp", "PPV", "NPV")),
                      t = em_out$t)
    # acc_em_out = list(Est = data.frame(Est = acc_values,
    #                                    row.names = c("Sn", "Sp")),
    #                   t = em_out$t)
  }
  # if detailed em_out is requested
  if (return_em_out == TRUE) {
    acc_em_out = list(Detail = em_out, Result = acc_em_out)
  }
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB (MNAR): EM-based Method\n\n")
  }
  return(acc_em_out)
}

# Bootstrap fun for EM-based method, internal use
acc_em_boot_fun = function(data_verified, indices, test, disease, covariate = NULL, mnar = TRUE,
                           show_t = TRUE, t_max = 500, cutoff = 0.0001,
                           t_print_freq = 100, data_unverified, r_print_freq = 100) {
  data_verified = data_verified[indices, ]
  data_all = rbind(data_verified, data_unverified)
  out = acc_em_point(data = data_all, test, disease, covariate, mnar = mnar,
                     show_t = show_t, show_boot = TRUE, description = FALSE,
                     t_max = t_max, cutoff = cutoff,
                     t_print_freq = t_print_freq, return_t = TRUE,  # need to return t to check convergence
                     return_em_out = FALSE,
                     r_print_freq = r_print_freq)
  out = as.matrix(rbind(out$acc_results, t = out$t))
  return(out)
}

# EM-based method, bootstrapped ci, internal use
#' @importFrom stats sd
acc_em_boot_ci = function(data, test, disease, covariate = NULL, mnar = TRUE,
                          show_boot = TRUE, description = TRUE,
                          show_t = TRUE, t_max = 500, cutoff = 0.0001,
                          ci_level = .95, ci_type = "basic",
                          R = 999, seednum = NULL,
                          t_print_freq = 10, r_print_freq = 100) {
  # check ci_type, based on 'boot' options
  ci_type_allowed = c("norm", "basic", "perc", "bca")
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

  # split data by verified & unverified
  data_verified = data[!is.na(data[, disease]), ]  # select disease != NA
  data_unverified = data[is.na(data[, disease]), ]  # select disease == NA

  # check for presence of seednum
  if (!is.null(seednum)) {
    set.seed(seednum)
  }  # else boot will select its own seednum

  # run em & get ci by bootstrap
  counter <<- 0  # starts counter for number of bootstrap iteration
  # set in .GlobalEnv, else R will clear it!
  acc_em_boot_data = boot::boot(data = data_verified, statistic = acc_em_boot_fun, R = R,
                          test = test, disease = disease, covariate = covariate, mnar = mnar,
                          show_t = show_t, t_max = t_max, cutoff = cutoff,
                          data_unverified = data_unverified,
                          t_print_freq = t_print_freq, r_print_freq = r_print_freq)
  if(show_boot == TRUE) {cat("[ Total Boot Iteration =", R, "]\n\n")}
  rm(counter, envir = .GlobalEnv)
  acc_em_boot_est = acc_em_boot_data$t0
  acc_em_boot_se = apply(acc_em_boot_data$t[, 1:4], 2, sd)
  acc_em_boot_ci_data = lapply(1:4, function(i) boot::boot.ci(acc_em_boot_data, conf = ci_level, type = ci_type, index = i))
  # in each list in ci_data, [[4]] refers to ci results, with last two values LL & UL
  acc_em_boot_ci = lapply(acc_em_boot_ci_data, function(list) list[[4]][(length(list[[4]])-1):length(list[[4]])])
  acc_em_boot_out = t(sapply(1:4, function(i) c(acc_em_boot_est[i, ], acc_em_boot_se[i], acc_em_boot_ci[[i]])))
  dimnames(acc_em_boot_out) = list(c("Sn", "Sp", "PPV", "NPV"), c("Est", "SE", "LowCI", "UppCI"))

  # output
  if (description == TRUE) {
    cat("Estimates of accuracy measures\nCorrected for PVB: EM-based Method\n\n")
  }
  #print(acc_em_boot_out)
  acc_em_boot_list = list(boot_data = acc_em_boot_data,
                          boot_ci_data = acc_em_boot_ci_data,
                          acc_results = acc_em_boot_out)  # allows inspection of boot data
  return(acc_em_boot_list)
}
