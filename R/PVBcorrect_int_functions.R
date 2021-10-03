# Functions for internal use by other functions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Author: Wan Nor Arifin
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# For lengthy & complicated functions,
# to be called by PVBcorrect_functions.R
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
  if (show_boot == TRUE) {
    show_fit = FALSE  # force set as FALSE
    if (exists("counter")) {
      if(counter %% r_print_freq == 0) {cat("=== Boot Iteration =", counter, "===\n")}
      counter <<- counter + 1
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
  counter <<- 0
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
  # weight_k must be declared outside this function
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
  weight_k <<- rep(1, nrow(data_pseudo))  # to global environment
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
    # P(T=1|D=1)
    acc_values[1] = predict(em_out$model_b, data.frame(D=1), type="response")
    # P(T=0|D=0)
    acc_values[2] = 1 - predict(em_out$model_b, data.frame(D=0), type="response")
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
    counter <<- counter + 1  # <<- so that counter won't reset with each boot iter
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
                     r_print_freq = TRUE)
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
