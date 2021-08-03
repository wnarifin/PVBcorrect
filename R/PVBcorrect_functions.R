# Functions for PVB correction
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Author: Wan Nor Arifin
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Functions ====

# Setup data
data_setup = function(test, disease, covariate = NULL) {
  # to include data checkup, i.e. T/D == 0/1 or not for ease of calculation & interpretation
  # covariate can include categorical / continuous variables
  data = data.frame(test = test, disease = disease )
  if (!is.null(covariate)) {
    data = cbind(data, data.frame(covariate))
  }
  return(data)
}

# View Table
view_table = function(data, show_unverified = FALSE) {
  test = data$test
  disease = data$disease
  if (show_unverified == TRUE) {
    tbl = table(Test = 2 - test, Disease = 2 - disease, useNA = "ifany")  # rearrange for epid view
    rownames(tbl) = c("yes","no")
    colnames(tbl) = c("yes","no","unverified")
  } else if (show_unverified == FALSE) {
    tbl = table(Test = 2 - test, Disease = 2 - disease)  # rearrange for epid view
    rownames(tbl) = colnames(tbl) = c("yes","no")
  }
  return(tbl)
}

# Sn & Sp
snsp = function(data, formatted = FALSE) {
  snsp_values = c(0, 0)
  tbl = table(data$test, data$disease)
  snsp_values[1] = tbl[2, 2] / sum(tbl[, 2])  # Sn
  snsp_values[2] = tbl[1, 1] / sum(tbl[, 1])  # Sp
  if (formatted == TRUE) {
    cat("Sensitivity = ", snsp_values[1], "\n")
    cat("Specificity = ", snsp_values[2])
  } else {
    return(snsp_values)
  }
}

# PPV & NPV, predictive values
pvs = function(data, formatted = FALSE) {
  pvs_values = c(0, 0)
  tbl = table(data$test, data$disease)
  pvs_values[1] = tbl[2,2]/sum(tbl[2,])  # PPV
  pvs_values[2] = tbl[1,1]/sum(tbl[1,])  # NPV
  if (formatted == TRUE) {
    cat("Positive Predictive Value = ", pvs_values[1], "\n")
    cat("Negative Predictive Value = ", pvs_values[2])
  } else {
    return(pvs_values)
  }
}

accuracy_measures = function(data, formatted = FALSE) {
  acc_values = c(0, 0)
  tbl = table(data$test, data$disease)
  acc_values[1] = tbl[2, 2] / sum(tbl[, 2])  # Sn
  acc_values[2] = tbl[1, 1] / sum(tbl[, 1])  # Sp
  acc_values[3] = tbl[2, 2] / sum(tbl[2, ])  # PPV
  acc_values[4] = tbl[1, 1] / sum(tbl[1, ])  # NPV
  if (formatted == TRUE) {
    cat("Sensitivity = ", acc_values[1], "\n")
    cat("Specificity = ", acc_values[2], "\n")
    cat("Positive Predictive Value = ", acc_values[3], "\n")
    cat("Negative Predictive Value = ", acc_values[4])
  } else {
    return(acc_values)
  }
}

# EBG using GLM
snsp_ebg = function(data, covariate = FALSE, formatted = FALSE) {
  snsp_values = c(0, 0)
  # data = complete case only for modeling, automatic by glm
  if (covariate == FALSE) {
    fit = glm(disease ~ test, data = data, family = "binomial")
  }
  else if (covariate == TRUE) {
    fit = glm(disease ~ test + b1 + b2 + c1 + c2 + c3, data = data, family = "binomial")
  }
  preds = prediction(fit, data)$fitted  # Predicted for complete & incomplete data
  snsp_values[1] = sum(data$test*preds) / sum(preds)  # P(T=1|D=1)
  snsp_values[2] = sum((1-data$test)*(1-preds)) / sum(1-preds)  # P(T=0|D=0)
  if (formatted == TRUE) {
    cat("Correction for PVB: Extended Begg and Greenes' Method\n\n")
    cat("Corrected Sensitivity = ", snsp_values[1], "\n")
    cat("Corrected Specificity = ", snsp_values[2])
  } else {
    return(snsp_values)
  }
}
