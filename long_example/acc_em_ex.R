# For sample run, test with low R boot number, will take very long time to finish for large R

# Here we save to R objects for detailed analysis later
em_out = acc_em(data = cad_pvb, test = "T", disease = "D", ci = TRUE, seednum = 12345,
                R = 2, t_max = 1000, cutoff = 0.0005)
em_out$acc_results
em_out$boot_data$t  # bootstrapped data, 1:5 columns are Sn, Sp, PPV, NPV,
                    # t (i.e. EM iteration taken for convergence)
em_out$boot_ci_data

# With covariate, will take some time
# Also check the time taken to finish
start_time = proc.time()
em_outx = acc_em(data = cad_pvb, test = "T", disease = "D", covariate = "X1", ci = TRUE,
                 seednum = 12345, R = 2, t_max = 10000, cutoff = 0.0005)
                 # with covariate, better set larger t_max
elapsed_time = proc.time() - start_time; elapsed_time  # view elapsed time in seconds
em_outx$acc_results
em_outx$boot_data$t
