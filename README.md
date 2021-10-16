# PVBcorrect

The package contains a number of functions to perform partial verification bias 
(PVB) correction for estimates of accuracy measures in diagnostic accuracy studies. The 
available methods are:

- Begg and Greenes' method (as extended by Alonzo & Pepe, 2005)
- Begg and Greenes' method 1 and 2 (with PPV and NPV as extended by deGroot et al, 2011)
- multiple imputation method by logistic regression (Harel & Zhou, 2006)
- EM-based logistic regression method (Kosinski & Barnhart, 2003)

## Prerequisites

The required packages are:

```
install.packages("boot", "mice")
```

## Installation

Install PVBcorrect package by running

```
install.packages("githubinstall")
githubinstall::githubinstall("PVBcorrect")
```

## License

This project is licensed under the GPL (>= 3) License - see the [LICENSE.md](LICENSE.md) file for details.

## Main References

1. Alonzo, T. A., & Pepe, M. S. (2005). Assessing accuracy of a continuous screening test in the presence of verification bias. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 173–190.
2. Begg, C. B., & Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. Biometrics, 207–215.
3. de Groot, J. A. H., Janssen, K. J. M., Zwinderman, A. H., Bossuyt, P. M. M., Reitsma, J. B., & Moons, K. G. M. (2011). Correcting for partial verification bias: a comparison of methods. Annals of Epidemiology, 21(2), 139–148.
4. Harel, O., & Zhou, X.-H. (2006). Multiple imputation for correcting verification bias. Statistics in Medicine, 25(22), 3769–3786.
5. He, H., & McDermott, M. P. (2012). A robust method using propensity score stratification for correcting verification bias for binary tests. Biostatistics, 13(1), 32–47.
6. Kosinski, A. S., & Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. Biometrics, 59(1), 163–171.
7. Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and specificity corrected for verification bias. Communications in Statistics-Theory and Methods, 22(11), 3177–3198.
8. Zhou, X.-H. (1994). Effect of verification bias on positive and negative predictive values. Statistics in Medicine, 13(17), 1737–1745.
9. Zhou, X.-H., Obuchowski, N. A., & McClish, D. K. (2011). Statistical Methods in Diagnostic Medicine (2nd ed.). John Wiley & Sons.
