# PVBcorrect

The package contains a number of functions to perform partial verification bias 
(PVB) correction for estimates of accuracy measures in diagnostic accuracy studies. The 
available methods are:

- Begg and Greenes' method (as extended by Alonzo & Pepe, 2005)
- Begg and Greenes' method 1 and 2 (with PPV and NPV as extended by deGroot et al, 2011)
- Multiple imputation method by logistic regression (Harel & Zhou, 2006)
- EM-based logistic regression method (Kosinski & Barnhart, 2003)
- Inverse Probability Bootstrap (IPB) sampling method (Arifin & Yusof, 2022; Nahorniak et al., 2015)
- Scaled Inverse Probability Resampling methods (Arifin & Yusof, 2023; Arifin & Yusof, 2025)

## Prerequisites

The required packages are:

```
install.packages("boot", "mice")
```

## Installation

Install PVBcorrect package by running

```
install.packages("devtools")
devtools::install_github("wnarifin/PVBcorrect")
```

## Usage, news and updates

Please view Wiki page: https://github.com/wnarifin/PVBcorrect/wiki

## License

This project is licensed under the GPL (>= 3) License - see the [LICENSE.md](LICENSE.md) file for details.

## References

1. Alonzo, T. A., & Pepe, M. S. (2005). Assessing accuracy of a continuous screening test in the presence of verification bias. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 173–190.
2. Arifin, W. N., & Yusof, U. K. (2025). Partial Verification Bias Correction Using Scaled Inverse Probability Resampling for Binary Diagnostic Tests. medRxiv. https://doi.org/10.1101/2025.03.09.25323631
3. Arifin, W. N. (2023). Partial verification bias correction in diagnostic accuracy studies using propensity score-based methods (PhD thesis, Universiti Sains Malaysia). https://erepo.usm.my/handle/123456789/19184
4. Arifin, W. N., & Yusof, U. K. (2022a). Correcting for partial verification bias in diagnostic accuracy studies: a tutorial using R. Statistics in Medicine, 41(9), 1709–1727.
5. Arifin, W. N., & Yusof, U. K. (2022b). Partial Verification Bias Correction Using Inverse Probability Bootstrap Sampling for Binary Diagnostic Tests. Diagnostics, 12, 2839.
6. Begg, C. B., & Greenes, R. A. (1983). Assessment of diagnostic tests when disease verification is subject to selection bias. Biometrics, 207–215.
7. de Groot, J. A. H., Janssen, K. J. M., Zwinderman, A. H., Bossuyt, P. M. M., Reitsma, J. B., & Moons, K. G. M. (2011). Correcting for partial verification bias: a comparison of methods. Annals of Epidemiology, 21(2), 139–148.
8. Harel, O., & Zhou, X.-H. (2006). Multiple imputation for correcting verification bias. Statistics in Medicine, 25(22), 3769–3786.
9. He, H., & McDermott, M. P. (2012). A robust method using propensity score stratification for correcting verification bias for binary tests. Biostatistics, 13(1), 32–47.
10. Kosinski, A. S., & Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. Biometrics, 59(1), 163–171.
11. Nahorniak, M., Larsen, D. P., Volk, C., & Jordan, C. E. (2015). Using Inverse Probability Bootstrap Sampling to Eliminate Sample Induced Bias in Model Based Analysis of Unequal Probability Samples. Plos One, 10(6), e0131765. https://doi.org/10.1371/journal.pone.0131765
12. Zhou, X.-H. (1993). Maximum likelihood estimators of sensitivity and specificity corrected for verification bias. Communications in Statistics-Theory and Methods, 22(11), 3177–3198.
13. Zhou, X.-H. (1994). Effect of verification bias on positive and negative predictive values. Statistics in Medicine, 13(17), 1737–1745.
13. Zhou, X.-H., Obuchowski, N. A., & McClish, D. K. (2011). Statistical Methods in Diagnostic Medicine (2nd ed.). John Wiley & Sons.
