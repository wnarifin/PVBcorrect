# PVBcorrect

The package contains a number of functions to perform partial verification bias 
(PVB) correction for estimates of accuracy measures in diagnostic accuracy studies. The 
available methods are Begg and Greenes' method (extended by Alonzo), 
multiple imputation method (logistic regression) and EM-based logistic regression method.

### Prerequisites

Required packages are

```
install.packages("prediction", "mice", "simstudy")
```

### Installing

Install the package by running

```
install.packages("devtools")
devtools::install_github("wnarifin/PVBcorrect")
```

or

```
install.packages("githubinstall")
githubinstall::githubinstall("PVBcorrect")
```

## License

This project is licensed under the GPL (>= 3) License - see the [LICENSE.md](LICENSE.md) file for details

