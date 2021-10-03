# PVBcorrect

The package contains a number of functions to perform partial verification bias 
(PVB) correction for estimates of accuracy measures in diagnostic accuracy studies. The 
available methods are:

- Begg and Greenes' method (as extended by Alonzo & Pepe, 2005)
- Begg and Greenes' method 1 and 2 (with PPV and NPV as extended by deGroot et al, 2011)
- multiple imputation method by logistic regression (Harel & Zhou, 2006)
- EM-based logistic regression method (Kosinski & Barnhart, 2003)

### Prerequisites

Required packages are:

```
install.packages("boot", "mice")
```

### Installing

Install PVBcorrect package by running:

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

This project is licensed under the GPL (>= 3) License - see the [LICENSE.md](LICENSE.md) file for details.

