**cosimmr - an R package for SIMMs with covariates!**

[cosimmr](https://github.com/emmagovan/cosimmr) is a Bayesian stable isotope mixing model
implemented in R using Fixed Form Variational Bayes. 


If you want the official stable version of the package from CRAN then go to R and type:

```
install.packages('cosimmr')
```

You can then load the package and view either the quick start or the full user manuals with:

```
library(cosimmr)
vignette("cosimmr")
vignette("quick_start")
```

As cosimmr is implemented using [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) when installing you may be warned that the package requires compilation of C/C++/Fortran.