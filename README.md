# qtl2feather

[![Build Status](https://travis-ci.org/byandell/qtl2feather.svg?branch=master)](https://travis-ci.org/byandell/qtl2feather)

[R/qtl2](http://kbroman.org/qtl2) package using [feather](https://github.com/wesm/feather) to store genotype probabilities

[Karl Broman](http://kbroman.org) & [Brian Yandell](http://www.stat.wisc.edu/~yandell)

R/qtl2feather is a new package that uses [feather](https://github.com/wesm/feather) to store a genotype probabilitis in a feather database for rapid access. It is fully integrated with 
[R/qtl2](http://kbroman.org/qtl2) (aka qtl2). See that package for the bigger story of the qtl2 suite of routines.

---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You first need to install the
[feather](https://github.com/wesm/feather) package and the 
[dplyr](http://dplyr.tidyverse.org/)

    install.packages(c("feather", "dplyr"))

Once you have installed these, install qtl2feather as

    install_github("byandell/qtl2feather")

---

### Vignettes

- [feather_genoprob](https://github.com/byandell/qtl2feather/blob/master/inst/doc/feather_genoprob.Rmd)
- [feather_scan1](https://github.com/byandell/qtl2feather/blob/master/inst/doc/feather_scan1.Rmd)

---

#### License

[Licensed](License.md) under [GPL-3](http://www.r-project.org/Licenses/GPL-3).
