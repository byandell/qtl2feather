% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cbind_feather_genoprob.R
\name{cbind.feather_genoprob}
\alias{cbind.feather_genoprob}
\title{Join genotype probabilities for different chromosomes}
\usage{
\method{cbind}{feather_genoprob}(..., fbase, fdir = NULL)
}
\arguments{
\item{...}{Genotype probability objects as produced by
\code{\link{feather_genoprob}}. Must have the same set of individuals.}

\item{fbase}{Base of fileame for feather database. 
Needed if objects have different feather databases.}

\item{fdir}{Directory for feather database.}
}
\value{
A single genotype probability object.
}
\description{
Join multiple genotype probability objects, as produced by
\code{\link{feather_genoprob}} for different individuals.
}
\examples{
library(qtl2)
grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
map <- insert_pseudomarkers(grav2$gmap, step=1)
probsA <- calc_genoprob(grav2[1:5,1:2], map, error_prob=0.002)
probsB <- calc_genoprob(grav2[1:5,3:4], map, error_prob=0.002)
dir <- tempdir()
fprobsA <- feather_genoprob(probsA, "exampleAc", dir)
fprobsB <- feather_genoprob(probsB, "exampleBc", dir)
fprobs <- cbind(fprobsA, fprobsB, fbase = "exampleABc")

}
