# feather2fst_genoprob
#' Migrate feather_genoprob to fst_genoprob database.
#'
#' See package \code{R/qtl2fst} for faster database version.
#'
#' @param fprobs object of class \code{feather_genoprob}.
#' @param verbose Show warning of feather creation if \code{TRUE} (default).
#'
#' @return An object of class \code{fst_genoprob}.
#'
#' @details
#' The genotype probabilities are stored using package \code{fst} rather package \code{feather} for faster access.
#' If the feather directory address is relative, then be sure to use \code{\link{setwd}} to get to the parent directory.
#'
#' @importFrom feather read_feather
#' @importFrom fst write_fst
#' @export
#' @keywords utilities
#'
feather2fst_genoprob <- function(fprobs, verbose = TRUE) {
  attrs <- attributes(fprobs)
  fp <- unclass(fprobs)
  fp$fst <- fp$feather
  fp$feather <- NULL
  
  # Read feather and write fst.
  # Don't use result$chr as this may be subset.
  for(chr in names(fp$dimnames)) {
    fname <- paste0(fp$fst, "_", chr)
    if(verbose) message("reading ", fname, " feather, writing fst")
    fst::write_fst(
      feather::read_feather(
        paste0(fname, ".feather")),
      paste0(fname, ".fst"))
  }
  
  # Set up attributes.
  ignore <- match(c("names","class"), names(attrs))
  for(a in names(attrs)[-ignore])
    attr(fp, a) <- attrs[[a]]
  
  class(fp) <- attrs$class
  
  fp
}