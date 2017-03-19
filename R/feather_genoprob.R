# feather_genoprob
#' Store genotype probabilities in feather database
#'
#' Uses package feather to convert R object created in R/qtl2geno for fast access.
#'
#' @param genoprob Object of class \code{"calc_genoprob"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}
#' and \code{\link[qtl2geno]{calc_genoprob}}.
#' @param filebase Base of fileame for feather database.
#' @param dirname Directory for feather database.
#'
#' @return A list containing the attributes of \code{genoprob} and the address for the created feather database.
#'
#' @details
#' The genotype probabilities are stored in 1-2 databases
#' as a table of (indivduals*genotypes) x (positions). The first database has all autosome positions, 
#' identified by marker or pseudomarker name. The optional second database is for the X chromosome if present.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' fprobs <- feather_genoprob(probs, "my.feather")

feather_genoprob <-
function(genoprob, filebase, dirname = ".")
{
  # Set up directory for feather objects.
  if(!dir.exists(dirname))
    stop(paste("directory", dirname, "does not exist"))
  
  # Get attributes from genoprob object.
  out <- list(attr = attributes(genoprob))
  
  # Get dimensions and dimnames for chromosome information
  out$chr_dim <- lapply(genoprob, function(x) attributes(x))
  
  # Identify individuals (for later subset use).
  out$ind <- out$chr_dim[[1]]$dimnames[[1]]
  
  is_x_chr <- attr(genoprob, "is_x_chr")
  
  # Add feather addresses
  if(missing(filebase))
    stop("need to supply filebase")
  out$feather <- file.path(dirname, paste(filebase, "feather", sep = "."))
  if(any(is_x_chr))
    out$featherX <- file.path(dirname, paste(filebase, "feather", sep = "_X."))
  
  # Turn list of 3D arrays into table
  # Need to handle X chr separately!
  tbl_array <- function(x) {
    out <- dplyr::tbl_df(matrix(x,, dim(x)[3]))
    names(out) = dimnames(x)[[3]]
    out
  }
  if(any(!is_x_chr)) {
    probs <- sapply(subset(genoprob, chr = !is_x_chr), tbl_array)
    probs <- dplyr::bind_cols(probs)
    feather::write_feather(probs, 
                           out$feather)
  }
  # X matrix probably different size
  if(any(is_x_chr))
    feather::write_feather(tbl_array(genoprob[[which(is_x_chr)]]),
                           out$featherX)
  
  for(a in names(out$attr)[-1])
    attr(out, a) <- out$attr[[a]]
  
  class(out) <- c("feather_genoprob", class(out))
  
  out
}
