# feather_genoprob
#' Store genotype probabilities in feather database
#'
#' Uses package feather to convert R object created in R/qtl2geno for fast access.
#'
#' @param genoprob Object of class \code{"calc_genoprob"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}
#' and \code{\link[qtl2geno]{calc_genoprob}}.
#' @param fbase Base of fileame for feather database.
#' @param fdir Directory for feather database.
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
function(genoprob, fbase, fdir = ".")
{
  # Set up directory for feather objects.
  if(!dir.exists(fdir))
    stop(paste("directory", fdir, "does not exist"))
  
  # Get attributes from genoprob object.
  attrs <- attributes(genoprob)
  
  # Get dimensions and dimnames for chromosome information
  result <- list(chr_dim = lapply(genoprob, function(x) attributes(x)),
              # Identify individuals (for later subset use).
              chr = names(genoprob))
  result$ind <- result$chr_dim[[1]]$dimnames[[1]]
  
  is_x_chr <- attr(genoprob, "is_x_chr")
  
  # Add feather addresses
  if(missing(fbase))
    stop("need to supply fbase")
  result$feather <- c(A = file.path(fdir, 
                                 paste(fbase, "feather", sep = ".")))
  if(any(is_x_chr))
    result$feather["X"] <- file.path(fdir, paste(fbase, "feather", sep = "_X."))
  
  # Turn list of 3D arrays into table
  # Need to handle X chr separately!
  tbl_array <- function(x) {
    result <- dplyr::tbl_df(matrix(x,, dim(x)[3]))
    names(result) = dimnames(x)[[3]]
    result
  }
  if(any(!is_x_chr)) {
    probs <- sapply(subset(genoprob, chr = !is_x_chr), tbl_array)
    probs <- dplyr::bind_cols(probs)
    feather::write_feather(probs, 
                           result$feather["A"])
  }
  # X matrix probably different size
  if(any(is_x_chr))
    feather::write_feather(tbl_array(genoprob[[which(is_x_chr)]]),
                           result$feather["X"])
  
  # Set up attributes.
  ignore <- match(c("names","class"), names(attrs))
  for(a in names(attrs)[-ignore])
    attr(result, a) <- attrs[[a]]
  
  class(result) <- c("feather_genoprob", attrs$class)
  
  result
}
#' @export
#' @export names.feather_genoprob
#' @method names feather_genoprob
#' 
names.feather_genoprob <- function(x)
  unclass(x)$chr