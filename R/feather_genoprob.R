# feather_genoprob
#' Store genotype probabilities in feather database
#'
#' Uses package feather to convert R object created in R/qtl2geno for fast access.
#'
#' @param genoprob Object of class \code{"calc_genoprob"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}
#' and \code{\link[qtl2geno]{calc_genoprob}}.
#' @param fbase Base of filename for feather database.
#' @param fdir Directory for feather database.
#' @param verbose Show warning of feather creation if \code{TRUE} (default).
#'
#' @return A list containing the attributes of \code{genoprob}
#' and the address for the created feather database.
#' Components are:
#' \itemize{
#' \item \code{dim} - List of all dimensions of 3-D arrays.
#' \item \code{dimnames} - List of all dimension names of 3-D arrays.
#' \item \code{is_x_chr} - Vector of all is_x_chr attributes.
#' \item \code{chr} - Vector of (subset of) chromosome names for this object.
#' \item \code{ind} - Vector of (subset of) individual names for this object.
#' \item \code{mar} - Vector of (subset of) marker names for this object.
#' \item \code{feather} - Names of feather databases: \code{A} = autosome database, \code{X} = X chromosome database (if needed).
#' }
#'
#' @details
#' The genotype probabilities are stored in separate databases for each chromosome
#' as tables of (indivduals*genotypes) x (positions) in directory \code{feather}. 
#' The \code{dim}, \code{dimnames} and \code{is_x_chr} elements of the object
#' have information about the entire feather database.
#' If a \code{feather_genoprob} object is a subset of another such object, 
#' the \code{chr}, \code{ind}, and \code{mar} contain information about what is in the subset.
#' However, the \code{feather} databases are not altered in a subset, and can be restored by
#' \code{\link{feather_genoprob_restore}}. The actual elements of a \code{feather_genoprob}
#' object are only accessible to the user after a call to \code{\link[base]{unclass}}; instead
#' the usual access to elements of the object invoke \code{\link{subset.feather_genoprob}}.
#'
#' @importFrom feather write_feather
#' @importFrom dplyr bind_cols
#' @export
#' @keywords utilities
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
#' dir <- tempdir()
#' fprobs <- feather_genoprob(probs, "grav2", dir)

feather_genoprob <-
function(genoprob, fbase, fdir = ".", verbose = TRUE)
{
  # Set up directory for feather objects.
  if(!dir.exists(fdir))
    stop(paste("directory", fdir, "does not exist"))

  # Get attributes from genoprob object.
  attrs <- attributes(genoprob)

  # Get dimensions and dimnames for chromosome information
  chr_dim <- lapply(genoprob, function(x) attributes(x))
  result <- list(dim = sapply(chr_dim, function(x) x$dim),
                 dimnames = lapply(chr_dim, function(x) x$dimnames),
                 is_x_chr = attr(genoprob, "is_x_chr"))

  # Identify chromosome, individuals, markers (for later subset use).
  result$chr <- names(genoprob)
  result$ind <- result$dimnames[[1]][[1]]
  tmp <- unlist(lapply(result$dimnames, function(x) x[[3]]))
  names(tmp) <- NULL
  result$mar <- tmp

  # Add feather addresses
  if(missing(fbase))
    stop("need to supply fbase")
  result$feather <- file.path(fdir, fbase)

  # Turn list of 3D arrays into table
  # Need to handle X chr separately!
  tbl_array <- function(x) {
    dims <- dim(x)
    dnames <- dimnames(x)
    dim(x) <- c(prod(dims[1:2]), dims[3])
    dimnames(x) <- list(NULL, dnames[[3]])
    as.data.frame(x)
  }
  for(chr in result$chr) {
    probs <- tbl_array(genoprob[[chr]])
    fname <- paste0(result$feather, "_", chr, ".feather")
    if(file.exists(fname))
      warning("writing over existing ", fname)
    if(verbose) message("writing ", fname)
    feather::write_feather(probs, fname)
  }

  # Set up attributes.
  ignore <- match(c("names","class"), names(attrs))
  for(a in names(attrs)[-ignore])
    attr(result, a) <- attrs[[a]]

  class(result) <- c("feather_genoprob", attrs$class)

  result
}
#' @export
names.feather_genoprob <- function(x)
  unclass(x)$chr
#' @export
length.feather_genoprob <- function(x)
  length(unclass(x)$chr)
#' @export
dim.feather_genoprob <- function(x) {
  x <- unclass(x)
  out <- x$dim[, x$chr, drop = FALSE]
  rownames(out) <- c("ind","gen","mar")
  out[1,] <- length(x$ind)
  out[2,] <- sapply(index_chr(x$dimnames[x$chr], 2),
                    length)
  out[3,] <- sapply(index_chr(x$dimnames[x$chr], 3, x$mar),
                    length)
  out
}
#' @export
dimnames.feather_genoprob <- function(x) {
  x <- unclass(x)
  dnames <- x$dimnames
  out <- list(ind = x$ind,
              gen = index_chr(dnames[x$chr], 2),
              mar = index_chr(dnames[x$chr], 3, x$mar))
  out
}
index_chr <- function(dnames, index, index_sub=NULL) {
  lapply(dnames, function(x, index, index_sub) {
      xnames <- x[[index]]
      if(!is.null(index_sub))
        xnames <- xnames[xnames %in% index_sub]
      xnames
    },
    index, index_sub)
}
