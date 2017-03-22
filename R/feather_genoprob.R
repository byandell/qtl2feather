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
#' @param verbose Show warning of feather creation if \code{TRUE} (default).
#'
#' @return A list containing the attributes of \code{genoprob}
#' and the address for the created feather database. 
#' Components are:
#' \itemize{
#' \item \code{dim} - List of all dimensions of 3-D arrays.
#' \item \code{dimnames} - List of all dimension names of 3-D arrays.
#' \item \code{chr} - Vector of (subset of) chromosome names for this object.
#' \item \code{ind} - Vector of (subset of) individual names for this object.
#' \item \code{mar} - Vector of (subset of) marker names for this object.
#' \item \code{feather} - Names of feather databases: \code{A} = autosome database, \code{X} = X chromosome database (if needed).
#' }
#'
#' @details
#' The genotype probabilities are stored in 1-2 databases
#' as a table of (indivduals*genotypes) x (positions). The first database has all autosome positions,
#' identified by marker or pseudomarker name. The optional second database is for the X chromosome if present.
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
#' fprobs <- feather_genoprob(probs, "my.feather")

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
                 dimnames = lapply(chr_dim, function(x) x$dimnames))

  # Identify chromosome, individuals, markers (for later subset use).
  result$chr <- names(genoprob)
  result$ind <- result$dimnames[[1]][[1]]
  result$mar <- as.vector(unlist(lapply(result$dimnames, function(x) x[[3]])))

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
    dims <- dim(x)
    dnames <- dimnames(x)
    dim(x) <- c(prod(dims[1:2]), dims[3])
    dimnames(x) <- list(NULL, dnames[[3]])
    as.data.frame(x)
  }
  if(any(!is_x_chr)) {
    probs <- lapply(subset(genoprob, chr = !is_x_chr), tbl_array)
    probs <- dplyr::bind_cols(probs)
    if(verbose) message("writing ", result$feather["A"])
    feather::write_feather(probs,
                           result$feather["A"])
  }
  # X matrix probably different size
  if(any(is_x_chr)) {
    probs <- tbl_array(genoprob[[which(is_x_chr)]])
    if(verbose) message("writing ", result$feather["X"])
    feather::write_feather(probs,
                           result$feather["X"])
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
