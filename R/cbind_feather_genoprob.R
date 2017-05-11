#' Join genotype probabilities for different chromosomes
#'
#' Join multiple genotype probability objects, as produced by
#' \code{\link{feather_genoprob}} for different individuals.
#'
#' @param ... Genotype probability objects as produced by
#' \code{\link{feather_genoprob}}. Must have the same set of individuals.
#' @param fbase Base of fileame for feather database. 
#' Needed if objects have different feather databases.
#' @param fdir Directory for feather database.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probsA <- calc_genoprob(grav2[1:5,1:2], map, error_prob=0.002)
#' probsB <- calc_genoprob(grav2[1:5,3:4], map, error_prob=0.002)
#' dir <- tempdir()
#' fprobsA <- feather_genoprob(probsA, "exampleAc", dir)
#' fprobsB <- feather_genoprob(probsB, "exampleBc", dir)
#' fprobs <- cbind(fprobsA, fprobsB, fbase = "exampleABc")
#'
#' @export
#' @export cbind.feather_genoprob
#' @method cbind feather_genoprob
#'
cbind.feather_genoprob <- function(..., fbase, fdir = NULL) {
  # to cbind: probs, is_x_chr
  # to pass through (must match): crosstype, alleles, alleleprobs
      
  bind_feather(list(...),
               check_cbind,
               append_chr,
               cbind,
               fbase, fdir)
}
bind_feather <- function(args, check_fn, append_fn, bind_fn,
                         fbase, fdir = NULL) {
    
    result <- args[[1]]
    if(!inherits(result, "feather_genoprob"))
      stop("argument ", 1, "is not of class feather_genoprob")
    if(length(args) == 1) 
      return(result)

    attrs <- attributes(result)
    result <- unclass(result)
    if(is.null(fdir))
      fdir <- dirname(result$feather)
    if(!dir.exists(fdir))
      stop("directory", fdir, "does not exist")
    
    check_fn(args)
    
    # paste stuff together
    diff_feather <- 0
    for(i in 2:length(args)) {
      # This requires some care, as need to combine feathers
      diff_feather <- (result$feather != unclass(args[[i]])$feather)
      if(diff_feather) {
        diff_feather <- i
        break
      }
      out <- append_fn(result, i, args[[i]], attrs)
      result <- out$result
      attrs <- out$attrs
    }

    if(diff_feather == 2) {
      # Result is just first argument.
      result <- args[[1]]
    } else {
      # Result has cbind of at least on other argument. 
      # Set up attributes.
      ignore <- match(c("names","class"), names(attrs))
      for(a in names(attrs)[-ignore])
        attr(result, a) <- attrs[[a]]
      
      class(result) <- attrs$class
    }

    # Done, unless some args have different feather files.
    if(!diff_feather)
      return(result)
    
    # Different feathers. Need to convert to calc_genoprob and back again.
    if(missing(fbase))
      stop("need to supply fbase to bind distinct feather_genoprob objects")
      
    result <- feather2calc_genoprob(result)
      
    # Convert rest to calc_genoprob and append in usual way.
    for(i in seq(diff_feather, length(args))) {
      argsi <- feather2calc_genoprob(args[[i]])
      result <- bind_fn(result, argsi)
    }

    feather_genoprob(result, fbase, fdir)
}

check_cbind <- function(args) {
  result <- unclass(args[[1]])
  
  for(i in 2:length(args)) {
    if(!inherits(args[[i]], "feather_genoprob"))
      stop("Input object ", i, "is not of class feather_genoprob")
    
    argsi <- unclass(args[[i]])
    if(length(result$ind) != length(argsi$ind) ||
       !all(result$ind == argsi$ind))
      stop("Input objects 1 and ", i, " have different individuals")
  }
}

append_chr <- function(result, i, argsi, attrs) {
  is_x_chr <- attr(argsi, "is_x_chr")
  argsi <- unclass(argsi)

  # Identify what chromosomes are new (if any).  
  new_chr <- is.na(match(argsi$chr, result$chr))
  if(!any(new_chr))
    stop("Input object ", i, "has no new chromosomes")
  if(any(!new_chr))
    warning("duplicate chr ", 
            paste(argsi$chr[!new_chr], collapse = ","),
            " in input object ", i, " ignored")
  new_chr <- argsi$chr[new_chr]
  
  # Identify new markers.
  marnamesi <- unlist(lapply(argsi$dimnames[new_chr], function(x) x[[3]]))
  new_mar <- argsi$mar[argsi$mar %in% marnamesi]

  # Append new components.
  result$chr <- c(result$chr, new_chr)
  result$mar <- c(result$mar, new_mar)
  result$dim <- cbind(result$dim, argsi$dim[, new_chr, drop = FALSE])
  result$dimnames <- c(result$dimnames, argsi$dimnames[new_chr])
  attrs$is_x_chr <- c(attrs$is_x_chr, is_x_chr[new_chr])
  
  list(result = result, attrs = attrs)
}