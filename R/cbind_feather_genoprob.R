#' Join genotype probabilities for different chromosomes
#'
#' Join multiple genotype probability objects, as produced by
#' \code{\link{feather_genoprob}} for different individuals.
#'
#' @param ... Genotype probability objects as produced by
#' \code{\link{feather_genoprob}}. Must have the same set of individuals.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probsA <- calc_genoprob(grav2[1:5,1:2], map, error_prob=0.002)
#' probsB <- calc_genoprob(grav2[1:5,3:4], map, error_prob=0.002)
#' probs <- cbind(probsA, probsB)
#'
#' @export
#' @export cbind.feather_genoprob
#' @method cbind feather_genoprob
#'
cbind.feather_genoprob <-
    function(...)
{
    args <- list(..., , basename, dirname = dirname(result$feather["A"]))
    
    # to cbind: probs, is_x_chr
    # to pass through (must match): crosstype, alleles, alleleprobs
    
    result <- args[[1]]
    if(length(args) == 1) return(result)

    attrs <- attributes(result)
    result <- as.list(result)
    
    # paste stuff together
    same_feather <- rep(TRUE, length(args))
    for(i in 2:length(args)) {
      if(!inherits(args[[i]], "feather_genoprob"))
        stop("argument ", i, "is not of class feather_genoprob")

      argsi <- as.list(args[[i]])
      
      if(length(result$ind) != length(argsi$ind) ||
         !all(result$ind == argsi$ind))
        stop("Input objects 1 and ", i, " have different individuals")
      
      # This requires some care, as need to combine feathers
      same_feather[i] <- (result$feather["A"] == argsi$feather["A"])
      if(same_feather[i]) {
        new_chr <- is.na(match(argsi$chr, result$chr))
        if(!any(new_chr))
          stop("argument ", i, "has no new chromosomes")
        new_chr <- argsi$chr[new_chr]
        result$chr <- c(result$chr, new_chr)
        result$chr_dim <- c(result$chr_dim, argsi$chr_dim[new_chr])
        attrs$is_x_chr <- c(attrs$is_x_chr, attr(args[[i]], "is_x_chr")[new_chr])
      }
    }
    
    # Set up attributes.
    ignore <- match(c("names","class"), names(attrs))
    for(a in names(attrs)[-ignore])
      attr(result, a) <- attrs[[a]]
    
    class(result) <- attrs$class

    # Done, unless some args have different feather files.
    if(any(!same_feather)) {
      if(missing(basename))
        stop("need to supply basename to bind distinct feather_genoprob objects")
      
      result <- feather2calc_genoprob(result)
      
      # *** Need to combine feathers here. ***
      # Want to do that once, not every loop.
      # read_feather using $ind and $chr_dim to subset.
      # DO I need feather2calc_genoprob?
      # Similar issue for rbind, but there we have to create new.
      for(i in which(!same_feather)) {
        argsi <- as.list(args[[i]])
        
        new_chr <- is.na(match(argsi$chr, result$chr))
        if(!any(new_chr))
          stop("argument ", i, "has no new chromosomes")
        new_chr <- argsi$chr[new_chr]

        # Append chr by chr
        result <- c(result, )
        
          result$chr <- c(result$chr, new_chr)
          result$chr_dim <- c(result$chr_dim, argsi$chr_dim[new_chr])
          attrs$is_x_chr <- c(attrs$is_x_chr, attr(args[[i]], "is_x_chr")[new_chr])
        }
      }
    }

    # Set up attributes.
    ignore <- match(c("names","class"), names(attrs))
    for(a in names(attrs)[-ignore])
      attr(result, a) <- attrs[[a]]
    
    class(result) <- attrs$class
    
    result
}
