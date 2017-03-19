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
    args <- list(...)
    
    # to cbind: probs, is_x_chr
    # to pass through (must match): crosstype, alleles, alleleprobs
    
    result <- args[[1]]
    if(length(args) == 1) return(result)
    
    # paste stuff together
    for(i in 2:length(args)) {
      if(length(args[[1]]$ind) != length(args[[i]]$ind) ||
         !all(args[[1]]$ind == args[[i]]$ind))
        stop("Input objects 1 and ", i, " have different individuals")
      
#     result <- c(result, args[[i]])
      
      # This requires some care, as need to combine feathers
      if(result$feather == args[[i]]$feather) {
        result$chr_dim <- c(result$chr_dim, args[[i]]$chr_dim)
        results$attr$is_x_chr <- c(results$attr$is_x_chr, args[[i]]$attr$is_x_chr)
      }
      else {
        stop("cannot handle this yet")
        
        # Need to combine feathers here.
        # Want to do that once, not every loop.
        # read_feather using $ind and $chr_dim to subset.
        # DO I need feather2calc_genoprob?
        # Similar issue for rbind, but there we have to create new.
      }
    }
    
    other_stuff <- c("is_x_chr")
    for(obj in other_stuff)
      attr(result, obj) <- attr(args[[1]], obj)
    for(i in 2:length(args)) {
      for(obj in other_stuff) {
        if(is.null(attr(args[[1]], obj)) && is.null(attr(args[[i]], obj))) next # not present
        if(is.null(attr(args[[1]], obj)) || is.null(attr(args[[i]], obj)))
          stop(obj, " not present in all inputs")
        attr(result, obj) <- c(attr(result, obj), attr(args[[i]], obj))
      }
    }
    
    # check that things match
    other_stuff <- c("crosstype", "alleles", "alleleprobs")
    for(i in 2:length(args)) {
      for(obj in other_stuff) {
        if(!is_same(attr(args[[1]], obj), attr(args[[i]], obj)))
          stop("Input objects 1 and ", i, " differ in their ", obj)
      }
    }
    for(obj in other_stuff)
      attr(result, obj) <- attr(args[[1]], obj)
    
    class(result) <- class(args[[1]])
    
    result
}
