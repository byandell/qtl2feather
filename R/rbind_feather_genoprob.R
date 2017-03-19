#' Join genotype probabilities for different individuals
#'
#' Join multiple genotype probability objects, as produced by
#' \code{\link{feather_genoprob}} for different individuals.
#'
#' @param ... Genotype probability objects as produced by
#' \code{\link{feather_genoprob}}. Must have the same set of markers and
#' genotypes.
#' @param basename Base of fileame for feather database. 
#' Needed if objects have different feather databases.
#' @param dirname Directory for feather database.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probsA <- calc_genoprob(grav2[1:5,], map, error_prob=0.002)
#' probsB <- calc_genoprob(grav2[6:12,], map, error_prob=0.002)
#' probs <- rbind(probsA, probsB)
#'
#' @export
#' @export rbind.feather_genoprob
#' @method rbind feather_genoprob
#' 
rbind.feather_genoprob <-
    function(..., basename, dirname = dirname(result$feather["A"]))
{
    args <- list(...)
    
    # to rbind: the data
    # to pass through (must match): crosstype, is_x_chr, alleles, alleleprobs
    
    result <- args[[1]]
    if(length(args) == 1) return(result)
    
    # check that things match
    other_stuff <- c("crosstype", "is_x_chr", "alleles", "alleleprobs")
    for(i in 2:length(args)) {
      if(!inherits(args[[i]], "feather_genoprob"))
        stop("argument ", i, "is not of class feather_genoprob")
      for(obj in other_stuff) {
        if(!is_same(attr(args[[1]], obj), attr(args[[i]], obj)))
          stop("Input objects 1 and ", i, " differ in their ", obj)
      }
    }
    
    attrs <- attributes(result)
    result <- unclass(result)
    
    # paste stuff together
    diff_feather <- FALSE
    for(i in 1:length(args)) {
      argsi <- unclass(args[[i]])
      if(!is_same(names(args[[1]]), names(args[[i]])))
        stop("Input objects 1 and ", i, " have different chromosome names")
      for(chr in names(args[[1]])) {
        dimn1 <- dimnames(args[[1]][[chr]])
        dimni <- dimnames(args[[i]][[chr]])
        if(!is_same(dimn1[-1], dimni[-1]))
          stop("Input objects 1 and ", i, " differ in shape on chromosome ", chr)
        
        result[[chr]][index[[i]],,] <- args[[i]][[chr]]
        rownames(result[[chr]])[index[[i]]] <- rownames(args[[i]][[chr]])
      }
    }
    
    # paste in the attributes
    for(obj in other_stuff)
      attr(result, obj) <- attr(args[[1]], obj)
    class(result) <- class(args[[1]])
    
    result
}
