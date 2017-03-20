# feather_genoprob_list
#' Extract list element(s) from feather_genoprob object
#'
#' Interprets object as an unclassed list to get components.
#'
#' @param object Object of class \code{\link{feather_genoprob}}.
#' 
#' @param element Element to extract
#' @param ... Additional elements. More than one causes return of list.
#'
#' @return Element or list of elements from \code{\link{feather_genoprob}} object.
#'
#' @details
#' Object is unclassed and elements are extracted. All elements are extracted if only object supplied.
#' See \code{\link{feather_genoprob}} return description for details.
#' 
#' @export
#' @keywords utilities
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
#' fprobs <- feather_genoprob(probs, "my.feather")
#' feather_genoprob_list(fprobs, "feather")
#' sapply(feather_genoprob_list(fprobs, "chr", "ind", "mar"), length)
#' 
feather_genoprob_list <- function(object, element, ...) {
  if(!inherits(object, "feather_genoprob"))
    stop("object must inherit class feather_genoprob")
  
  object <- unclass(object)
  
  if(missing(element))
    return(object)
  
  args <- list(...)
  if(length(args))
    return(object[c(element, unlist(args))])

  object[[element]]
}