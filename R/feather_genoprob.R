# feather_genoprob
#' Store genotype probabilities in feather database
#'
#' Uses package feather to convert R object created in R/qtl2geno for fast access.
#'
#' @param genoprob Object of class \code{"calc_genoprob"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}
#' and \code{\link[qtl2geno]{calc_genoprob}}.
#' @param filename Name of file for feather database.
#'
#' @return A list containing the following
#' \itemize{
#' \item \code{probs} - An object with information to access feather database.
#' \item \code{map} - The genetic map as a list of vectors of marker positions.
#' \item \code{grid} - A list of logical vectors, indicating which
#'     positions correspond to a grid of markers/pseudomarkers. (may be
#'     absent)
#' \item \code{crosstype} - The cross type of the input \code{cross}.
#' \item \code{is_x_chr} - Logical vector indicating whether chromosomes
#'     are to be treated as the X chromosome or not, from input \code{cross}.
#' \item \code{is_female} - Vector of indicators of which individuals are female, from input
#'     \code{cross}.
#' \item \code{cross_info} - Matrix of cross information for the
#'     individuals, from input \code{cross}.
#' \item \code{alleles} - Vector of allele codes, from input
#'     \code{cross}.
#' \item \code{alleleprobs} - Logical value (\code{FALSE}) that
#'     indicates whether the probabilities are compressed to allele
#'     probabilities, as from \code{\link{genoprob_to_alleleprob}}.
#' \item \code{step} - the value of the \code{step} argument.
#' \item \code{off_end} - the value of the \code{off_end} argument.
#' \item \code{stepwidth} - the value of the \code{stepwidth} argument.
#' \item \code{error_prob} - the value of the \code{error_prob} argument.
#' \item \code{map_function} - the value of the \code{map_function} argument.
#' }
#'
#' @details
#' To be added here
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
function(genoprob, filename)
{
  ## fill in code here
  genoprob
}
