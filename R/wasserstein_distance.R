#' Compute Lp Wasserstein Distance 
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param q internal parameter
#' @param p The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @return The value for the Lp Wassesterstein computed in (\code{d1}, \code{d2})
#' @export
wasserstein.distance = function(d1, d2, q, p = 2) {
  
  wasserstein_distance(Diag1 = d1, Diag2 = d2, q = q, internal_p = p, 
                       delta = 0.01, initial_eps = 1.0, eps_factor = 5.0)
  
}