L2.wass <- function(d1, d2, h, dimension, pp){
  
  d1 = d1[d1[,1]==dimension, 2:3, drop = F]
  d2 = d2[d2[,1]==dimension, 2:3, drop = F]
  
  return(exp(-1/h * wasserstein_distance(d1, d2, q = pp,internal_p = 2)))
  
}

#' Compute L_{2} Geodesic Gaussian Kernel (GGK)
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time) or a list of diagrams
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param h bandwidth of the kernel
#' @param dimension The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @param p order of the p-Wasserstein distance
#' @return If \code{d1} is a list of Persistence Diagrams, this function returns a matrix whose (i,j) entry is the GGK computed in (\code{d1}[[i]], \code{d2}[[j]]), 
#' otherwise it returns the value for the GGK computed in (\code{d1}, \code{d2}).
#' @references 
#' @export
gaus.kernel <- function (d1, d2 = NULL, h, dimension, p)
{
  if(!is.null(d2)) {
    out = L2.wass(d1, d2, h = h, dimension=dimension, pp = p)
  }
  else{
    k.fun = function(x, y) L2.wass(d1[[x]], d1[[y]], h = h, dimension=dimension, pp = p)
    k.fun = Vectorize(k.fun)
    d.idx = seq_along(d1)
    out   = outer(d.idx,d.idx, k.fun)
  }
  return(out)
}


#' Compute L_{2} Geodesic Laplacian Kernel (GLK)
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time) or a list of diagrams
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param h bandwidth of the kernel
#' @param dimension The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @param p order of the p-Wasserstein distance
#' @return If \code{d1} is a list of Persistence Diagrams, this function returns a matrix whose (i,j) entry is the GLK computed in (\code{d1}[[i]], \code{d2}[[j]]),
#' otherwise it returns the value for the GLK computed in (\code{d1}, \code{d2}).
#' @export
lapl.kernel <- function (d1, d2 = NULL, h, dimension, p)
{
  if(!is.null(d2)) {
    out = L2.wass(d1, d2, h = h, dimension=dimension, pp = p)
  }
  else{
    k.fun = function(x, y) L2.wass(d1[[x]], d1[[y]], h = h, dimension=dimension, pp = p)
    k.fun = Vectorize(k.fun)
    d.idx = seq_along(d1)
    out   = outer(d.idx,d.idx, k.fun)
  }
  return(out)
}



# Linfty Wasserstein Kernel -----------------------------------------------

#' Compute L_{\infty}} Wasserstein Kernel (WK)
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time) or a list of diagrams
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param h bandwidth of the kernel
#' @param dimension The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @param p order of the p-Wasserstein distance
#' @return If \code{d1} is a list of Persistence Diagrams, this function returns a matrix whose (i,j) entry is the WK computed in (\code{d1}[[i]], \code{d2}[[j]]),
#' otherwise it returns the value for the WK computed in (\code{d1}, \code{d2}).
wass.kernel <- function (d1, d2 = NULL, h, dimension, p)
{
  
  
  if(!is.null(d2)) {
    out = exp(-1/h*TDA::wasserstein(d1, d2, p=p, dimension = dimension)^2)
  }
  else{
    k.fun = function(x, y) exp(-1/h*TDA::wasserstein(d1[[x]], d1[[y]], p=p, dimension = dimension)^2)
    k.fun = Vectorize(k.fun)
    d.idx = seq_along(d1)
    out   = outer(d.idx,d.idx, k.fun)
  }
  return(out)
}



# Persistence Scale Space Kernel ------------------------------------------


pss.k <- function(x, y, h, dimension)
{
  d1 = x[x[,1]==dimension, drop = F ]
  d2 = y[y[,1]==dimension, drop = F ]
  
  kvect = 0
  
  if(is.null(nrow(d1))) d1 = matrix(d1, ncol = 3)
  if(is.null(nrow(d2))) d2 = matrix(d2, ncol = 3)
  
  if(nrow(d1)==0 || nrow(d2)==0) return(0)
  else {
    for(i in 1:nrow(d1)){
      for(j in 1:nrow(d2)){
        exp.dist1 = sum((d1[i,2:3]-d2[j,2:3])^2)
        exp.dist2 = sum((d1[i,2:3]-d2[j,3:2])^2)
        sum.elem = exp(-1/(8*h) * exp.dist1) - exp(-1/(8*h) * exp.dist2)
        kvect = kvect + sum.elem
      }
    }
    return(1/(pi*8*h) * kvect)
  }
}

#' Compute Persistence Scale Space Kernel (PSSK)
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time) or a list of diagrams
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param h bandwidth of the kernel
#' @param dimension The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @return If \code{d1} is a list of Persistence Diagrams, this function returns a matrix whose (i,j) entry is the PSSK computed in (\code{d1}[[i]], \code{d2}[[j]]),
#' otherwise it returns the value for the PSSK computed in (\code{d1}, \code{d2})
#' @export
pss.kernel <- function (d1, d2 = NULL, h, dimension)
{
  if(!is.null(d2)) {
    out = pss.k(d1, d2, dimension=dimension, h=h)
  }
  else{
    k.fun = function(x, y) pss.k(d1[[x]], d1[[y]], h = h, dimension = dimension)
    k.fun = Vectorize(k.fun)
    d.idx = seq_along(d1)
    out   = outer(d.idx,d.idx, k.fun)
  }
  return(out)
}



# Sliced Wasserstein Kernel -----------------------------------------------


sw.k <- function(x, y, dimension, h, M){
  
  
  d1 = matrix(x[x[,1]==dimension,], ncol = 3)
  d2 = matrix(y[y[,1]==dimension,], ncol = 3)
  
  
  p1 = (d1[,3]-d1[,2])/2
  p2 = (d2[,3]-d2[,2])/2
  
  n1 = nrow(d1)
  n2 = nrow(d2)
  
  dd1 = rbind( d1, cbind(rep(dimension, n2), p2, p2))
  dd2 = rbind( d2, cbind(rep(dimension, n1), p1, p1))
  
  angles = seq(-pi/2, pi/2, length.out = M +1)
  step.theta = pi/M
  
  theta = cbind(cos(angles), sin(angles))
  
  
  V1 = matrix(NA, nrow = n1+n2, ncol = M)
  V2 = matrix(NA, nrow = n1+n2, ncol = M)
  
  for(m in 1:M){
    V1[, m] = apply(dd1[,2:3], 1, FUN= crossprod, y = theta[m,])
    V2[, m] = apply(dd2[,2:3], 1, FUN= crossprod, y = theta[m,])
  }
  
  VV1 = apply(V1, 2, sort, decreasing = F)
  VV2 = apply(V2, 2, sort, decreasing = F)
  
  SW = step.theta * sum(abs(VV1 - VV2))
  
  return(exp( - pi*SW/h))
}


#' Compute Persistence Sliced Wasserstein Kernel (SWK)
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time) or a list of diagrams
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param h bandwidth of the kernel
#' @param dimension The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @param M number of directions on which to approximate the Sliced Wasserstein Distance 
#' @return If \code{d1} is a list of Persistence Diagrams, this function returns a matrix whose (i,j) entry is the SWK computed in (\code{d1}[[i]], \code{d2}[[j]]),
#' otherwise it returns the value for the SWK computed in (\code{d1}, \code{d2})
#' @export
sw.kernel <- function(d1, d2 = NULL, h, dimension, M=10){
  
  if(!is.null(d2)) {
    out = sw.k(d1, d2, dimension=dimension, h=h, M=M)
  }
  else{
    k.fun = function(x, y) sw.k(d1[[x]], d1[[y]], h = h, dimension = dimension, M = M)
    k.fun = Vectorize(k.fun)
    d.idx = seq_along(d1)
    out   = outer(d.idx,d.idx, k.fun)
  }
  return(out)
  
}


# Persistence Fisher Kernel -----------------------------------------------

pf.k <- function(x, y, dimension, h, sigma)
{
  d1 = x[x[,1]==dimension, 2:3]
  d2 = y[y[,1]==dimension, 2:3]
  
  
  if(is.null(nrow(d1))) d1 = matrix(d1, ncol = 2)
  if(is.null(nrow(d2))) d2 = matrix(d2, ncol = 2)
  
  projd1 = apply(d1, 1, mean)
  projd2 = apply(d2, 1, mean)
  
  dd1 = rbind(d1, cbind(projd2,projd2))
  dd2 = rbind(d2, cbind(projd1,projd1))
  dd  = rbind(dd1, dd2)
  
  wght = 1/nrow(dd1)
  
  f1 = apply(dd1, 1, function(x) mvtnorm::dmvnorm(dd, mean = x, sigma = sigma*diag(2) ))
  f2 = apply(dd2, 1, function(x) mvtnorm::dmvnorm(dd, mean = x, sigma = sigma*diag(2) ))
  
  ff1 = apply(wght * f1, 1, sum)
  ff1 = sqrt(ff1/sum(ff1))
  
  ff2 = apply(wght * f2, 1, sum)
  ff2 = sqrt(ff2/sum(ff2))
  
  fd = acos(ff1%*%ff2)

 return(exp(- 1/h * fd))

}


#' Compute Persistence Fisher Kernel (PFK)
#'
#' @param d1 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time) or a list of diagrams
#' @param d2 A persistence diagram (matrix with 3 col where the first one is the dimension, the second is the birth-time and the third is the death-time).
#' @param h bandwidth of the PFK
#' @param dimension The dimension of the topological feature (0 for connected components, 1 for cycles etc)
#' @param sigma standard deviation of Gaussian Kernel used to 
#' @return If \code{d1} is a list of Persistence Diagrams, this function returns a matrix whose (i,j) entry is the PFK computed in (\code{d1}[[i]], \code{d2}[[j]]),
#' otherwise it returns the value for the PFK computed in (\code{d1}, \code{d2})
#' @export
pf.kernel <- function(d1, d2 = NULL, h, dimension, sigma){
  if(!is.null(d2)) {
    out = pf.k(d1, d2, h = h, dimension=dimension, sigma = sigma)
  }
  else{
    k.fun = function(x, y) pf.k(d1[[x]], d1[[y]], h = h, dimension=dimension, sigma = sigma)
    k.fun = Vectorize(k.fun)
    d.idx = seq_along(d1)
    out   = outer(d.idx,d.idx, k.fun)
  }
  return(out)
  
  
}


