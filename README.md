
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kernelTDA

<!-- badges: start -->

<!-- badges: end -->

This `R`-package provides an implementation of the most famous kernels
found in the framework of Topological Data Analysis (TDA), more
specifically:

  - [Persistence Scale Space
    Kernel](http://openaccess.thecvf.com/content_cvpr_2015/papers/Reininghaus_A_Stable_Multi-Scale_2015_CVPR_paper.pdf)
  - [Sliced Wasserstein
    Kernel](https://dl.acm.org/citation.cfm?id=3305450)
  - [Persistence Fisher
    Kernel](http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams.pdf)
  - [Geodesic Wasserstein Kernel(s)](https://arxiv.org/abs/1709.07100)
  - [Persistence
    Images](http://www.jmlr.org/papers/volume18/16-337/16-337.pdf)

In addition, it also provides an `R` interface to the C++ libray
[HERA](https://bitbucket.org/grey_narn/hera/src/master/), which contains
an efficient implementation of the \(L_p\) q-Wasserstein distance
between persistence diagrams.

This package is not yet on CRAN, yet you can install it from this
repository with:

``` r
# install.packages("remotes")

remotes::install_github("tulliapadellini/kernelTDA")
```
