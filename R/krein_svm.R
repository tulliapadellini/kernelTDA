#' Krein Support Vector Machine
#' 
#' Solves a kernelized Support Vector Machine in the case where the kernel used may not be positive semidefinite.
#'
#' @param kernelmat the kernel matrix computed for all observations
#' @param ... additional parameters, see \code{\link{krein.svm.default}} for more details on how to use this function
#' @author Tullia Padellini, Francesco Palini, David Meyer. The included C++ library LIBSVM is authored by Chih-Chung Chang and Chih-Jen Lin)
#' @details This function implements the Krein Support Vector Machine solver as defined by Loosli et al. (2015). 
#' The implementation of the solver is a modified version of the popular C++ library `LIBSVM`, while the connection to `R` 
#' heavily relies on the `R`-package \pkg{e1701}. 
#' @examples 
#' ## DO NOT RUN:
#' # library(TDA)
#' # set.seed(123)
#' # foo.data = list()
#' # for(i in 1:20){
#' #    foo = circleUnif(100)
#' #    foo.data[[i]] = ripsDiag(foo, 1,1)$diagram}
#' #    for(i in 21:40){   
#' #     foo = cbind(runif(100), runif(100))
#' #     foo.data[[i]] = ripsDiag(foo, 1,1)$diagram
#' #     }
#' # GSWkernel = gaus.kernel(foo.data, h =1, dimension = 1,  q = 2)
#' # GGKclass = krein.svm(kernelmat = GSWkernel, y = rep(c(1,2), c(20,20)))
#' @return An object of class \code{krein.svm} containing the fitted model, including: 
#' \describe{
#'   \item{\code{SV}}{a matrix containing the Support Vectors}
#'   \item{\code{index}}{index of the resulting support vectors in the data matrix}
#'   \item{\code{coefs}}{a matrix containing corresponding coefficients times the training labels}
#'   \item{\code{rho}}{value of the (negative) intercept}
#' }
#' @references 
#' \insertRef{loosli2015learning}{kernelTDA}
#' 
#' \insertRef{chang2011libsvm}{kernelTDA}
#' 
#' \insertRef{dimitriadou2008misc}{kernelTDA}
#' @export
krein.svm <-
function (kernelmat, ...)
    UseMethod ("krein.svm")

#' Krein Support Vector Machine
#' 
#' Solves a kernelized Support Vector Machine in the case where the kernel used may not be positive semidefinite.
#'
#' @param kernelmat the kernel matrix computed for all observations
#' @param y a vector of labels
#' @param cost cost of violating the constraint
#' @param class.weights a named vector of weights for the different classes, used for asymmetric class sizes. Not all factor levels have to be supplied (default weight: 1). All components have to be named. Specifying "inverse" will choose the weights inversely proportional to the class distribution.
#' @param cross number of fold in a k-fold cross validation 
#' @param probability logical indicating whether the model should allow for probability predictions (default: \code{FALSE}).
#' @param fitted logical indicating whether the fitted values should be computed and included in the model or not (default: \code{TRUE})
#' @param subset an index vector specifying the cases to be used in the training sample. (NOTE: If given, this argument must be named.)
#' @param ... additional parameters
#' @author Tullia Padellini, Francesco Palini, David Meyer. The included C++ library LIBSVM is authored by Chih-Chung Chang and Chih-Jen Lin)
#' @details This function implements the Krein Support Vector Machine solver as defined by Loosli et al. (2015). 
#' The implementation of the solver is a modified version of the popular C++ library `LIBSVM`, while the connection to `R` 
#' heavily relies on the `R`-package \pkg{e1701}. 
#' @examples 
#' ## DO NOT RUN:
#' # library(TDA)
#' # set.seed(123)
#' # foo.data = list()
#' # for(i in 1:20){
#' #    foo = circleUnif(100)
#' #    foo.data[[i]] = ripsDiag(foo, 1,1)$diagram}
#' #    for(i in 21:40){   
#' #     foo = cbind(runif(100), runif(100))
#' #     foo.data[[i]] = ripsDiag(foo, 1,1)$diagram
#' #     }
#' # GSWkernel = gaus.kernel(foo.data, h =1, dimension = 1,  q = 2)
#' # GGKclass = krein.svm(kernelmat = GSWkernel, y = rep(c(1,2), c(20,20)))
#' @return An object of class \code{krein.svm} containing the fitted model, including: 
#' \describe{
#'   \item{\code{SV}}{a matrix containing the Support Vectors}
#'   \item{\code{index}}{index of the resulting support vectors in the data matrix}
#'   \item{\code{coefs}}{a matrix containing corresponding coefficients times the training labels}
#'   \item{\code{rho}}{value of the (negative) intercept}
#' }
#' @references 
#' \insertRef{loosli2015learning}{kernelTDA}
#' 
#' \insertRef{chang2011libsvm}{kernelTDA}
#' 
#' \insertRef{dimitriadou2008misc}{kernelTDA}
#' @export
krein.svm.default <-
function (kernelmat   = NULL,
          y           = NULL,
          cost        = 1,
          class.weights = NULL,
          cross       = 0,
          probability = FALSE,
          fitted      = TRUE,
          subset,
          ...)
{
    
    scale = FALSE
    sparse = FALSE
    shrinking   = FALSE
    na.action = na.omit
    
    degree      = 0
    gamma       = 1
    coef0       = 0
    nu          = 0.5
    cachesize   = 40
    tolerance   = 0.001
    epsilon     = 0.1
    
    
    
    x = kernelmat
    if(!is.factor(y)) y = as.factor(y)
    yorig <- y
    
    if(inherits(x, "Matrix")) {
        loadNamespace("SparseM")
        loadNamespace("Matrix")
        x <- as(x, "matrix.csr")
    }
    if (sparse <- inherits(x, "matrix.csr"))
        loadNamespace("SparseM")

    ## NULL parameters?
    if(is.null(cost)) stop(sQuote("cost"), " must not be NULL!")
    if(is.null(tolerance)) stop(sQuote("tolerance"), " must not be NULL!")

    x.scale <- y.scale <- NULL

    type = 0
    ## determine model type
    # if (is.null(type)) type <-
    #     if (is.null(y)) "one-classification"
    #     else if (is.factor(y)) "C-classification"
    #     else "eps-regression"

    # for now only C-classification is implemented for indefinite kernels
    
    # type <- pmatch(type, c("C-classification",
    #                        "nu-classification",
    #                        "one-classification",
    #                        "eps-regression",
    #                        "nu-regression"), 99) - 1
    # if (type > 10) stop("wrong type specification!")
    
    
    kernel = 4

    # kernel <- pmatch(kernel, c("linear",
    #                            "polynomial",
    #                            "radial",
    #                            "sigmoid",
    #                            "precomputed"), 99) - 1
    # 
    # if (kernel > 10) stop("wrong kernel specification!")

#    if (nrow(x) != ncol(x)) stop("kernel matrix must be squared")
    if(kernel == 4) x = cbind(1:nrow(x), x)
    xhold   <- if (fitted) x else NULL
    
    nac <- attr(x, "na.action")

    ## further parameter checks
    nr <- nrow(x)
    if (cross > nr)
        stop(sQuote("cross"), " cannot exceed the number of observations!")

    ytmp <- y
    attributes(ytmp) <- NULL
    if (!is.vector(ytmp) && !is.factor(y) && type != 2)
        stop("y must be a vector or a factor.")
    if (type != 2 && length(y) != nr)
        stop("x and y don't match.")

    if (cachesize < 0.1)
        cachesize <- 0.1

    if (type > 2 && !is.numeric(y))
        stop("Need numeric dependent variable for regression.")
    
    if(is.factor(y)){
        lev = levels(y)
        class.weights = rep(1, length(lev))
        names(class.weights) <- lev
        weightlabels <- match (names(class.weights), lev)
        }
    if(is.numeric(y)){
        lev <- 1:nr
        weightlabels <- 1:nr
    }

    
    ## in case of classification: transform factors into integers
    if (type == 2) # one class classification --> set dummy
        y <- rep(1, nr)
    else
        if (is.factor(y)) {
            lev <- levels(y)
            y <- as.integer(y)
        } else {
            if (type < 3) {
                if(any(as.integer(y) != y))
                    stop("dependent variable has to be of factor or integer type for classification mode.")
                y <- as.factor(y)
                lev <- levels(y)
                y <- as.integer(y)
            } else lev <- unique(y)
        }

    if (type < 3 && !is.null(class.weights)) {
        if (is.character(class.weights) && class.weights == "inverse")
            class.weights <- 1 / table(y)
        if (is.null(names(class.weights)))
            stop("Weights have to be specified along with their according level names !")
        weightlabels <- match (names(class.weights), lev)
        if (any(is.na(weightlabels)))
            stop("At least one level name is missing or misspelled.")
    }
    
    if(is.null(class.weights)){
        if(is.factor(y)){
            class.weights <- rep(1, length(levels(y)))
            names(class.weights) <- levels(y)
            weightlabels <- match (names(class.weights), levels(y))
        } else { 
            class.weights <- rep(1, nr)
            names(class.weights) <- y
        }
    }

    
    nclass <- 2
    if (type < 2) nclass <- length(lev)

    
    if (is.null(type)) stop("type argument must not be NULL!")
    if (is.null(kernel)) stop("kernel argument must not be NULL!")
    if (is.null(degree)) stop("degree argument must not be NULL!")
    if (is.null(gamma)) stop("gamma argument must not be NULL!")
    if (is.null(coef0)) stop("coef0 seed argument must not be NULL!")
    if (is.null(cost)) stop("cost argument must not be NULL!")
    if (is.null(nu)) stop("nu argument must not be NULL!")
    if (is.null(cachesize)) stop("cachesize argument must not be NULL!")
    if (is.null(tolerance)) stop("tolerance argument must not be NULL!")
    if (is.null(epsilon)) stop("epsilon argument must not be NULL!")
    if (is.null(shrinking)) stop("shrinking argument must not be NULL!")
    if (is.null(cross)) stop("cross argument must not be NULL!")
    if (is.null(sparse)) stop("sparse argument must not be NULL!")
    if (is.null(probability)) stop("probability argument must not be NULL!")

    cret <- svmtrain_R(
                ## data
                if (sparse) x@ra else t(x),
                as.integer (nr), as.integer(ncol(x)),
                as.double  (y),
                ## sparse index info
                as.integer (if (sparse) x@ia else 0),
                as.integer (if (sparse) x@ja else 0),

                ## parameters
                as.integer (type),
                as.integer (kernel),
                as.integer (degree),
                as.double  (gamma),
                as.double  (coef0),
                as.double  (cost),
                as.double  (nu),
                as.integer (weightlabels),
                as.double  (class.weights),
                as.integer (length (class.weights)),
                as.double  (cachesize),
                as.double  (tolerance),
                as.double  (epsilon),
                as.integer (shrinking),
                as.integer (cross),
                as.integer (sparse),
                as.integer (probability),

                ## results
                nclasses = integer  (1),
                nr       = integer  (1), # nr of support vectors
                index    = integer  (nr),
                labels   = integer  (nclass),
                nSV      = integer  (nclass),
                rho      = double   (nclass * (nclass - 1) / 2),
                coefs    = double   (nr * (nclass - 1)),
                sigma    = double   (1),
                probA    = double   (nclass * (nclass - 1) / 2),
                probB    = double   (nclass * (nclass - 1) / 2),

                cresults = double   (cross),
                ctotal1  = double   (1),
                ctotal2  = double   (1)

                )

    if (cret$error != "")
        stop(paste(cret$error, "!", sep=""))

    cret$index <- cret$index[1:cret$nr]

    ret <- list (
                 call     = match.call(),
                 type     = type,
                 kernel   = kernel,
                 cost     = cost,
                 degree   = degree,
                 gamma    = gamma,
                 coef0    = coef0,
                 nu       = nu,
                 epsilon  = epsilon,
                 sparse   = sparse,
                 scaled   = scale,
                 x.scale  = x.scale,
                 y.scale  = y.scale,

                 nclasses = cret$nclasses, #number of classes
                 levels   = lev,
                 tot.nSV  = cret$nr, #total number of sv
                 nSV      = cret$nSV[1:cret$nclasses], #number of SV in diff. classes
                 labels   = cret$labels[1:cret$nclasses], #labels of the SVs.
                 SV       = if (sparse) SparseM::t(SparseM::t(x[cret$index]))
                 else t(t(x[cret$index,,drop = FALSE])), #copy of SV
                 index    = cret$index,  #indexes of sv in x
                 ##constants in decision functions
                 rho      = cret$rho[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
                 ##probabilites
                 compprob = probability,
                 probA    = if (!probability) NULL else
                 cret$probA[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
                 probB    = if (!probability) NULL else
                 cret$probB[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
                 sigma    = if (probability) cret$sigma else NULL,
                 ##coefficiants of sv
                 coefs    = if (cret$nr == 0) NULL else
                 t(matrix(cret$coefs[1:((cret$nclasses - 1) * cret$nr)],
                          nrow = cret$nclasses - 1,
                          byrow = TRUE)),
                 na.action = nac
                 )

    ## cross-validation-results
    if (cross > 0)
        if (type > 2) {
            scale.factor     <- if (any(scale)) crossprod(y.scale$"scaled:scale") else 1;
            ret$MSE          <- cret$cresults * scale.factor;
            ret$tot.MSE      <- cret$ctotal1  * scale.factor;
            ret$scorrcoeff   <- cret$ctotal2;
        } else {
            ret$accuracies   <- cret$cresults;
            ret$tot.accuracy <- cret$ctotal1;
        }

    class (ret) <- "krein.svm"

    if (fitted) {
        ret$fitted <- na.action(predict(ret, xhold,
                                        decision.values = TRUE))
        ret$decision.values <- attr(ret$fitted, "decision.values")
        attr(ret$fitted, "decision.values") <- NULL
        if (type > 1) ret$residuals <- yorig - ret$fitted
    }

    ret
}

#' @export
predict.krein.svm <-
function (object, newdata,
          decision.values = FALSE,
          probability = FALSE,
          ...,
          na.action = na.omit)
{
    if (missing(newdata))
        return(fitted(object))

    if (object$tot.nSV < 1)
        stop("Model is empty!")


    if(inherits(newdata, "Matrix")) {
        loadNamespace("SparseM")
        loadNamespace("Matrix")
        newdata <- as(newdata, "matrix.csr")
    }
    if(inherits(newdata, "simple_triplet_matrix")) {
       loadNamespace("SparseM")
       ind <- order(newdata$i, newdata$j)
       newdata <- new("matrix.csr",
                      ra = newdata$v[ind],
                      ja = newdata$j[ind],
                      ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))),
                      dimension = c(newdata$nrow, newdata$ncol))
   }

    sparse <- inherits(newdata, "matrix.csr")
    if (object$sparse || sparse)
        loadNamespace("SparseM")

    act <- NULL
    if ((is.vector(newdata) && is.atomic(newdata)))
        newdata <- t(t(newdata))
    if (sparse)
        newdata <- SparseM::t(SparseM::t(newdata))
    preprocessed <- !is.null(attr(newdata, "na.action"))
    rowns <- if (!is.null(rownames(newdata)))
        rownames(newdata)
    else
        1:nrow(newdata)
    if (!object$sparse) {
        if (inherits(object, "krein.svm.formula")) {
            if(is.null(colnames(newdata)))
                colnames(newdata) <- colnames(object$SV)
            newdata <- na.action(newdata)
            act <- attr(newdata, "na.action")
            newdata <- model.matrix(delete.response(terms(object)),
                                    as.data.frame(newdata))
        } else {
            newdata <- na.action(as.matrix(newdata))
            act <- attr(newdata, "na.action")
        }
    }

    if (!is.null(act) && !preprocessed)
        rowns <- rowns[-act]

    if (any(object$scaled))
        newdata[,object$scaled] <-
            scale(newdata[,object$scaled, drop = FALSE],
                  center = object$x.scale$"scaled:center",
                  scale  = object$x.scale$"scaled:scale"
                  )

    if (ncol(object$SV) != ncol(newdata))
        stop ("test data does not match model !")

    ret <- svmpredict_R(
               as.integer (decision.values),
               as.integer (probability),

               ## model
               if (object$sparse) object$SV@ra else t(object$SV),
               as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
               as.integer (if (object$sparse) object$SV@ia else 0),
               as.integer (if (object$sparse) object$SV@ja else 0),
               as.double  (as.vector(object$coefs)),
               as.double  (object$rho),
               as.integer (object$compprob),
               as.double  (if (object$compprob) object$probA else 0),
               as.double  (if (object$compprob) object$probB else 0),
               as.integer (object$nclasses),
               as.integer (object$tot.nSV),
               as.integer (object$labels),
               as.integer (object$nSV),
               as.integer (object$sparse),

               ## parameter
               as.integer (object$type),
               as.integer (object$kernel),
               as.integer (object$degree),
               as.double  (object$gamma),
               as.double  (object$coef0),

               ## test matrix
               if (sparse) newdata@ra else t(newdata),
               as.integer (nrow(newdata)),
               as.integer (if (sparse) newdata@ia else 0),
               as.integer (if (sparse) newdata@ja else 0),
               as.integer (sparse),

               ## decision-values
               ret = double(nrow(newdata)),
               dec = double(nrow(newdata) * object$nclasses * (object$nclasses - 1) / 2),
               prob = double(nrow(newdata) * object$nclasses)


               )

    ret2 <- if (is.character(object$levels)) # classification: return factors
        factor (object$levels[ret$ret], levels = object$levels)
    else if (object$type == 2) # one-class-classification: return TRUE/FALSE
        ret$ret == 1
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
        ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
    else
        ret$ret

    names(ret2) <- rowns
    ret2 <- napredict(act, ret2)

    if (decision.values) {
        colns = c()
        for (i in 1:(object$nclasses - 1))
            for (j in (i + 1):object$nclasses)
                colns <- c(colns,
                           paste(object$levels[object$labels[i]],
                                 "/", object$levels[object$labels[j]],
                                 sep = ""))
        attr(ret2, "decision.values") <-
            napredict(act,
                      matrix(ret$dec, nrow = nrow(newdata), byrow = TRUE,
                             dimnames = list(rowns, colns)
                             )
                      )
    }

    if (probability && object$type < 2) {
        if (!object$compprob)
            warning("SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
        else
            attr(ret2, "probabilities") <-
                napredict(act,
                          matrix(ret$prob, nrow = nrow(newdata), byrow = TRUE,
                                 dimnames = list(rowns, object$levels[object$labels])
                                 )
                          )
    }

    ret2
}

#' @export
print.krein.svm <-
function (x, ...)
{
    cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
    cat("Parameters:\n")
    cat("   SVM-Type: ", c("C-classification",
                           "nu-classification",
                           "one-classification",
                           "eps-regression",
                           "nu-regression")[x$type+1], "\n")
    cat(" SVM-Kernel: ", c("linear",
                           "polynomial",
                           "radial",
                           "sigmoid")[x$kernel+1], "\n")
    if (x$type==0 || x$type==3 || x$type==4)
        cat("       cost: ", x$cost, "\n")
    cat("\nNumber of Support Vectors: ", x$tot.nSV)
    cat("\n\n")

}

#' @export
summary.krein.svm <-
function(object, ...)
    structure(object, class="summary.krein.svm")

#' @export
print.summary.krein.svm <-
function (x, ...)
{
    print.krein.svm(x)
    if (x$type<2) {
        cat(" (", x$nSV, ")\n\n")
        cat("\nNumber of Classes: ", x$nclasses, "\n\n")
        cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels)
    }
    cat("\n\n")
    if (x$type==2) cat("\nNumber of Classes: 1\n\n\n")

    if ("MSE" %in% names(x)) {
        cat(length (x$MSE), "-fold cross-validation on training data:\n\n", sep="")
        cat("Total Mean Squared Error:", x$tot.MSE, "\n")
        cat("Squared Correlation Coefficient:", x$scorrcoef, "\n")
        cat("Mean Squared Errors:\n", x$MSE, "\n\n")
    }
    if ("accuracies" %in% names(x)) {
        cat(length (x$accuracies), "-fold cross-validation on training data:\n\n", sep="")
        cat("Total Accuracy:", x$tot.accuracy, "\n")
        cat("Single Accuracies:\n", x$accuracies, "\n\n")
    }
    cat("\n\n")
}

#' @export
scale.data.frame <-
function(x, center = TRUE, scale = TRUE)
{
    i <- sapply(x, is.numeric)
    if (ncol(x[, i, drop = FALSE])) {
        x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
        if(center || !is.logical(center))
            attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
        if(scale || !is.logical(scale))
            attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
}
