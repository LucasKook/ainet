### ainet

#' ainet main function
#'
#' @param formula formula; model formula
#' @param data data.frame; data
#' @param imp_data data.frame; optional data to compute importance out-of-sample
#' @param pen.f function; penalty function (default is \code{NULL} and
#'     \code{ainet:::.vimp} is used)
#' @param plot logical; plot variance importance of \code{rf}
#' @param ... additional arguments to \code{fglmnet()}
#'
#' @return object of class \code{"glmnet"}
#'
#' @export
ainet <- function(formula, data, imp_data = NULL, pen.f = NULL, plot = FALSE,
                  ...) {
  if (is.null(pen.f))
    pen.f <- .vimp(formula, ifelse(is.null(imp_data), list(data),
                                   list(imp_data))[[1]])
  fglmnet(formula, data, penalty.factor = pen.f, ... = ...)
}

#' glmnet with formula interface
#'
#' @param formula formula; model formula
#' @param data data.frame; data
#' @param ... additional args passed to \code{glmnet()}
#'
#' @return object of class \code{"glmnet"}
#'
#' @export
fglmnet <- function(formula, data, ...) {
  dat <- .preproc(formula, data)
  glmnet(x = dat$X, y = dat$Y, ... = ...)
}

#' cv.glmnet with formula interface
#' @inheritParams fglmnet
#' @inheritParams ainet
#' @param nfolds integer; number of folds
#' @param ... additional args to \code{cv.glmnet()}
#'
#' @return object of class \code{"cv.glmnet"}
#'
#' @export
cv.fglmnet <- function(formula, data, imp_data = NULL, pen.f = NULL, nfolds = 5,
                       ...) {
  dat <- .preproc(formula, data)
  if (is.null(pen.f))
    pen.f <- .vimp(formula, ifelse(is.null(imp_data), list(data),
                                   list(imp_data))[[1]])
  cv.glmnet(x = dat$X, y = dat$Y, penalty.factor = pen.f, nfolds = nfolds,
            ... = ...)
}

### Helpers

# Remove intercept
.rm_int <- function(X) {
  if (all(X[, 1] == 1))
    return(X[, -1, drop = FALSE])
  return(X)
}

# preprocess data for use in glmnet
.preproc <- function(formula, data) {
  X <- .rm_int(model.matrix(formula, data))
  Y <- model.response(model.frame(formula, data))
  if (is.factor(Y) & length(levels(Y)) == 2L) {
    Y <- as.numeric(Y) - 1
  }
  return(list(Y = Y, X = X))
}

# Compute variable importance
.vimp <- function(formula, data, which = c("impurity", "adaptive lasso"), ...) {
  which <- match.arg(which)
  if (which == "impurity") {
    rf <- ranger(formula, data, importance = which)
  } else if (which == "adaptive lasso") {
    rf <- glm(formula, data, family = binomial)
  }
  .importance_penalty(rf, which = which, ... = ...)
}

# compute importance penalty
.importance_penalty <- function(rf, gamma = 1, which = c("impurity", "adaptive lasso"),
                                renorm = c("shift", "trunc")) {
  which <- match.arg(which)
  renorm <- match.arg(renorm)
  if (which == "impurity") {
    imp <- importance(rf)
    imp <- switch(renorm, "trunc" = pmax(imp, 0),
                  "shift" = imp - min(imp))
    ret <- 1 - (imp / sum(imp))^gamma
  } else if (which == "adaptive lasso") {
    ret <- 1 / abs(coef(rf))^gamma
  }
  return(ret)
}
