### ainet

#' ainet main function
#' @export
ainet <- function(formula, data, imp_data = NULL, pen.f = NULL, plot = FALSE,
									 ...) {
	if (plot)
		varImpPlot(rf)
	if (is.null(pen.f))
		pen.f <- .vimp(formula, ifelse(is.null(imp_data), list(data),
																	 list(imp_data))[[1]])
	fglmnet(formula, data, penalty.factor = pen.f, ... = ...)
}

#' glmnet with formula interface
#' @export
fglmnet <- function(formula, data, ...) {
	dat <- .preproc(formula, data)
	glmnet(x = dat$X, y = dat$Y, ... = ...)
}

#' cv.glmnet with formula interface
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
