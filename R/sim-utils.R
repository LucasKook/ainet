
### Data generating process

#' Generate data from logistic model
#' @examples
#' dat <- generateData(p = 5, prev = 0.01, rho = 0.5)
#' cor(dat[,-1])
#' @export
generateData <- function(n = 1e2, p = 20, b = rnorm(p), prev = 0.5, rho = 0,
                         nonlin = FALSE) {
  if (is.null(nonlin))
    nonlin <- FALSE
  b0 <- qlogis(p = prev)
  Sigma <- matrix(data = rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
	X <- tX <- rmvn(n = n, mu = rep(0, p), sigma = Sigma)
	if (nonlin)
	  tX[, 1] <- ifelse(X[, 1] > 0, X[, 1], 0)
	Y <- factor(x = rbinom(n = n, size = 1, prob = plogis(q = b0 + tX %*% b)),
              levels = c(0, 1))
	return(data.frame(Y = Y, X = X))
}

### Evaluation

#' Compute accuracy
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' acc(ndat$Y, preds)
#' @export
acc <- function(y_true, y_pred) {
  y_true <- .assert_numeric_binary(y_true)
  y_pred <- y_pred >= 0.5
	mean(y_true == y_pred, na.rm = TRUE)
}

#' Compute negative log-likelihood / n (log-score)
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' nll(ndat$Y, preds)
#' @export
nll <- function(y_true, y_pred) {
  y_true <- .assert_numeric_binary(y_true)
	- mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred), na.rm = TRUE)
}

#' Compute AUC
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' auroc(ndat$Y, preds)
#' @export
auroc <- function(y_true, y_pred) {
  suppressMessages({
    tauc <- auc(y_true, as.vector(y_pred))
    attributes(tauc) <- NULL
    tauc
  })
}

#' Compute Brier score
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' brier(ndat$Y, preds)
#' @export
brier <- function(y_true, y_pred) {
  y_true <- .assert_numeric_binary(y_true)
  mean((y_true - y_pred)^2, na.rm = TRUE)
}

#' Compute variance of Brier score contributions
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' brierVar(ndat$Y, preds)
#' @export
brierVar <- function(y_true, y_pred) {
  y_true <- .assert_numeric_binary(y_true)
  var((y_true - y_pred)^2, na.rm = TRUE)
}

#' Compute scaled Brier score
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' scaledBrier(ndat$Y, preds)
#' @export
scaledBrier <- function(y_true, y_pred) {
  y_true <- .assert_numeric_binary(y_true)
  prev <- mean(y_true, na.rm = TRUE)
  1 - brier(y_true, y_pred) / brier(y_true, prev)
}

#' Compute calibration slope
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' calibrationSlope(ndat$Y, preds)
#' @export
calibrationSlope <- function(y_true, y_pred) {
  logits <- qlogis(p = y_pred)
  m <- try(glm(y_true ~ 1 + logits, family = binomial), silent = TRUE)
  ifelse(inherits(m, "try-error"), NA, unname(coef(m)[2]))
}

#' Compute calibration in the large
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- glm(Y ~ ., data = dat, family = binomial)
#' preds <- predict(m, type = "response", newdata = ndat)
#' calibrationInTheLarge(ndat$Y, preds)
#' @export
calibrationInTheLarge <- function(y_true, y_pred) {
  logits <- qlogis(p = y_pred)
  m <- try(glm(y_true ~ 1 + offset(logits), family = binomial), silent = TRUE)
  ifelse(inherits(m, "try-error"), NA, unname(coef(m)[1]))
}

#' Evaluate model at new observations with `loss`
#' @export
evaluateModel <- function(m, newx, y_true, loss, ...) {
  UseMethod("evaluateModel")
}

#' Evaluate ranger
#' @examples
#' dat <- generateData()
#' ndat <- generateData()
#' m <- ranger(Y ~ ., data = dat, probability = TRUE)
#' evaluateModel(m, newx = ndat, y_true = ndat$Y, loss = nll)
#' @method evaluateModel ranger
#' @export
evaluateModel.ranger <- function(m, newx, y_true, loss, ...) {
	y_pred <- predict(m, data = newx, ... = ...)$predictions[, 2]
	loss(y_true, y_pred)
}

#' Evaluate glmnet
#' @examples
#' set.seed(10)
#' dat <- generateData()
#' ndat <- generateData()
#' m <- fglmnet(Y ~ ., data = dat, alpha = 0.5, lambda = 0.01, family = "binomial")
#' evaluateModel(m, newx = as.matrix(ndat[,-1]), y_true = ndat$Y, loss = nll)
#' @method evaluateModel glmnet
#' @export
evaluateModel.glmnet <- function(m, newx, y_true, loss, ...) {
	y_pred <- predict(m, newx = newx, type = "response", ... = ...)
	loss(y_true, y_pred)
}

### Sim Design Functions

#' SimDesign function for generating the data
#' @examples
#' condition <- data.frame(n = 100, EPV = 20, sparsity = 0.9, prev = 0.01,
#' sigma2 = 1, rho = 0.9, p = 2, nonlin = TRUE)
#' generate(condition)
#' @export
generate <- function(condition, fixed_objects = list(ntest = 1e4)) {
  ## Condition
  n <- condition$n
  epv <- condition$epv
  sigma2 <- condition$sigma2
  p <- condition$p
  rho <- condition$rho
  prev <- condition$prev
  sparsity <- condition$sparsity
  nonlin <- condition$nonlin
  fixed <- ifelse(is.null(condition$fixed), FALSE, condition$fixed)

  ## Simulation of coefficients
  if (fixed)
    betas <- c(1.1, rnorm(p - 1))
  else
    betas <- rnorm(p)

  if (!is.null(sparsity)) {
    idx <- 1 + sample.int(p - 1, ceiling(sparsity * (p - 1)))
    betas[idx] <- 0
  }

  ## Simulate training and test data
  train <- generateData(n = n, p = p, b = betas, prev = prev, rho = rho, nonlin = nonlin)
  test <- generateData(n = fixed_objects$ntest, p = p, b = betas, prev = prev, rho = rho,
                       nonlin = nonlin)
  list(train = train, test = test, beta = betas)
}

#' SimDesign function for analyzing simulated data
#' @examples
#' condition <- data.frame(n = 100, epv = 0.1, sigma2 = 1, p = 2, rho = 0.9,
#' prev = 0.1)
#' dat <- generate(condition)
#' analyze(condition, dat)
#' @export
analyze <- function(condition, dat, fixed_objects = list(ntest = 1e4)) {
  ## Data
  train <- dat$train
  test <- dat$test
  fml <- Y ~ .
  newx <- ainet:::.rm_int(model.matrix(fml, test))
  y_true <- test$Y
  ncoef <- length(dat$beta)
  nmcoef <- c("X.0", paste0("X.", seq_len(ncoef)))
  ocoef <- c(qlogis(condition$prev), dat$beta)

  ## Checks
  # Only events/non-events skipped
  if (all(train$Y == train$Y[1]))
    stop("No events.")

  ## Logistic regression
  if (condition$p < condition$n) { # Fit unpenalized glm if low-dim
    GLM <- try(fglmnet(fml, data = train, alpha = 0, lambda = 0, family = "binomial"),
               silent = TRUE)
  } else { # Tune a ridge regression if high-dim
    cvGLM <- try(cv.fglmnet(fml, data = train, alpha = 0, family = "binomial"),
                 silent = TRUE)
    if (inherits(cvGLM, "try-error")) {
      GLM <- NULL
    } else {
      GLM <- fglmnet(fml, data = train, alpha = 0, lambda = cvGLM$lambda.1se,
                     family = "binomial")
    }
  }
  if (inherits(GLM, "try-error"))
    GLM <- NULL

  ## AINET
	pen.f <- ainet:::.vimp(fml, train, which = "impurity", gamma = 1, renorm = "trunc")
	cvAINET <- try(cv.fglmnet(fml, train, pen.f = pen.f, family = "binomial", relax = TRUE),
	               silent = TRUE)
	if (inherits(cvAINET, "try-error")) {
	  AINET <- NULL
	} else {
	  AINET <- ainet(fml, data = train, pen.f = pen.f, lambda = cvAINET$relaxed$lambda.1se,
	                 alpha = cvAINET$relaxed$gamma.1se, family = "binomial")
	}

	## Elastic net
	cvEN <- try(cv.fglmnet(fml, data = train, relax = TRUE, family = "binomial"),
	            silent = TRUE)
	if (inherits(cvEN, "try-error")) {
	  EN <- NULL
	} else {
	  EN <- fglmnet(fml, data = train, alpha = cvEN$relaxed$gamma.1se,
	                lambda = cvEN$relaxed$lambda.1se, family = "binomial")
	}

	## Adaptive elastic net
	if (!is.null(GLM)) {
	  ENpenf <- 1 / abs(coef(GLM)[-1])
	  cvAEN <- try(cv.fglmnet(fml, data = train, relax = TRUE, pen.f = ENpenf, family = "binomial"),
	               silent = TRUE)
	  if (inherits(cvAEN, "try-error")) {
	    AEN <- NULL
	  } else {
	    AEN <- fglmnet(fml, data = train, pen.f = ENpenf, alpha = cvAEN$relaxed$gamma.1se,
	                   lambda = cvAEN$relaxed$lambda.1se, family = "binomial")
	  }
	} else {
	  AEN <- NULL
	}

	## Random forest
	RF <- try(ranger(fml, data = train, probability = TRUE), silent = TRUE)
	if (inherits(RF, "try-error")) {
	  RF <- NULL
	}

  ## List estimands and models
  metrics <- list(brier = brier, scaledBrier = scaledBrier, brierVar = brierVar, nll = nll,
                  acc = acc, auc = auroc, cslope = calibrationSlope, clarge = calibrationInTheLarge)
  models <- list(AINET = AINET, GLM = GLM, EN = EN, AEN = AEN, RF = RF)

  ## Coefs of all models but RF
  # TODO: coef(NULL) throws NULL, so return rep(NA, length(coefs))
  coefs <- lapply(models[-length(models)],
                  function(mod) {
                    if (is.null(mod)) rep(NA, ncoef) else as.vector(coef(mod))
                  }
  )
  coefs <- do.call("cbind", coefs)
  coefs <- data.frame(condition, coef = nmcoef, coefs, oracle = ocoef)

  ## Compute oracle versions of the estimands
  oracle_predictions <- plogis(qlogis(condition$prev) + newx %*% dat$beta)
  oracles <- lapply(metrics, function(met) met(y_true, oracle_predictions))
  oracles <- do.call("c", oracles)
  names(oracles) <- paste0(names(metrics), "_oracle")

  ## Return (compute all estimands to all predictions from all models)
  res <- sapply(models, function(mod) {
    sapply(metrics, function(met) {
      ifelse(is.null(mod), NA, evaluateModel(mod, newx = newx, y_true = y_true, loss = met))
    })
  })
  ret <- data.frame(condition, t(res), t(oracles))
  ret$model <- names(models)
  list(estimands = ret, coefs = coefs)
}

#' SimDesign function for summarizing simulation results
#' @examples
#' set.seed(123)
#' condition <- data.frame(n = 100, epv = 10, sigma2 = 1, p = 10, rho = 0,
#' prev = 0.5)
#' dat <- generate(condition)
#' res <- analyze(condition, dat)
#' summarize(condition, res)
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows mutate summarise group_by
#' @importFrom magrittr `%>%`
#' @export
summarize <- function(condition, results, fixed_objects = NULL) {
  if (!is.null(names(results[1]))) { # only one sim, for testing
    estimands <- results$estimands
    coefs <- results$coefs
  } else {
    estimands <- bind_rows(lapply(results, function(x) x[["estimands"]]),
                                  .id = "run")
    coefs <- bind_rows(lapply(results, function(x) x[["coefs"]]),
                              .id = "run")
  }

  ## Exceptions (catch exceptions from analyze() function)
  if (nrow(estimands) == 0)
    return(list(estimands = NULL, coefs = NULL))

  ## Summaries
  sumFUN <- function(x, FUNs = list(mean = mean, median = median, sd = sd, iqr = IQR)) {
    ret <- unlist(lapply(FUNs, function(fun) fun(x, na.rm = TRUE)))
    names(ret) <- names(FUNs)
    data.frame(t(ret))
  }

  ## Estimands
  estimands <- gather(estimands, key = "estimand", value = "value", brier:clarge_oracle) %>%
    separate(estimand, into = c("estimand", "oracle"), sep = "_", fill = "right") %>%
    spread(key = oracle, value = value) %>%
    rename(value = `<NA>`)

  estimands_summary <- estimands %>%
    mutate(oracle_adj = value - oracle) %>%
    group_by(model, estimand) %>%
    summarise(value = sumFUN(value), oracle_adj = sumFUN(oracle_adj), oracle = oracle)

  ## Coefficients
  coefs_summary <- gather(coefs, key = "model", value = "estimate", AINET:AEN) %>%
    mutate(bias = estimate - oracle) %>%
    group_by(model, coef) %>%
    summarise(estimate = sumFUN(estimate), bias = sumFUN(bias))

  ## Return
  list(estimands = estimands_summary, coefs = coefs_summary)
}

### Helpers

.assert_numeric_binary <- function(vec) {
	if (is.factor(vec) & length(levels(vec)) == 2)
		vec <- as.numeric(vec) - 1
  if(any(!vec %in% c(0, 1)) | is.factor(vec))
    stop("Wrong formatting for binary outcome. Need numeric vector of 0's and 1's.")
	vec
}
