
#' Read simulation results
#'
#' @examples
#' res <- read_results("inst/simResults-test/")
#' ano <- run_anova(data = res)
#' vis_results(pdat = ano, save = FALSE)
#' nadat <- res %>%
#'    group_by(n, EPV, prev, rho, sparsity, model) %>%
#'    summarise(frac_na = round(100 * mean(is.na(brier)), 1))
#' vis_na(nadat, save = FALSE)
#' vis_calibration(res, save = FALSE)
#'
#' res <- read_results("inst/simResults-test/", "coefs")
#' cres <- filter_coef(res, "X.0")
#' sumr <- summarize_coef(cres)
#' vis_coefs(cres, save = FALSE)
#' nadat <- cres %>%
#'    group_by(n, EPV, prev, rho, sparsity, method, coef) %>%
#'    summarise(frac_na = round(100 * mean(is.na(estimate)), 1))
#' vis_coefs_na(nadat, save = FALSE)
#'
#' @export
#'
read_results <- function(path, which = c("estimands", "coefs")) {

  which <- match.arg(which)

  files <- list.files(path, pattern = ".rds", full.names = TRUE)

  simres <- lapply(X = files, FUN = function(filename) {
    resList <- readRDS(file = filename)
    lapply(X = resList$results, FUN = function(resListSimi) {
      resListSimi[[which]]
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()

  if (is.null(simres$sparsity))
    simres <- simres %>% mutate(sparsity = 0)

  if (which == "coefs")
    return(simres)

  adat <- simres %>%
    mutate(inputp = ceiling(n * prev / EPV)) %>%
    filter(inputp != 1) %>%
    mutate_at(c("n", "EPV", "prev", "rho", "sparsity"),
              ~ factor(.x, levels = sort(unique(.x)))) %>%
    mutate(fct = factor(paste0(model, "n", n, "EPV", EPV, "prev",
                               prev, "rho", rho, "sparsity", sparsity)))
}


# ANOVA -------------------------------------------------------------------

#' Analyze simulation results via anova
#'
#' @export
#'
run_anova <- function(
  formula = brier ~ 0 + fct, data,
  conds = with(data, unique(paste0("n", n, "EPV", EPV, "prev", prev,
                                   "rho", rho, "sparsity", sparsity))),
  models = c("GLM", "EN", "AEN", "RF"), compare_against = "AINET"
) {
  m <- lm(formula, data = data)
  out <- list()

  pb <- txtProgressBar(min = 1, max = length(conds), style = 3)
  for (cond in seq_along(conds)) {
    setTxtProgressBar(pb, cond)
    g1 <- paste0("fct", compare_against, conds[cond])
    g2 <- paste0("fct", models, conds[cond])
    lfct <- paste(g1, "-", g2, "== 0")
    res <- try(glht(m, linfct = lfct))
    if (inherits(res, "try-error"))
      next
    pval <- summary(res)$test$pvalues
    cf <- confint(res)$confint
    nms <- str_split(conds[cond], pattern = "[0-9]|\\.")[[1]]
    nms <- nms[nms != ""]
    nums <- str_split(conds[cond], pattern = "[a-z]|[A-Z]")[[1]]
    nums <- nums[nums != ""]
    names(nums) <- nms
    out[[cond]] <- data.frame(cf, pval = pval, contrast = models, t(nums))
  }

  return(out)
}

#' Visualize anova results
#'
#' @export
#'
vis_results <- function(
  pdat, xlab = "brier", save = TRUE,
  lim = range(pdat$Estimate, na.rm = TRUE)
) {
  xxlab <- switch(xlab, "brier" = "Brier score",
                  "scaledBrier" = "scaled Brier score",
                  "nll" = "log score",
                  "acc" = "accuracy",
                  "auc" = "AUC")

  out2 <- pdat %>%
    bind_rows() %>%
    mutate_at(c("n", "EPV", "prev", "rho", "sparsity"),
              ~ factor(.x, levels = sort(unique(as.numeric(as.character(.x))))))

  rho_plot <- function(trho, tsparse) {
    ggplot(out2 %>% filter(rho == trho, sparsity == tsparse),
           aes(x = contrast, y = Estimate, ymin = lwr, ymax = upr,
               color = ordered(EPV))) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
      geom_pointrange(fatten = 0.75, position = position_dodge(width = 0.7)) +
      geom_errorbar(width = 0.35, position = position_dodge(width = 0.7)) +
      facet_grid(prev ~ n, labeller = label_both) +
      theme_bw() +
      theme(legend.position = "top", panel.grid.major.y = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 7)) +
      labs(y = paste("Difference in", xxlab, "(AINET - other method)"),
           x = "Contrast", subtitle = bquote(rho==~.(trho)~sparsity==.(tsparse)),
           color = "EPV") +
      geom_vline(xintercept = seq(1.5, 3.5, 1), alpha = 0.1, size = 0.8) +
      coord_flip()
  }

  lapply(unique(as.numeric(as.character(out2$sparsity))), function(sparse) {

    ps <- lapply(unique(as.numeric(as.character(out2$rho))), function(rho) {
      rho_plot(trho = rho, tsparse = sparse)
    })

    pfs <- ggarrange(plotlist = ps, common.legend = TRUE, ncol = 2, nrow = 2)

    if (save) {
      pnm <- file.path(outdir, paste0("tie-fighter_", xlab, "_sparsity", sparse, ".pdf"))
      ggsave(pnm, plot = pfs, height = 1.5 * 8.3, width = 1.5 * 11.7)
      write.csv(out2, file.path(outdir, paste0("anova_", xlab, ".csv")),
                row.names = FALSE, quote = FALSE)
    }

    pfs
  })

}

#' Visualize missing estimates
#'
#' @export
#'
vis_na <- function(
  pdat, rhos = unique(as.numeric(as.character(pdat$rho))),
  sparsities = unique(as.numeric(as.character(pdat$sparsity))),
  xlab = "brier", save = TRUE
) {
  lapply(sparsities, function(tsparse) {
    lop <- lapply(rhos, function(trho) {
      pdat %>%
        filter(rho == trho, sparsity == tsparse) %>%
        ggplot(aes(x = model, y = EPV, fill = frac_na)) +
        geom_tile(alpha = 0.3) +
        geom_text(aes(label = paste0(frac_na, "%"))) +
        facet_grid(prev ~ n, labeller = label_both) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(y = "EPV",
             x = element_blank(),
             subtitle = bquote(rho==~.(trho)~sparsity==.(tsparse)),
             fill = paste0("Percent missing in ", xlab)) +
        coord_flip() +
        scale_fill_viridis_c()
    })
    pf <- ggarrange(plotlist = lop, common.legend = TRUE, ncol = 2, nrow = 2)
    if (save) {
      pnm <- file.path(outdir, paste0("missing_", xlab, "_sparsity", tsparse, ".pdf"))
      ggsave(pnm, plot = pf, height = 1.5 * 8.3, width = 1.5 * 11.7)
    }
    pf
  })
}


# Coefs -------------------------------------------------------------------

#' Filter coefficients from simulation results
#'
#' @export
#'
filter_coef <- function(dat, cf = "X.0") {
    dat %>%
        filter(coef == cf) %>%
        gather("method", "estimate", AINET:AEN) %>%
        group_by(ID, n, EPV, prev, sigma2, rho, sparsity, p, coef, method) %>%
        mutate(bias = estimate - oracle)
}

#' Summarize coefficients from simulation results
#'
#' @export
#'
summarize_coef <- function(dat, svars = c("estimate", "bias"),
                           funs = list(mean = mean, sd = sd,
                                       n = function(.x, ...) length(na.omit(.x)))) {
    dat %>%
        summarize_at(svars, funs, na.rm = TRUE) %>%
        mutate(bias_lwr = bias_mean - 1.96 * bias_sd / sqrt(bias_n),
               bias_upr = bias_mean + 1.96 * bias_sd / sqrt(bias_n))
}

#' Viz coefficients from simulation results
#'
#' @export
#'
vis_coefs <- function(pdat, cf = "X.0", save = TRUE,
                        lim = range(pdat$estimate, na.rm = TRUE)) {

    out <- pdat %>%
        mutate_at(c("n", "EPV", "prev", "rho", "sparsity"),
                  ~ factor(.x, levels = sort(unique(as.numeric(as.character(.x))))))

    rho_plot <- function(trho, tsparse) {
        ggplot(out %>% filter(rho == trho, sparsity == tsparse),
               aes(x = method, y = bias,
                   color = ordered(EPV))) +
            geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
            geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA) +
            stat_mean(shape = 4, position = position_dodge(width = 1)) +
            facet_grid(prev ~ n, labeller = label_both) +
            theme_bw() +
            theme(legend.position = "top", panel.grid.major.y = element_blank(),
                  axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 7)) +
            labs(y = "Bias", x = "Method", subtitle = bquote(rho==~.(trho)),
                 color = "EPV", caption = cf) +
            geom_vline(xintercept = seq(1.5, 3.5, 1), alpha = 0.1, size = 0.8) +
            coord_flip(ylim = lim)
    }

    lapply(unique(as.numeric(as.character(out$sparsity))), function(sparse) {
      ps <- lapply(unique(as.numeric(as.character(out$rho))), function(rho) {
            rho_plot(trho = rho, tsparse = sparse)
        })

        pf <- ggarrange(plotlist = ps, common.legend = TRUE, ncol = 2, nrow = 2)

        if (save) {
            pnm <- file.path(outdir, paste0("coef", cf, "_sparsity", sparse, ".pdf"))
            ggsave(pnm, plot = pf, height = 1.5 * 8.3, width = 1.5 * 11.7)
        }

        pf
    })

}

#' Viz missing coefficients from simulation results
#'
#' @export
#'
vis_coefs_na <- function(pdat, rhos = unique(as.numeric(as.character(pdat$rho))),
                   sparsities = unique(as.numeric(as.character(pdat$sparsity))),
                   cf = "X.0", save = TRUE) {
  lapply(sparsities, function(tsparse) {
    lop <- lapply(rhos, function(trho) {
            pdat %>%
                filter(rho == trho, sparsity == tsparse) %>%
                ggplot(aes(x = method, y = ordered(EPV), fill = frac_na)) +
                geom_tile(alpha = 0.3) +
                geom_text(aes(label = paste0(frac_na, "%"))) +
                facet_grid(prev ~ n, labeller = label_both) +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                labs(y = "EPV",
                     x = element_blank(),
                     subtitle = bquote(rho==~.(trho)~sparsity==.(tsparse)),
                     fill = paste0("Percent missing in ", cf)) +
                coord_flip() +
                scale_fill_viridis_c()
        })
        pf <- ggarrange(plotlist = lop, common.legend = TRUE, ncol = 2, nrow = 2)
        if (save) {
          pnm <- file.path(outdir, paste0("coef-missing_", cf, "_sparisty", tsparse, ".pdf"))
          ggsave(pnm, plot = pf, height = 1.5 * 8.3, width = 1.5 * 11.7)
        }
        pf
    })
}

# Calibration -------------------------------------------------------------

#' Viz calibration from simulation results
#'
#' @export
#'
vis_calibration <- function(pdat, metric = c("cslope", "clarge"), save = TRUE,
                            lim = c(-100, 100), only_one = FALSE) {

	metric <- match.arg(metric)
	yint <- switch(metric, "cslope" = 1, "clarge" = 0)
	xxlab <- switch(metric, "cslope" = "calibration slope",
									"clarge" = "calibration in the large")

	out2 <- pdat %>%
		bind_rows() %>%
		mutate_at(c("n", "EPV", "prev", "rho", "sparsity"),
							~ factor(.x, levels = sort(unique(as.numeric(as.character(.x))))))

	nadat <- out2 %>%
		group_by(n, EPV, prev, rho, sparsity, model) %>%
		summarise(frac_na = round(100 * mean(is.na(!!sym(metric))), 1),
							frac_na = paste0(frac_na, "%"))

	rho_plot <- function(trho, tsparse) {
		ggplot(out2 %>% filter(rho == trho, sparsity == tsparse),
					 aes(x = model, y = !!sym(metric), color = ordered(EPV))) +
			geom_hline(yintercept = yint, linetype = 2, alpha = 0.5) +
			geom_boxplot(position = position_dodge(width = 0.7), outlier.size = 0.1) +
			stat_mean(shape = 4, position = position_dodge(width = 0.7)) +
			facet_grid(prev ~ n, labeller = label_both) +
			geom_text(aes(y = lim[1] * 0.8, label = frac_na),
			          data = nadat %>% filter(rho == trho, sparsity == tsparse),
								position = position_dodge(width = 0.7)) +
			theme_bw() +
			theme(legend.position = "top", panel.grid.major.y = element_blank(),
						axis.text.x = element_text(size = 7)) +
			labs(y = xxlab, x = element_blank(),
					 subtitle = bquote(rho==~.(trho)~sparsity==~.(tsparse)), color = "EPV") +
			coord_flip(ylim = lim)
	}

	if (only_one)
		return(rho_plot(0.95, 0.9))


	lapply(unique(as.numeric(as.character(out2$sparsity))), function(sparse) {
	  ps <- lapply(unique(as.numeric(as.character(out2$rho))), function(rho) {
			rho_plot(rho, sparse)
		})
		pf <- ggarrange(plotlist = ps, common.legend = TRUE, ncol = 2, nrow = 2)

		if (save) {
			pnm <- file.path(outdir, paste0("calibration-", metric, "_sparsity", sparse, ".pdf"))
			ggsave(pnm, plot = pf, height = 1.5 * 8.3, width = 1.5 * 11.7)
		}

		pf
	})

}
