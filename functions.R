## ===========================================================================
##  Statistics of multivariate extremes under random censoring
##  Martin Bladt
##
##  functions.R -- estimator, standardizations and simulation models.
##
##  This file defines every function used by simulations.R and application.R.
##  It is organized to follow the paper: the directional reduction and the
##  estimator first, then the two marginal standardizations, then the
##  competitors used for comparison, then the data-generating models, and
##  finally the Monte Carlo driver and the common figure style.
##
##  Notation follows the paper throughout.  Z is the n x d matrix of observed
##  coordinates Z^(j) = X^(j) ^ C^(j), delta the matching matrix of censoring
##  indicators, q a direction in (0, infinity)^d, and k the marginal threshold
##  count.
## ===========================================================================
## ---------------------------------------------------------------------------
## 0.  Package check
## ---------------------------------------------------------------------------

required_packages <- function(extra = character()) {
  packages <- unique(c("copula", "ggplot2", "scales", extra))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}


## -------------------------------------------------------------------------
## Tail models and exact targets
## -------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## 1.  The directional reduction and the estimator
##
##     For a direction q the marks A^(j) = psi_j(Z^(j))/q_j are formed, and
##       W(q)   = min_j A^(j)(q),
##       Delta(q) = max{ delta^(j) : A^(j)(q) = W(q) },
##     so that (W, Delta) is a univariate right-censored observation of
##     T(q) = min_j psi_j(X^(j))/q_j with censoring variable D(q).  The
##     Kaplan-Meier product limit of these pairs, read at v = 1, estimates
##     p_n(q) = P(T(q) > 1), and the estimator is the ratio
##       Fhat(q) = phat_n(q) / phat_n(1).
## ---------------------------------------------------------------------------

row_minimum <- function(x) {
  do.call(pmin, as.data.frame(x))
}

ray_product_limit <- function(y, delta, q, influence = FALSE) {
  n <- nrow(y)
  a <- sweep(y, 2L, q, "/")
  w <- row_minimum(a)
  mark <- as.integer(rowSums((a == w) & (delta == 1L)) > 0L)

  ord <- order(w)
  ws <- w[ord]
  ds <- mark[ord]
  runs <- rle(ws)
  last <- cumsum(runs$lengths)
  first <- last - runs$lengths + 1L
  cumulative_events <- cumsum(ds)
  events <- cumulative_events[last] - c(0, cumulative_events[head(last, -1L)])
  risk <- n - first + 1L
  keep <- runs$values <= 1 & events > 0

  if (!any(keep)) {
    return(list(survival = 1, influence = if (influence) numeric(n) else NULL,
                effective = sum(w >= 1)))
  }
  if (any(risk[keep] <= events[keep])) {
    return(list(survival = 0, influence = if (influence) rep(NA_real_, n) else NULL,
                effective = sum(w >= 1)))
  }

  survival <- prod(1 - events[keep] / risk[keep])
  if (!influence) {
    return(list(survival = survival, influence = NULL, effective = sum(w >= 1)))
  }

  failure_time <- runs$values[keep]
  failure_events <- events[keep]
  failure_risk <- risk[keep]
  base <- n * failure_events / (failure_risk * (failure_risk - failure_events))
  index <- findInterval(w, failure_time, rightmost.closed = TRUE)
  cumulative_base <- cumsum(base)
  infl <- numeric(n)
  infl[index > 0L] <- cumulative_base[index[index > 0L]]

  is_failure <- mark == 1L & w <= 1
  failure_index <- match(w[is_failure], failure_time)
  infl[is_failure] <- infl[is_failure] -
    n / (failure_risk[failure_index] - failure_events[failure_index])

  list(survival = survival, influence = infl, effective = sum(w >= 1))
}


directional_km <- function(y, delta, directions, conf_level = NULL,
                           keep_influence = FALSE) {
  y <- as.matrix(y)
  delta <- as.matrix(delta)
  directions <- as.matrix(directions)
  stopifnot(all(dim(y) == dim(delta)), ncol(directions) == ncol(y),
            all(directions > 0))

  need_influence <- !is.null(conf_level) || keep_influence
  reference <- ray_product_limit(y, delta, rep(1, ncol(y)), need_influence)
  if (!is.finite(reference$survival) || reference$survival <= 0) {
    stop("The product limit at the reference direction is zero.")
  }

  fits <- lapply(seq_len(nrow(directions)), function(i) {
    ray_product_limit(y, delta, directions[i, ], need_influence)
  })
  estimate <- vapply(fits, `[[`, numeric(1), "survival") / reference$survival
  effective <- vapply(fits, `[[`, numeric(1), "effective")
  result <- data.frame(directions, estimate = estimate, effective = effective,
                       check.names = FALSE)

  if (need_influence) {
    influence <- vapply(fits, function(fit) fit$influence - reference$influence,
                        numeric(nrow(y)))
    se_log <- apply(influence, 2L, stats::sd) / sqrt(nrow(y))
    result$se_log <- se_log
    if (!is.null(conf_level)) {
      critical <- stats::qnorm(1 - (1 - conf_level) / 2)
      result$lower <- estimate * exp(-critical * se_log)
      result$upper <- estimate * exp(critical * se_log)
    }
    if (keep_influence) attr(result, "influence") <- influence
  }
  result
}


## ---------------------------------------------------------------------------
## 2.  Marginal standardizations
##
##     Multiplicative:  psi_j(x) = x / u_j,  u_j an intermediate marginal
##                      quantile, estimated by a censored Weissman rule.
##     Kaplan-Meier:    psi_j(x) = (k/n) / Fbar_j(x), with Fbar_j the marginal
##                      Kaplan-Meier survival function.  No marginal model.
## ---------------------------------------------------------------------------

km_fit <- function(z, delta) {
  ord <- order(z)
  z <- z[ord]
  delta <- delta[ord]
  runs <- rle(z)
  last <- cumsum(runs$lengths)
  first <- last - runs$lengths + 1L
  cumulative_events <- cumsum(delta)
  events <- cumulative_events[last] - c(0, cumulative_events[head(last, -1L)])
  risk <- length(z) - first + 1L
  factors <- 1 - events / risk
  survival <- cumprod(factors)
  list(time = runs$values, events = events, risk = risk, survival = survival)
}

km_survival <- function(fit, x) {
  index <- findInterval(x, fit$time, rightmost.closed = TRUE)
  out <- rep(1, length(x))
  positive <- index > 0L
  out[positive] <- fit$survival[index[positive]]
  out
}

censored_weissman <- function(z, delta, k) {
  n <- length(z)
  stopifnot(k >= 2L, k < n)
  ord <- order(z)
  threshold <- z[ord[n - k]]
  upper <- ord[(n - k + 1L):n]
  hill_z <- mean(log(z[upper] / threshold))
  p_uncensored <- mean(delta[upper])
  if (p_uncensored <= 0) stop("No uncensored observations above the threshold.")

  gamma <- hill_z / p_uncensored
  survival <- km_survival(km_fit(z, delta), threshold)
  u <- threshold * (survival / (k / n))^gamma
  list(u = unname(u), gamma = gamma, p = p_uncensored,
       threshold_z = threshold, survival = unname(survival))
}

weissman_thresholds <- function(z, delta, k) {
  z <- as.matrix(z)
  delta <- as.matrix(delta)
  fits <- lapply(seq_len(ncol(z)), function(j) censored_weissman(z[, j], delta[, j], k))
  list(
    threshold = vapply(fits, `[[`, numeric(1), "u"),
    gamma = vapply(fits, `[[`, numeric(1), "gamma"),
    p = vapply(fits, `[[`, numeric(1), "p")
  )
}

standardize_multiplicative <- function(z, threshold) {
  sweep(as.matrix(z), 2L, rep_len(threshold, ncol(z)), "/")
}

standardize_km <- function(z, delta, k) {
  z <- as.matrix(z)
  delta <- as.matrix(delta)
  n <- nrow(z)
  out <- matrix(NA_real_, nrow = n, ncol = ncol(z))
  for (j in seq_len(ncol(z))) {
    survival <- km_survival(km_fit(z[, j], delta[, j]), z[, j])
    out[, j] <- (k / n) / survival
  }
  out
}


## -------------------------------------------------------------------------
## Directional product limit and Greenwood influence calculation
## -------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 3.  Competing estimators used for comparison
## ---------------------------------------------------------------------------

naive_tail <- function(y, directions) {
  y <- as.matrix(y)
  directions <- as.matrix(directions)
  reference <- mean(rowSums(sweep(y, 2L, rep(1, ncol(y)), ">")) == ncol(y))
  if (reference <= 0) return(rep(NA_real_, nrow(directions)))
  apply(directions, 1L, function(q) {
    mean(rowSums(sweep(y, 2L, q, ">")) == ncol(y)) / reference
  })
}

dabrowska_tail <- function(y, delta, directions) {
  required_packages("mhazard")
  y <- as.matrix(y)
  delta <- as.matrix(delta)
  directions <- as.matrix(directions)
  stopifnot(ncol(y) == 2L)
  keep <- y[, 1L] >= 1 & y[, 2L] >= 1
  yy <- y[keep, , drop = FALSE]
  dd <- delta[keep, , drop = FALSE]
  if (nrow(yy) < 5L) return(rep(NA_real_, nrow(directions)))

  apply(directions, 1L, function(q) {
    if (all(abs(q - 1) < 1e-12)) return(1)
    t1 <- sort(unique(yy[dd[, 1L] == 1L & yy[, 1L] <= q[1L], 1L]))
    t2 <- sort(unique(yy[dd[, 2L] == 1L & yy[, 2L] <= q[2L], 2L]))
    fit <- try(mhazard:::jointHazLambda(
      yy[, 1L], yy[, 2L], t1, t2, dd[, 1L], dd[, 2L], 1L
    ), silent = TRUE)
    if (inherits(fit, "try-error")) NA_real_ else
      as.numeric(fit[length(t1) + 1L, length(t2) + 1L, 1L])
  })
}


## -------------------------------------------------------------------------
## Monte Carlo and figure conventions
## -------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 4.  Data-generating models and their exact targets
## ---------------------------------------------------------------------------

q_burr <- function(u, tau, lambda) {
  ((1 - u)^(-1 / lambda) - 1)^(1 / tau)
}

gumbel_tail_copula <- function(x, theta) {
  x <- as.numeric(x)
  d <- length(x)
  out <- 0
  for (m in seq_len(d)) {
    terms <- utils::combn(seq_len(d), m, simplify = FALSE)
    out <- out + (-1)^(m + 1L) * sum(vapply(
      terms, function(a) sum(x[a]^theta)^(1 / theta), numeric(1)
    ))
  }
  out
}

gumbel_truth <- function(directions, theta, gamma = NULL) {
  directions <- as.matrix(directions)
  d <- ncol(directions)
  if (is.null(gamma)) gamma <- rep(1, d)
  gamma <- rep_len(gamma, d)
  denominator <- gumbel_tail_copula(rep(1, d), theta)
  apply(directions, 1L, function(q) {
    gumbel_tail_copula(q^(-1 / gamma), theta) / denominator
  })
}

## Dependence structures used by the event vector.  The Gumbel family is
## asymptotically dependent for theta > 1; the Gaussian family is asymptotically
## independent for every |rho| < 1 and is used to violate R^X(1,1) > 0.
event_copula <- function(n, d, family = "gumbel", parameter = 2) {
  required_packages()
  object <- switch(
    family,
    gumbel = copula::gumbelCopula(parameter, dim = d),
    gaussian = copula::normalCopula(parameter, dim = d, dispstr = "ex"),
    stop("Unknown copula family: ", family)
  )
  copula::rCopula(n, object)
}

censoring_uniforms <- function(n, d, censor_theta = 1) {
  if (censor_theta <= 1) return(matrix(stats::runif(n * d), ncol = d))
  required_packages()
  copula::rCopula(n, copula::gumbelCopula(censor_theta, dim = d))
}

simulate_burr <- function(n, theta, censoring, gamma_x = c(0.2, 0.4),
                          tau_x = c(10, 5), lambda_c = 0.5,
                          censor_theta = 1, family = "gumbel") {
  required_packages()
  d <- length(gamma_x)
  censoring <- rep_len(censoring, d)
  lambda_x <- 1 / (gamma_x * tau_x)
  gamma_c <- (1 - censoring) / censoring * gamma_x
  tau_c <- 1 / (gamma_c * lambda_c)

  ux <- event_copula(n, d, family, theta)
  uc <- censoring_uniforms(n, d, censor_theta)
  x <- vapply(seq_len(d), function(j) q_burr(ux[, j], tau_x[j], lambda_x[j]), numeric(n))
  censor <- vapply(seq_len(d), function(j) q_burr(uc[, j], tau_c[j], lambda_c), numeric(n))
  z <- pmin(x, censor)

  list(x = x, censor = censor, z = z, delta = 1L * (x <= censor),
       gamma_x = gamma_x, tau_x = tau_x, lambda_x = lambda_x,
       gamma_c = gamma_c, tau_c = tau_c)
}

## Exact finite-sample target p_n(q)/p_n(1,1) under the multiplicative
## standardization, for Burr margins and a bivariate copula.  It is defined
## whether or not the tail copula is degenerate, and is what the directional
## product limit estimates for every n.
burr_finite_target <- function(directions, family, parameter, gamma_x, tau_x, n, k) {
  required_packages()
  d <- 2L
  lambda_x <- 1 / (gamma_x * tau_x)
  object <- switch(
    family,
    gumbel = copula::gumbelCopula(parameter, dim = d),
    gaussian = copula::normalCopula(parameter, dim = d, dispstr = "ex"),
    stop("Unknown copula family: ", family)
  )
  threshold <- vapply(seq_len(d), function(j) {
    q_burr(1 - k / n, tau_x[j], lambda_x[j])
  }, numeric(1))
  survival <- function(x, j) (1 + x^tau_x[j])^(-lambda_x[j])
  joint <- function(q) {
    s <- vapply(seq_len(d), function(j) survival(q[j] * threshold[j], j), numeric(1))
    1 - sum(1 - s) + copula::pCopula(1 - s, object)
  }
  directions <- as.matrix(directions)
  reference <- joint(rep(1, d))
  list(
    target = apply(directions, 1L, joint) / reference,
    reference = reference, threshold = threshold
  )
}

simulate_weibull <- function(n, theta, censoring, shape_x = c(0.75, 1.5),
                             scale_x = c(1, 1.4)) {
  required_packages()
  d <- length(shape_x)
  censoring <- rep_len(censoring, d)
  beta <- censoring / (1 - censoring)
  scale_c <- scale_x * beta^(-1 / shape_x)

  ux <- copula::rCopula(n, copula::gumbelCopula(theta, dim = d))
  uc <- matrix(stats::runif(n * d), ncol = d)
  x <- vapply(seq_len(d), function(j) {
    stats::qweibull(ux[, j], shape = shape_x[j], scale = scale_x[j])
  }, numeric(n))
  censor <- vapply(seq_len(d), function(j) {
    stats::qweibull(uc[, j], shape = shape_x[j], scale = scale_c[j])
  }, numeric(n))
  z <- pmin(x, censor)

  list(x = x, censor = censor, z = z, delta = 1L * (x <= censor),
       shape_x = shape_x, scale_x = scale_x, scale_c = scale_c)
}

simulate_pareto <- function(n, d, theta, censoring, gamma_x = 0.5,
                            censor_theta = 1) {
  required_packages()
  gamma_x <- rep_len(gamma_x, d)
  censoring <- rep_len(censoring, d)
  gamma_c <- (1 - censoring) / censoring * gamma_x

  ux <- copula::rCopula(n, copula::gumbelCopula(theta, dim = d))
  uc <- censoring_uniforms(n, d, censor_theta)
  x <- sweep(1 - ux, 2L, -gamma_x, "^")
  censor <- sweep(1 - uc, 2L, -gamma_c, "^")
  z <- pmin(x, censor)

  list(x = x, censor = censor, z = z, delta = 1L * (x <= censor),
       gamma_x = gamma_x, gamma_c = gamma_c)
}


## -------------------------------------------------------------------------
## Marginal product limits and censor-adapted Weissman thresholds
## -------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 5.  Monte Carlo driver and figure style
## ---------------------------------------------------------------------------

mc_lapply <- function(n_rep, fun, seed, n_cores = 1L) {
  ids <- seq_len(n_rep)
  worker <- function(i) {
    set.seed(seed + i)
    fun(i)
  }
  if (.Platform$OS.type != "windows" && n_cores > 1L) {
    parallel::mclapply(ids, worker, mc.cores = min(n_cores, n_rep), mc.preschedule = TRUE)
  } else {
    lapply(ids, worker)
  }
}

## Panel labels: one size and weight for every figure in the paper.
panel_label_size <- 8

simulation_palette <- c(
  proposed = "#0072B2",
  dabrowska = "#D55E00",
  naive = "#6B6B6B",
  oracle = "#009E73",
  fixed = "#0072B2",
  matched = "#CC79A7"
)

theme_simulation <- function(base_size = 9) {
  required_packages()
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.22, colour = "#E4E4E4"),
      panel.border = ggplot2::element_rect(linewidth = 0.45, colour = "#303030"),
      axis.title = ggplot2::element_text(colour = "#202020"),
      axis.text = ggplot2::element_text(colour = "#303030"),
      strip.background = ggplot2::element_rect(fill = "#F1F1F1", colour = "#303030",
                                               linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1,
                                         margin = ggplot2::margin(b = 5)),
      plot.tag = ggplot2::element_text(face = "plain", size = panel_label_size)
    )
}

save_simulation_figure <- function(plot, filename, width = 7.0, height = 4.6) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = width, height = height,
                  device = grDevices::pdf, units = "in", useDingbats = FALSE)
}
