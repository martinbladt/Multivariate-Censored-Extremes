## ===========================================================================
##  Statistics of multivariate extremes under random censoring
##  Martin Bladt
##
##  simulations.R -- every Monte Carlo study reported in the paper and in the
##  Supplementary Material.
##
##  Output (written to figures/ and results/):
##
##    Figure 1   simulation_effective_sample.pdf   paper, Section 4.1
##    Figure 2   simulation_dimension.pdf          paper, Section 4.2
##    Figure S1  simulation_directional.pdf        supplement, Section S2.1
##    Figure S2  simulation_directional_additional.pdf   supplement, S2.2
##    Figure S3  simulation_dimension_coverage.pdf       supplement, S2.3
##    Figure S4  simulation_dimension_additional.pdf     supplement, S2.3
##    Figure S5  simulation_robustness.pdf                supplement, S2.4
##    Table S1   study4_*.csv                             supplement, S2.4
##
##  Every study uses B = 500 replications and the common random numbers of
##  mc_lapply, so the results are reproducible on any platform.  Set the
##  environment variable MC_REPS to a smaller number for a quick pass, and
##  REUSE_RESULTS=1 to re-use cached raw output and only redraw the figures.
##
##  Running time is a few hours at the default settings.
## ===========================================================================

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else {
  getwd()
}
source(file.path(root, "functions.R"))
required_packages(c("mhazard", "cowplot"))

figure_dir <- file.path(root, "figures")
result_dir <- file.path(root, "results")
dir.create(figure_dir, showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("MC_REPS", "500"))
n_cores <- as.integer(Sys.getenv(
  "MC_CORES", as.character(max(1L, min(8L, parallel::detectCores() - 1L)))
))
reuse_results <- identical(Sys.getenv("REUSE_RESULTS", "0"), "1")



## ===========================================================================
## 1.  Effective sample size and regularly varying margins (Figure 1)
## ===========================================================================

theta <- 2
gamma_x <- c(0.2, 0.4)

angle <- seq(-log(2), log(2), length.out = 7L)
directions <- cbind(exp(pmax(angle, 0)), exp(pmax(-angle, 0)))
truth <- gumbel_truth(directions, theta, gamma_x)

scenario <- expand.grid(
  n = c(3000L, 6000L, 12000L, 24000L),
  censoring = c(0.10, 0.25, 0.40),
  target_effective = c(25, 50),
  KEEP.OUT.ATTRS = FALSE
)
tail_coefficient <- gumbel_tail_copula(c(1, 1), theta)
scenario$k <- mapply(function(n, censoring, target) {
  beta <- censoring / (1 - censoring)
  round(((target / tail_coefficient) * n^(2 * beta))^(1 / (1 + 2 * beta)))
}, scenario$n, scenario$censoring, scenario$target_effective)

run_scenario <- function(row, scenario_id) {
  n <- row$n
  k <- row$k
  censoring <- row$censoring

  one_replication <- function(replication) {
    dat <- simulate_burr(n, theta, censoring, gamma_x = gamma_x)
    threshold_fit <- try(weissman_thresholds(dat$z, dat$delta, k), silent = TRUE)
    if (inherits(threshold_fit, "try-error")) return(NULL)
    y <- standardize_multiplicative(dat$z, threshold_fit$threshold)

    proposed <- try(directional_km(y, dat$delta, directions), silent = TRUE)
    dabrowska <- try(dabrowska_tail(y, dat$delta, directions), silent = TRUE)
    naive <- naive_tail(y, directions)
    if (inherits(proposed, "try-error")) return(NULL)

    estimates <- cbind(
      proposed = proposed$estimate,
      dabrowska = if (inherits(dabrowska, "try-error")) NA_real_ else dabrowska,
      naive = naive
    )
    out <- data.frame(
      replication = replication,
      direction = seq_along(angle),
      angle = angle,
      truth = truth,
      effective = proposed$effective[which.min(abs(angle))],
      gamma_1 = threshold_fit$gamma[1L],
      gamma_2 = threshold_fit$gamma[2L]
    )
    do.call(rbind, lapply(colnames(estimates), function(method) {
      transform(out, method = method, estimate = estimates[, method])
    }))
  }

  cat(sprintf(
    "Study 1, scenario %02d/%02d: n=%d, k=%d, censoring=%.0f%%, target m=%d\n",
    scenario_id, nrow(scenario), n, k, 100 * censoring, row$target_effective
  ))
  result <- mc_lapply(n_rep, one_replication, seed = 110000 + 1000 * scenario_id,
                      n_cores = n_cores)
  result <- result[!vapply(result, is.null, logical(1))]
  if (!length(result)) stop("All replications failed in scenario ", scenario_id)
  result <- do.call(rbind, result)
  transform(result, scenario = scenario_id, n = n, k = k,
            censoring = censoring, target_effective = row$target_effective)
}

raw_path <- file.path(result_dir, "study1_raw.rds")
if (reuse_results && file.exists(raw_path)) {
  raw <- readRDS(raw_path)
} else {
  raw <- do.call(rbind, lapply(seq_len(nrow(scenario)), function(i) {
    run_scenario(scenario[i, ], i)
  }))
  saveRDS(raw, raw_path)
}

## Integrated performance over non-reference directions.
nonreference <- abs(raw$angle) > 1e-12
split_key <- interaction(raw$scenario, raw$method, raw$replication, drop = TRUE)
replicate_ise <- tapply((raw$estimate[nonreference] - raw$truth[nonreference])^2,
                        split_key[nonreference], mean, na.rm = TRUE)
key_frame <- unique(raw[c("scenario", "method", "replication", "n", "k",
                          "censoring", "target_effective", "effective")])
key_frame$key <- interaction(key_frame$scenario, key_frame$method,
                             key_frame$replication, drop = TRUE)
key_frame$ise <- unname(replicate_ise[as.character(key_frame$key)])

groups <- split(key_frame, interaction(key_frame$scenario, key_frame$method, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) data.frame(
  scenario = x$scenario[1L], method = x$method[1L], n = x$n[1L], k = x$k[1L],
  censoring = x$censoring[1L], target_effective = x$target_effective[1L],
  effective = mean(x$effective, na.rm = TRUE),
  rmse = sqrt(mean(x$ise, na.rm = TRUE)),
  q10 = sqrt(stats::quantile(x$ise, 0.10, na.rm = TRUE)),
  q90 = sqrt(stats::quantile(x$ise, 0.90, na.rm = TRUE)),
  usable = sum(is.finite(x$ise))
)))
summary$method <- factor(summary$method,
                         levels = c("proposed", "dabrowska", "naive"),
                         labels = c("Directional KM", "Dabrowska", "Naive"))
utils::write.csv(summary, file.path(result_dir, "study1_summary.csv"), row.names = FALSE)

method_colours <- c(
  "Directional KM" = unname(simulation_palette["proposed"]),
  "Dabrowska" = unname(simulation_palette["dabrowska"]),
  "Naive" = unname(simulation_palette["naive"])
)
summary$censor_label <- factor(sprintf("%.0f%%", 100 * summary$censoring),
                               levels = c("10%", "25%", "40%"))

proposed_summary <- subset(summary, method == "Directional KM")
scale_colours <- c("10%" = "#0072B2", "25%" = "#E69F00", "40%" = "#D55E00")
scale_limits <- range(proposed_summary$rmse)
scale_limits <- scale_limits + c(-1, 1) * 0.12 * diff(scale_limits)

make_scale_panel <- function(censoring_value, target_value, first_column) {
  x <- subset(proposed_summary,
              censoring == censoring_value & target_effective == target_value)
  x <- x[order(x$n), ]
  x$position <- seq_len(nrow(x))
  ggplot2::ggplot(x, ggplot2::aes(position, rmse)) +
    ggplot2::geom_line(linewidth = 0.6,
                       colour = unname(scale_colours[as.character(x$censor_label[1L])])) +
    ggplot2::geom_point(size = 1.9,
                        colour = unname(scale_colours[as.character(x$censor_label[1L])])) +
    ggplot2::scale_x_continuous(
      breaks = x$position,
      labels = sprintf("%.1f", sqrt(x$effective)),
      sec.axis = ggplot2::dup_axis(
        breaks = x$position, labels = sprintf("%.1f", sqrt(x$k)),
        name = NULL
      )
    ) +
    ggplot2::coord_cartesian(ylim = scale_limits) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = sprintf("%d%%; nominal count %d", round(100 * censoring_value), target_value)
    ) +
    theme_simulation(base_size = 8) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 8.4, hjust = 0.5),
      axis.text.x.top = ggplot2::element_text(colour = "#4D4D4D"),
      axis.text.y = if (first_column) NULL else ggplot2::element_blank(),
      axis.ticks.y = if (first_column) NULL else ggplot2::element_blank()
    )
}

scale_panels <- list()
for (target_value in c(25, 50)) {
  for (censoring_value in c(0.10, 0.25, 0.40)) {
    scale_panels[[length(scale_panels) + 1L]] <-
      make_scale_panel(censoring_value, target_value,
                       first_column = censoring_value == 0.10)
  }
}
shared_axis_label <- function(label) {
  cowplot::ggdraw() + cowplot::draw_label(label, size = 8.5)
}
p_scales <- cowplot::plot_grid(plotlist = scale_panels, ncol = 3,
                              rel_widths = c(1.09, 1, 1))
p_scales <- cowplot::plot_grid(
  shared_axis_label(expression(sqrt(k))), p_scales,
  shared_axis_label(expression(hat(r)[n])),
  ncol = 1, rel_heights = c(0.06, 1, 0.06)
)
p_scales <- cowplot::plot_grid(NULL, p_scales, ncol = 2, rel_widths = c(0.022, 1))
p_scales <- cowplot::ggdraw(p_scales) +
  cowplot::draw_label("Integrated RMSE", x = 0.011, y = 0.5, angle = 90,
                      size = 8.5, vjust = 0.5)

reference_rmse <- subset(summary, method == "Directional KM", c(scenario, rmse))
names(reference_rmse)[2L] <- "reference_rmse"
relative <- merge(subset(summary, method != "Directional KM"), reference_rmse,
                  by = "scenario")
relative$ratio <- relative$rmse / relative$reference_rmse
relative_groups <- split(
  relative,
  interaction(relative$method, relative$censoring,
              relative$target_effective, drop = TRUE)
)
relative_summary <- do.call(rbind, lapply(relative_groups, function(x) data.frame(
  method = x$method[1L], censoring = x$censoring[1L],
  target_effective = x$target_effective[1L],
  mean = mean(x$ratio), lower = min(x$ratio), upper = max(x$ratio),
  configurations = nrow(x)
)))
utils::write.csv(relative_summary, file.path(result_dir, "study1_relative.csv"),
                 row.names = FALSE)
relative_summary$target_label <- factor(
  relative_summary$target_effective,
  levels = c(25, 50), labels = c("Nominal count 25", "Nominal count 50")
)
p_methods <- ggplot2::ggplot(
  relative_summary,
  ggplot2::aes(100 * censoring, mean, colour = method, fill = method, group = method)
) +
  ggplot2::geom_hline(yintercept = 1, linewidth = 0.4, colour = "#4D4D4D") +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::scale_fill_manual(values = method_colours) +
  ggplot2::scale_x_continuous(breaks = c(10, 25, 40), labels = c("10%", "25%", "40%")) +
  ggplot2::scale_y_continuous(breaks = c(1, 1.5, 2, 2.5)) +
  ggplot2::coord_cartesian(ylim = c(0.9, 2.55)) +
  ggplot2::facet_wrap(~target_label, nrow = 1) +
  ggplot2::labs(x = "Tail censoring", y = "Relative integrated RMSE",
                colour = NULL, fill = NULL,
                title = "Competitor error relative to directional KM") +
  theme_simulation() +
  ggplot2::guides(fill = "none") +
  ggplot2::theme(legend.position = "bottom")

figure <- cowplot::plot_grid(
  p_scales, p_methods,
  labels = c("(a)", "(b)"), ncol = 1, rel_heights = c(1.58, 1),
  label_size = panel_label_size, label_fontface = "plain"
)
save_simulation_figure(figure, file.path(figure_dir, "simulation_effective_sample.pdf"),
                       width = 7.2, height = 5.65)

cat("Study 1 complete.\n")
print(summary, row.names = FALSE, digits = 3)


## ===========================================================================
## 2.  Arbitrary margins and the Gaussian approximation (Figure S1)
## ===========================================================================

n <- 8000L
k <- 800L
theta <- 2
shape_x <- c(0.75, 1.5)
scale_x <- c(1, 1.4)

angle <- seq(-log(2.5), log(2.5), length.out = 25L)
directions <- cbind(exp(pmax(angle, 0)), exp(pmax(-angle, 0)))
truth <- gumbel_truth(directions, theta)

scenarios <- list(
  light = list(censoring = c(0.10, 0.10), label = "Symmetric: 10%, 10%"),
  moderate = list(censoring = c(0.25, 0.25), label = "Symmetric: 25%, 25%"),
  asymmetric = list(censoring = c(0.10, 0.40), label = "Asymmetric: 10%, 40%")
)

run_scenario <- function(name, scenario_id) {
  specification <- scenarios[[name]]

  one_replication <- function(replication) {
    dat <- simulate_weibull(n, theta, specification$censoring,
                            shape_x = shape_x, scale_x = scale_x)
    y <- standardize_km(dat$z, dat$delta, k)
    oracle_survival <- vapply(seq_along(shape_x), function(j) {
      exp(-(dat$z[, j] / scale_x[j])^shape_x[j])
    }, numeric(n))
    y_oracle <- (k / n) / oracle_survival

    fit <- try(directional_km(y, dat$delta, directions, conf_level = 0.95), silent = TRUE)
    oracle <- try(directional_km(y_oracle, dat$delta, directions), silent = TRUE)
    if (inherits(fit, "try-error") || inherits(oracle, "try-error")) return(NULL)

    data.frame(
      replication = replication,
      scenario = name,
      scenario_label = specification$label,
      direction = seq_along(angle),
      angle = angle,
      truth = truth,
      proposed = fit$estimate,
      oracle = oracle$estimate,
      naive = naive_tail(y, directions),
      lower = fit$lower,
      upper = fit$upper,
      covered = fit$lower <= truth & truth <= fit$upper,
      effective = fit$effective,
      effective_reference = fit$effective[which.min(abs(angle))],
      sup_plugin_oracle = max(abs(fit$estimate - oracle$estimate), na.rm = TRUE)
    )
  }

  cat(sprintf("Study 2, scenario %d/%d: %s\n", scenario_id, length(scenarios),
              specification$label))
  result <- mc_lapply(n_rep, one_replication, seed = 220000 + 1000 * scenario_id,
                      n_cores = n_cores)
  result <- result[!vapply(result, is.null, logical(1))]
  if (!length(result)) stop("All replications failed in scenario ", name)
  do.call(rbind, result)
}

raw_path <- file.path(result_dir, "study2_raw.rds")
if (reuse_results && file.exists(raw_path)) {
  raw <- readRDS(raw_path)
} else {
  raw <- do.call(rbind, Map(run_scenario, names(scenarios), seq_along(scenarios)))
  saveRDS(raw, raw_path)
}

curve_groups <- split(raw, interaction(raw$scenario, raw$direction, drop = TRUE))
curves <- do.call(rbind, lapply(curve_groups, function(x) data.frame(
  scenario = x$scenario[1L], scenario_label = x$scenario_label[1L],
  direction = x$direction[1L], angle = x$angle[1L], truth = x$truth[1L],
  proposed = mean(x$proposed, na.rm = TRUE),
  proposed_q10 = stats::quantile(x$proposed, 0.10, na.rm = TRUE),
  proposed_q90 = stats::quantile(x$proposed, 0.90, na.rm = TRUE),
  oracle = mean(x$oracle, na.rm = TRUE),
  naive = mean(x$naive, na.rm = TRUE),
  coverage = mean(x$covered, na.rm = TRUE),
  width = mean(x$upper - x$lower, na.rm = TRUE),
  effective = mean(x$effective, na.rm = TRUE),
  usable = sum(is.finite(x$proposed))
)))
scenario_levels <- vapply(scenarios, `[[`, character(1), "label")
curves$scenario_label <- factor(curves$scenario_label, levels = scenario_levels)
utils::write.csv(curves, file.path(result_dir, "study2_summary.csv"), row.names = FALSE)

line_data <- rbind(
  data.frame(curves[c("scenario_label", "angle")], method = "Directional KM",
             estimate = curves$proposed),
  data.frame(curves[c("scenario_label", "angle")], method = "Naive",
             estimate = curves$naive),
  data.frame(curves[c("scenario_label", "angle")], method = "Truth",
             estimate = curves$truth)
)
line_data$method <- factor(line_data$method,
                           levels = c("Truth", "Directional KM", "Naive"))
method_colours <- c(
  "Truth" = "#151515",
  "Directional KM" = unname(simulation_palette["proposed"]),
  "Naive" = unname(simulation_palette["naive"])
)

p_curve <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = curves,
    ggplot2::aes(angle, ymin = proposed_q10, ymax = proposed_q90),
    fill = unname(simulation_palette["proposed"]), alpha = 0.14
  ) +
  ggplot2::geom_line(
    data = line_data,
    ggplot2::aes(angle, estimate, colour = method, linetype = method),
    linewidth = 0.65
  ) +
  ggplot2::facet_wrap(~scenario_label, nrow = 1) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::scale_linetype_manual(values = c("solid", "solid", "22")) +
  ggplot2::scale_x_continuous(
    breaks = log(c(1 / 2, 1, 2)), labels = c("1/2", "1", "2")
  ) +
  ggplot2::labs(x = expression(s / t), y = expression(F^degree(s, t)),
                title = "Directional recovery under Weibull margins") +
  theme_simulation() +
  ggplot2::theme(legend.position = "bottom")

mc_half_width <- 1.96 * sqrt(0.95 * 0.05 / n_rep)
coverage_curves <- subset(curves, abs(angle) > 1e-12)
coverage_curves$side <- factor(sign(coverage_curves$angle))
p_coverage <- ggplot2::ggplot(
  coverage_curves, ggplot2::aes(angle, coverage, group = side)
) +
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                    ymin = 0.95 - mc_half_width, ymax = 0.95 + mc_half_width,
                    fill = "#BDBDBD", alpha = 0.20) +
  ggplot2::geom_hline(yintercept = 0.95, linewidth = 0.45, colour = "#303030") +
  ggplot2::geom_line(colour = unname(simulation_palette["proposed"]), linewidth = 0.65) +
  ggplot2::geom_point(ggplot2::aes(size = effective),
                      colour = unname(simulation_palette["proposed"]), alpha = 0.80) +
  ggplot2::facet_wrap(~scenario_label, nrow = 1) +
  ggplot2::scale_size_continuous(range = c(0.8, 2.3), breaks = NULL) +
  ggplot2::scale_x_continuous(
    breaks = log(c(1 / 2, 1, 2)), labels = c("1/2", "1", "2")
  ) +
  ggplot2::coord_cartesian(ylim = c(0.75, 1.005)) +
  ggplot2::labs(x = expression(s / t), y = "Empirical coverage",
                title = "Pointwise 95% Greenwood intervals (reference omitted)") +
  theme_simulation() +
  ggplot2::theme(legend.position = "none")

figure <- cowplot::plot_grid(
  p_curve, p_coverage, labels = c("(a)", "(b)"), ncol = 1,
  rel_heights = c(1.12, 1), label_size = panel_label_size, label_fontface = "plain"
)
save_simulation_figure(figure, file.path(figure_dir, "simulation_directional.pdf"),
                       width = 7.2, height = 6.1)

discrepancy <- unique(raw[c("replication", "scenario_label", "effective_reference",
                            "sup_plugin_oracle")])
discrepancy$scenario_label <- factor(discrepancy$scenario_label, levels = scenario_levels)
discrepancy$scaled_discrepancy <- sqrt(discrepancy$effective_reference) *
  discrepancy$sup_plugin_oracle
p_additional <- ggplot2::ggplot(
  discrepancy, ggplot2::aes(scenario_label, scaled_discrepancy)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.45,
                       fill = "#D9EAF4", colour = "#303030") +
  ## A fixed jitter seed keeps the figure reproducible.
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.11, height = 0,
                                                          seed = 20260101L),
                      alpha = 0.09, size = 0.65,
                      colour = unname(simulation_palette["proposed"])) +
  ggplot2::labs(x = NULL, y = "Scaled sup-norm discrepancy",
                title = "Plug-in discrepancy on the joint scale") +
  theme_simulation() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1),
                 legend.position = "none")

## -------------------------------------------------------------------------
## Growth in n at a shrinking tail fraction
##
## By (KMrate) and (rate), r_n * zeta_n depends on n only through n/k, so the
## sample fraction has to shrink for the plug-in--oracle separation to vanish.
## -------------------------------------------------------------------------

growth_ns <- c(2000L, 8000L, 32000L, 128000L)
growth_path <- file.path(result_dir, "study2_growth.rds")
if (reuse_results && file.exists(growth_path)) {
  growth <- readRDS(growth_path)
} else {
  ## One cell per (censoring regime, sample size).  Seeds are set per
  ## replication inside mc_lapply, so the result does not depend on the
  ## number of cores.
  growth_design <- expand.grid(n = growth_ns, scenario = names(scenarios),
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  growth <- do.call(rbind, lapply(seq_len(nrow(growth_design)), function(i) {
    nm <- growth_design$scenario[i]
    nn <- growth_design$n[i]
    cc <- scenarios[[nm]]$censoring
    kk <- as.integer(ceiling(nn^0.75))
    D <- unlist(mc_lapply(n_rep, function(b) {
      dat <- simulate_weibull(nn, theta, cc, shape_x = shape_x, scale_x = scale_x)
      y <- standardize_km(dat$z, dat$delta, kk)
      oracle_survival <- vapply(seq_along(shape_x), function(j) {
        exp(-(dat$z[, j] / scale_x[j])^shape_x[j])
      }, numeric(nn))
      fit <- try(directional_km(y, dat$delta, directions), silent = TRUE)
      orc <- try(directional_km((kk / nn) / oracle_survival, dat$delta, directions),
                 silent = TRUE)
      if (inherits(fit, "try-error") || inherits(orc, "try-error")) return(NA_real_)
      sqrt(sum(y[, 1L] >= 1 & y[, 2L] >= 1)) *
        max(abs(fit$estimate - orc$estimate), na.rm = TRUE)
    }, seed = 450000 + 1000L * i, n_cores = n_cores))
    data.frame(scenario_label = scenarios[[nm]]$label, n = nn, k = kk,
               median = stats::median(D, na.rm = TRUE),
               q25 = unname(stats::quantile(D, 0.25, na.rm = TRUE)),
               q75 = unname(stats::quantile(D, 0.75, na.rm = TRUE)),
               q90 = unname(stats::quantile(D, 0.90, na.rm = TRUE)))
  }))
  saveRDS(growth, growth_path)
}
growth$scenario_label <- factor(growth$scenario_label, levels = scenario_levels)
utils::write.csv(growth, file.path(result_dir, "study2_growth.csv"), row.names = FALSE)

p_growth <- ggplot2::ggplot(growth, ggplot2::aes(n, median)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75),
                       fill = unname(simulation_palette["proposed"]), alpha = 0.14) +
  ggplot2::geom_line(colour = unname(simulation_palette["proposed"]), linewidth = 0.65) +
  ggplot2::geom_point(colour = unname(simulation_palette["proposed"]), size = 1.5) +
  ggplot2::facet_wrap(~scenario_label, nrow = 1) +
  ggplot2::scale_x_log10(breaks = growth_ns, labels = c("2", "8", "32", "128")) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(x = "Sample size n (thousands)", y = "Scaled sup-norm discrepancy",
                title = "Error as the sample grows") +
  theme_simulation()

figure_additional <- cowplot::plot_grid(
  p_additional, p_growth, labels = c("(a)", "(b)"), ncol = 1,
  rel_heights = c(1, 1.02), label_size = panel_label_size, label_fontface = "plain"
)
save_simulation_figure(figure_additional,
                       file.path(figure_dir, "simulation_directional_additional.pdf"),
                       width = 7.2, height = 5.9)

cat("Study 2 complete.\n")
print(curves[curves$direction %in% c(5L, 13L, 21L), ], row.names = FALSE, digits = 3)


## ===========================================================================
## 3.  Dimension, joint information and directional paths (Figures 2, S3, S4)
## ===========================================================================

n <- 10000L
dimensions <- 2L:10L
theta <- 3
gamma_x <- 0.5
censoring <- 0.05
fixed_k <- 600L
matched_effective <- 100
maximum_direction <- 2

atlas_grid <- function(d, n_lambda) {
  subset_size <- unique(c(1L, max(1L, floor(d / 2)), d - 1L))
  lambda <- seq(0, 1, length.out = n_lambda)
  out <- lapply(subset_size, function(size) {
    q <- matrix(maximum_direction^lambda, nrow = length(lambda), ncol = d)
    q[, seq_len(size)] <- 1
    data.frame(subset_size = size, lambda = lambda, q, check.names = FALSE)
  })
  do.call(rbind, out)
}

scenario <- expand.grid(d = dimensions, design = c("fixed", "matched"),
                        KEEP.OUT.ATTRS = FALSE)
scenario$k <- mapply(function(d, design) {
  if (design == "fixed") return(fixed_k)
  beta <- censoring / (1 - censoring)
  coefficient <- gumbel_tail_copula(rep(1, d), theta)
  round(((matched_effective / coefficient) * n^(d * beta))^(1 / (1 + d * beta)))
}, scenario$d, scenario$design)

build_atlas_specification <- function(d) {
  atlas <- atlas_grid(d, 21L)
  direction_columns <- paste0("q", seq_len(d))
  names(atlas)[3 + seq_len(d) - 1L] <- direction_columns
  direction_key <- apply(atlas[direction_columns], 1L, paste, collapse = ":")
  keep <- !duplicated(direction_key)
  directions <- as.matrix(atlas[keep, direction_columns])
  map <- match(direction_key, direction_key[keep])
  truth_unique <- gumbel_truth(directions, theta, rep(gamma_x, d))
  truth <- truth_unique[map]
  list(atlas = atlas, directions = directions, map = map, truth = truth)
}
atlas_specification <- setNames(lapply(dimensions, build_atlas_specification),
                                as.character(dimensions))

one_replication <- function(replication) {
  ## The same ten-dimensional latent draw is used throughout the dimension grid.
  ## Taking its first d coordinates gives a nested common-random-number design.
  dat <- simulate_pareto(n, max(dimensions), theta, censoring, gamma_x)
  fixed_threshold_fit <- try(
    weissman_thresholds(dat$z, dat$delta, fixed_k), silent = TRUE
  )
  if (inherits(fixed_threshold_fit, "try-error")) return(NULL)
  fixed_y <- standardize_multiplicative(dat$z, fixed_threshold_fit$threshold)

  result <- vector("list", nrow(scenario))
  for (scenario_id in seq_len(nrow(scenario))) {
    row <- scenario[scenario_id, ]
    d <- row$d
    k <- row$k
    specification <- atlas_specification[[as.character(d)]]
    z <- dat$z[, seq_len(d), drop = FALSE]
    delta <- dat$delta[, seq_len(d), drop = FALSE]

    if (row$design == "fixed") {
      threshold_fit <- list(
        threshold = fixed_threshold_fit$threshold[seq_len(d)],
        gamma = fixed_threshold_fit$gamma[seq_len(d)]
      )
      y <- fixed_y[, seq_len(d), drop = FALSE]
    } else {
      threshold_fit <- try(weissman_thresholds(z, delta, k), silent = TRUE)
      if (inherits(threshold_fit, "try-error")) return(NULL)
      y <- standardize_multiplicative(z, threshold_fit$threshold)
    }
    proposed <- try(
      directional_km(y, delta, specification$directions), silent = TRUE
    )
    if (inherits(proposed, "try-error")) return(NULL)
    estimates <- cbind(
      proposed = proposed$estimate[specification$map],
      naive = naive_tail(y, specification$directions)[specification$map]
    )
    out <- data.frame(
      replication = replication,
      subset_size = specification$atlas$subset_size,
      lambda = specification$atlas$lambda,
      truth = specification$truth,
      effective_reference = proposed$effective[
        which.min(rowSums((specification$directions - 1)^2))
      ],
      mean_gamma = mean(threshold_fit$gamma)
    )
    result[[scenario_id]] <- do.call(rbind, lapply(colnames(estimates), function(method) {
      transform(out, method = method, estimate = estimates[, method])
    }))
    result[[scenario_id]] <- transform(
      result[[scenario_id]], scenario = scenario_id, d = d,
      design = row$design, k = k
    )
  }
  do.call(rbind, result)
}

raw_path <- file.path(result_dir, "study3_raw.rds")
if (reuse_results && file.exists(raw_path)) {
  raw <- readRDS(raw_path)
} else {
  cat(sprintf(
    "Study 3: %d nested replications, dimensions %d--%d, %d configurations\n",
    n_rep, min(dimensions), max(dimensions), nrow(scenario)
  ))
  raw <- mc_lapply(n_rep, one_replication, seed = 330000, n_cores = n_cores)
  raw <- raw[!vapply(raw, is.null, logical(1))]
  if (!length(raw)) stop("All Study 3 replications failed.")
  raw <- do.call(rbind, raw)
  saveRDS(raw, raw_path)
}

nonreference <- raw$lambda > 0
split_key <- interaction(raw$scenario, raw$method, raw$replication, drop = TRUE)
replicate_ise <- tapply((raw$estimate[nonreference] - raw$truth[nonreference])^2,
                        split_key[nonreference], mean, na.rm = TRUE)
key_frame <- unique(raw[c("scenario", "method", "replication", "d", "design", "k",
                          "effective_reference")])
key_frame$key <- interaction(key_frame$scenario, key_frame$method,
                             key_frame$replication, drop = TRUE)
key_frame$ise <- unname(replicate_ise[as.character(key_frame$key)])

groups <- split(key_frame, interaction(key_frame$scenario, key_frame$method, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(x) data.frame(
  scenario = x$scenario[1L], method = x$method[1L], d = x$d[1L],
  design = x$design[1L], k = x$k[1L],
  effective = mean(x$effective_reference, na.rm = TRUE),
  rmse = sqrt(mean(x$ise, na.rm = TRUE)),
  q10 = sqrt(stats::quantile(x$ise, 0.10, na.rm = TRUE)),
  q90 = sqrt(stats::quantile(x$ise, 0.90, na.rm = TRUE)),
  usable = sum(is.finite(x$ise))
)))
summary$method <- factor(summary$method, levels = c("proposed", "naive"),
                         labels = c("Directional KM", "Naive"))
summary$design_label <- factor(summary$design, levels = c("fixed", "matched"),
                               labels = c("Common marginal k", "Matched effective count"))
utils::write.csv(summary, file.path(result_dir, "study3_summary.csv"), row.names = FALSE)

trajectory <- subset(key_frame, method == "proposed" & is.finite(ise))
trajectory$root_ise <- sqrt(trajectory$ise)
trajectory$scaled_error <- sqrt(trajectory$effective_reference) * trajectory$root_ise
trajectory$design_label <- factor(
  trajectory$design, levels = c("fixed", "matched"),
  labels = c("Common marginal k", "Matched effective count")
)
utils::write.csv(trajectory, file.path(result_dir, "study3_trajectories.csv"),
                 row.names = FALSE)

method_colours <- c(
  "Directional KM" = unname(simulation_palette["proposed"]),
  "Naive" = unname(simulation_palette["naive"])
)
p_scale <- ggplot2::ggplot(summary,
                           ggplot2::aes(d, rmse, colour = method, group = method)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q10, ymax = q90, fill = method),
                       colour = NA, alpha = 0.10) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_text(
    data = subset(summary, method == "Directional KM" &
                    ((design == "fixed" & d %% 2L == 0L) |
                     (design == "matched" & d %in% c(2L, 10L)))),
    ggplot2::aes(label = paste0("hat(r)[n]^2==", round(effective))),
    parse = TRUE,
    colour = "#303030", vjust = -0.9, size = 2.25, show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~design_label, nrow = 1) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::scale_fill_manual(values = method_colours) +
  ggplot2::scale_x_continuous(breaks = dimensions,
                              expand = ggplot2::expansion(mult = 0.07)) +
  ggplot2::labs(x = "Dimension d", y = "Integrated RMSE",
                title = "Dimension and effective sample size") +
  theme_simulation() +
  ggplot2::theme(legend.position = "bottom")

atlas_raw <- subset(raw, d == max(dimensions) & design == "matched")
atlas_groups <- split(atlas_raw,
                      interaction(atlas_raw$method, atlas_raw$subset_size,
                                  atlas_raw$lambda, drop = TRUE))
atlas_summary <- do.call(rbind, lapply(atlas_groups, function(x) data.frame(
  method = x$method[1L], subset_size = x$subset_size[1L], lambda = x$lambda[1L],
  truth = x$truth[1L], estimate = mean(x$estimate, na.rm = TRUE),
  q10 = stats::quantile(x$estimate, 0.10, na.rm = TRUE),
  q90 = stats::quantile(x$estimate, 0.90, na.rm = TRUE)
)))
atlas_summary$method <- factor(atlas_summary$method,
                               levels = c("proposed", "naive"),
                               labels = c("Directional KM", "Naive"))
atlas_summary$path <- factor(
  atlas_summary$subset_size,
  levels = c(1L, 5L, 9L),
  labels = c("1 coordinate fixed", "5 coordinates fixed", "9 coordinates fixed")
)
truth_atlas <- unique(atlas_summary[c("subset_size", "lambda", "truth", "path")])

p_atlas <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = subset(atlas_summary, method == "Directional KM"),
    ggplot2::aes(lambda, ymin = q10, ymax = q90),
    fill = unname(simulation_palette["proposed"]), alpha = 0.14
  ) +
  ggplot2::geom_line(
    data = atlas_summary,
    ggplot2::aes(lambda, estimate, colour = method, linetype = method), linewidth = 0.65
  ) +
  ggplot2::geom_line(data = truth_atlas, ggplot2::aes(lambda, truth),
                     colour = "#151515", linewidth = 0.65) +
  ggplot2::facet_wrap(~path, nrow = 1) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::scale_linetype_manual(values = c("solid", "22")) +
  ggplot2::labs(x = expression(lambda), y = expression(F^degree(q[S](lambda))),
                title = "Directional paths for d = 10") +
  theme_simulation() +
  ggplot2::theme(legend.position = "bottom")

figure <- cowplot::plot_grid(
  p_scale, p_atlas, labels = c("(a)", "(b)"), ncol = 1,
  rel_heights = c(1, 1.05), label_size = panel_label_size, label_fontface = "plain"
)
save_simulation_figure(figure, file.path(figure_dir, "simulation_dimension.pdf"),
                       width = 7.2, height = 6.15)

trajectory_mean <- stats::aggregate(
  scaled_error ~ d + design_label, trajectory, mean
)
p_additional <- ggplot2::ggplot(
  trajectory,
  ggplot2::aes(d, scaled_error, group = interaction(replication, design_label))
) +
  ggplot2::geom_line(colour = "#6B6B6B", alpha = 0.035, linewidth = 0.20) +
  ggplot2::geom_line(
    data = trajectory_mean,
    ggplot2::aes(d, scaled_error, colour = design_label, group = design_label),
    linewidth = 0.85, inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = trajectory_mean,
    ggplot2::aes(d, scaled_error, colour = design_label),
    size = 1.8, inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~design_label, nrow = 1) +
  ggplot2::scale_colour_manual(values = c(
    "Common marginal k" = unname(simulation_palette["fixed"]),
    "Matched effective count" = unname(simulation_palette["matched"])
  )) +
  ggplot2::scale_x_continuous(breaks = dimensions) +
  ggplot2::labs(
    x = "Dimension d",
    y = expression(hat(r)[b*","*d] * " replication-level integrated error"),
    title = "Nested replication paths on the effective scale"
  ) +
  theme_simulation() +
  ggplot2::theme(legend.position = "none")
save_simulation_figure(p_additional,
                       file.path(figure_dir, "simulation_dimension_additional.pdf"),
                       width = 7.0, height = 3.35)


## -------------------------------------------------------------------------
## Coverage of the Gaussian intervals across dimension
##
## Study 2 validates (greenwoodCI) at d = 2.  Here the joint count is matched
## across d, so the information and the 21-point grid are held fixed and only
## the dimension varies.
## -------------------------------------------------------------------------

coverage_path <- file.path(result_dir, "study3_coverage.rds")
if (reuse_results && file.exists(coverage_path)) {
  coverage_raw <- readRDS(coverage_path)
} else {
  matched_k <- setNames(scenario$k[scenario$design == "matched"],
                        as.character(scenario$d[scenario$design == "matched"]))
  one_coverage <- function(replication) {
    do.call(rbind, lapply(dimensions, function(d) {
      k <- matched_k[[as.character(d)]]
      specification <- atlas_specification[[as.character(d)]]
      dat <- simulate_pareto(n, d, theta, censoring, gamma_x)
      fit <- try(weissman_thresholds(dat$z, dat$delta, k), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      y <- standardize_multiplicative(dat$z, fit$threshold)
      est <- try(directional_km(y, dat$delta, specification$directions,
                                conf_level = 0.95), silent = TRUE)
      if (inherits(est, "try-error")) return(NULL)
      truth <- specification$truth
      lo <- est$lower[specification$map]
      hi <- est$upper[specification$map]
      eff <- est$effective[specification$map]
      data.frame(
        replication = replication, d = d, k = k,
        subset_size = specification$atlas$subset_size,
        lambda = specification$atlas$lambda,
        covered = as.integer(truth >= lo & truth <= hi),
        width = hi - lo,
        effective = eff
      )
    }))
  }
  cat(sprintf("Study 3 coverage: %d replications, dimensions %d--%d\n",
              n_rep, min(dimensions), max(dimensions)))
  coverage_raw <- mc_lapply(n_rep, one_coverage, seed = 770000, n_cores = n_cores)
  coverage_raw <- do.call(rbind, coverage_raw[!vapply(coverage_raw, is.null, logical(1))])
  saveRDS(coverage_raw, coverage_path)
}

## Pointwise coverage at each direction of the path grid, the reference omitted.
reference_effective <- vapply(
  split(coverage_raw$effective[coverage_raw$lambda == 0],
        coverage_raw$d[coverage_raw$lambda == 0]), mean, numeric(1))
coverage_raw <- coverage_raw[coverage_raw$lambda > 0, ]
coverage_points <- do.call(rbind, lapply(
  split(coverage_raw, list(coverage_raw$d, coverage_raw$subset_size, coverage_raw$lambda),
        drop = TRUE),
  function(x) data.frame(d = x$d[1L], subset_size = x$subset_size[1L],
                         lambda = x$lambda[1L], coverage = mean(x$covered),
                         effective = mean(x$effective))
))
coverage_summary <- do.call(rbind, lapply(split(coverage_raw, coverage_raw$d), function(x) {
  per_rep <- vapply(split(x$covered, x$replication), mean, numeric(1))
  data.frame(d = x$d[1L], k = x$k[1L], replications = length(per_rep),
             coverage = mean(per_rep), se = stats::sd(per_rep) / sqrt(length(per_rep)),
             width = mean(x$width),
             reference_effective = reference_effective[[as.character(x$d[1L])]])
}))
utils::write.csv(coverage_summary, file.path(result_dir, "study3_coverage.csv"),
                 row.names = FALSE)

mc_half_width <- 1.96 * sqrt(0.95 * 0.05 / n_rep)
coverage_points$dimension_label <- factor(paste0("d = ", coverage_points$d),
                                          levels = paste0("d = ", dimensions))
## The three paths of the grid hold 1, floor(d/2) and d-1 coordinates at the
## reference, so they carry different amounts of joint information and must be
## distinguishable rather than drawn on top of one another.
coverage_points$path_role <- with(coverage_points, ifelse(
  subset_size == 1L, "1 coordinate fixed",
  ifelse(subset_size == d - 1L, "d - 1 coordinates fixed",
         "d/2 coordinates fixed")))
coverage_points$path_role <- factor(
  coverage_points$path_role,
  levels = c("1 coordinate fixed", "d/2 coordinates fixed", "d - 1 coordinates fixed"))
p_coverage <- ggplot2::ggplot(
  coverage_points,
  ggplot2::aes(lambda, coverage, group = subset_size, linetype = path_role,
               colour = path_role)
) +
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                    ymin = 0.95 - mc_half_width, ymax = 0.95 + mc_half_width,
                    fill = "#BDBDBD", alpha = 0.20) +
  ggplot2::geom_hline(yintercept = 0.95, linewidth = 0.45, colour = "#303030") +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(ggplot2::aes(size = effective), alpha = 0.80) +
  ggplot2::facet_wrap(~dimension_label, nrow = 3) +
  ggplot2::scale_size_continuous(range = c(0.8, 2.3), breaks = NULL, guide = "none") +
  ggplot2::scale_linetype_manual(values = c("solid", "22", "42")) +
  ggplot2::scale_colour_manual(values = c(
    "1 coordinate fixed" = unname(simulation_palette["proposed"]),
    "d/2 coordinates fixed" = unname(simulation_palette["oracle"]),
    "d - 1 coordinates fixed" = unname(simulation_palette["dabrowska"]))) +
  ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1)) +
  ggplot2::coord_cartesian(ylim = c(0.75, 1.005)) +
  ggplot2::labs(x = expression(lambda), y = "Empirical coverage", linetype = NULL,
                colour = NULL,
                title = "Pointwise 95% Greenwood intervals (reference omitted)") +
  theme_simulation() +
  ggplot2::theme(legend.position = "bottom")
save_simulation_figure(p_coverage,
                       file.path(figure_dir, "simulation_dimension_coverage.pdf"),
                       width = 7.2, height = 5.4)

print(coverage_summary, row.names = FALSE, digits = 3)

cat("Study 3 complete.\n")
print(summary, row.names = FALSE, digits = 3)


## ===========================================================================
## 4.  Departures from the prime case (Figure S5, Table S1)
## ===========================================================================


theta <- 2
gamma_x <- c(0.2, 0.4)
tau_x <- c(10, 5)
angle <- seq(-log(2), log(2), length.out = 7L)
directions <- cbind(exp(pmax(angle, 0)), exp(pmax(-angle, 0)))
nonreference <- abs(angle) > 1e-12
reference_index <- which.min(abs(angle))
tail_coefficient <- gumbel_tail_copula(c(1, 1), theta)

## Marginal threshold count matching a nominal joint count at the reference.
matched_k <- function(n, censoring, target) {
  beta <- mean(censoring) / (1 - mean(censoring))
  round(((target / tail_coefficient) * n^(2 * beta))^(1 / (1 + 2 * beta)))
}

integrated <- function(estimate, truth) {
  mean((estimate[nonreference] - truth[nonreference])^2, na.rm = TRUE)
}


## -------------------------------------------------------------------------
## (a) Dependent censoring coordinates
## -------------------------------------------------------------------------

censoring_dependence <- c(1, 1.25, 1.5, 2, 2.5, 3, 4)
design_a <- expand.grid(
  censor_theta = censoring_dependence, censoring = c(0.25, 0.40),
  KEEP.OUT.ATTRS = FALSE
)
design_a$n <- 12000L
design_a$k <- mapply(matched_k, design_a$n, design_a$censoring, 50)
truth_rv <- gumbel_truth(directions, theta, gamma_x)

run_a <- function(row, id) {
  cat(sprintf("Study 4(a) %02d/%02d: censoring=%.0f%%, censoring copula theta=%.1f\n",
              id, nrow(design_a), 100 * row$censoring, row$censor_theta))
  out <- mc_lapply(n_rep, function(b) {
    dat <- simulate_burr(row$n, theta, row$censoring, gamma_x = gamma_x,
                         tau_x = tau_x, censor_theta = row$censor_theta)
    fit <- try(weissman_thresholds(dat$z, dat$delta, row$k), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    y <- standardize_multiplicative(dat$z, fit$threshold)
    proposed <- try(directional_km(y, dat$delta, directions, conf_level = 0.95),
                    silent = TRUE)
    if (inherits(proposed, "try-error")) return(NULL)
    covered <- proposed$lower <= truth_rv & truth_rv <= proposed$upper
    c(proposed = integrated(proposed$estimate, truth_rv),
      naive = integrated(naive_tail(y, directions), truth_rv),
      coverage = mean(covered[nonreference], na.rm = TRUE),
      effective = proposed$effective[reference_index],
      censored_pairs = mean(rowSums(dat$delta == 0L) == 2L))
  }, seed = 410000 + 1000 * id, n_cores = n_cores)
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  data.frame(
    censoring = row$censoring, censor_theta = row$censor_theta, n = row$n, k = row$k,
    irmse_proposed = sqrt(mean(out[, "proposed"])),
    irmse_naive = sqrt(mean(out[, "naive"])),
    coverage = mean(out[, "coverage"], na.rm = TRUE),
    effective = mean(out[, "effective"]),
    censored_pairs = mean(out[, "censored_pairs"])
  )
}

path_a <- file.path(result_dir, "study4_dependent_censoring.csv")
if (reuse_results && file.exists(path_a)) {
  result_a <- utils::read.csv(path_a)
} else {
  result_a <- do.call(rbind, lapply(seq_len(nrow(design_a)), function(i) {
    run_a(design_a[i, ], i)
  }))
  result_a$scaled <- sqrt(result_a$effective) * result_a$irmse_proposed
  utils::write.csv(result_a, path_a, row.names = FALSE)
}


## -------------------------------------------------------------------------
## (b) Heavy censoring and the bivariate competitor
## -------------------------------------------------------------------------

design_b <- data.frame(censoring = c(0.10, 0.25, 0.40, 0.50, 0.60))
design_b$n <- 12000L
design_b$k <- mapply(matched_k, design_b$n, design_b$censoring, 50)

run_b <- function(row, id) {
  cat(sprintf("Study 4(b) %02d/%02d: censoring=%.0f%%\n", id, nrow(design_b),
              100 * row$censoring))
  out <- mc_lapply(n_rep, function(b) {
    dat <- simulate_burr(row$n, theta, row$censoring, gamma_x = gamma_x, tau_x = tau_x)
    fit <- try(weissman_thresholds(dat$z, dat$delta, row$k), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    y <- standardize_multiplicative(dat$z, fit$threshold)
    proposed <- try(directional_km(y, dat$delta, directions, conf_level = 0.95),
                    silent = TRUE)
    if (inherits(proposed, "try-error")) return(NULL)
    dabrowska <- try(dabrowska_tail(y, dat$delta, directions), silent = TRUE)
    covered <- truth_rv >= proposed$lower & truth_rv <= proposed$upper
    c(proposed = integrated(proposed$estimate, truth_rv),
      dabrowska = if (inherits(dabrowska, "try-error")) NA_real_ else
        integrated(dabrowska, truth_rv),
      naive = integrated(naive_tail(y, directions), truth_rv),
      coverage = mean(covered[nonreference], na.rm = TRUE),
      effective = proposed$effective[reference_index])
  }, seed = 420000 + 1000 * id, n_cores = n_cores)
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  data.frame(
    censoring = row$censoring, n = row$n, k = row$k,
    irmse_proposed = sqrt(mean(out[, "proposed"])),
    irmse_dabrowska = sqrt(mean(out[, "dabrowska"], na.rm = TRUE)),
    irmse_naive = sqrt(mean(out[, "naive"])),
    coverage = mean(out[, "coverage"], na.rm = TRUE),
    effective = mean(out[, "effective"]),
    dabrowska_failures = mean(is.na(out[, "dabrowska"]))
  )
}

path_b <- file.path(result_dir, "study4_heavy_censoring.csv")
if (reuse_results && file.exists(path_b)) {
  result_b <- utils::read.csv(path_b)
} else {
  result_b <- do.call(rbind, lapply(seq_len(nrow(design_b)), function(i) {
    run_b(design_b[i, ], i)
  }))
  utils::write.csv(result_b, path_b, row.names = FALSE)
}


## -------------------------------------------------------------------------
## (c) Asymptotic independence
##
## The estimator is compared with the exact finite-sample ratio
## p_n(q)/p_n(1,1), which it targets for every n, and with the tail-copula
## limit F^circ, which exists only in the asymptotically dependent cases.
## -------------------------------------------------------------------------

design_c <- data.frame(
  family = c(rep("gumbel", 4), rep("gaussian", 4)),
  parameter = c(2, 1.75, 1.5, 1.25, 0.75, 0.65, 0.55, 0.5),
  stringsAsFactors = FALSE
)
design_c$n <- 12000L
design_c$k <- 600L

run_c <- function(row, id) {
  cat(sprintf("Study 4(c) %02d/%02d: %s(%.2f)\n", id, nrow(design_c), row$family,
              row$parameter))
  finite <- burr_finite_target(directions, row$family, row$parameter, gamma_x,
                               tau_x, row$n, row$k)
  limit <- if (row$family == "gumbel" && row$parameter > 1) {
    gumbel_truth(directions, row$parameter, gamma_x)
  } else {
    rep(NA_real_, nrow(directions))
  }
  out <- mc_lapply(n_rep, function(b) {
    dat <- simulate_burr(row$n, row$parameter, 0.25, gamma_x = gamma_x,
                         tau_x = tau_x, family = row$family)
    fit <- try(weissman_thresholds(dat$z, dat$delta, row$k), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    y <- standardize_multiplicative(dat$z, fit$threshold)
    proposed <- try(directional_km(y, dat$delta, directions, conf_level = 0.95),
                    silent = TRUE)
    if (inherits(proposed, "try-error")) return(NULL)
    covered <- finite$target >= proposed$lower & finite$target <= proposed$upper
    c(finite = integrated(proposed$estimate, finite$target),
      limit = if (all(is.na(limit))) NA_real_ else integrated(proposed$estimate, limit),
      bias = mean((proposed$estimate - finite$target)[nonreference]),
      coverage = mean(covered[nonreference], na.rm = TRUE),
      effective = proposed$effective[reference_index])
  }, seed = 430000 + 1000 * id, n_cores = n_cores)
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  data.frame(
    family = row$family, parameter = row$parameter, n = row$n, k = row$k,
    reference_probability = finite$reference,
    irmse_finite = sqrt(mean(out[, "finite"])),
    irmse_limit = sqrt(mean(out[, "limit"])),
    bias_finite = mean(out[, "bias"]),
    coverage = mean(out[, "coverage"], na.rm = TRUE),
    effective = mean(out[, "effective"]),
    separation = max(abs(finite$target - limit))
  )
}

path_c <- file.path(result_dir, "study4_asymptotic_independence.csv")
if (reuse_results && file.exists(path_c)) {
  result_c <- utils::read.csv(path_c)
} else {
  result_c <- do.call(rbind, lapply(seq_len(nrow(design_c)), function(i) {
    run_c(design_c[i, ], i)
  }))
  utils::write.csv(result_c, path_c, row.names = FALSE)
}


## -------------------------------------------------------------------------
## (d) Heavy censoring in dimensions two to ten
## -------------------------------------------------------------------------

dimension_theta <- 3
dimension_gamma <- 0.5
dimension_censoring <- 0.25
dimensions <- c(2L, 4L, 6L, 8L, 10L)
design_d <- data.frame(d = dimensions, n = 10000L)
## Matched joint counts: k solves k R_{theta,d}(1,...,1) (k/n)^{d beta} = target.
design_d$k <- vapply(dimensions, function(d) {
  beta <- dimension_censoring / (1 - dimension_censoring)
  coefficient <- gumbel_tail_copula(rep(1, d), dimension_theta)
  round(((100 / coefficient) * 10000^(d * beta))^(1 / (1 + d * beta)))
}, numeric(1))

run_d <- function(row, id) {
  d <- row$d
  paths <- lapply(c(1L, max(1L, floor(d / 2)), d - 1L), function(size) {
    lambda <- seq(0, 1, length.out = 11L)
    do.call(rbind, lapply(lambda, function(l) {
      q <- rep(2^l, d)
      q[seq_len(size)] <- 1
      q
    }))
  })
  grid <- unique(do.call(rbind, paths))
  keep <- rowSums(abs(grid - 1) < 1e-12) < d
  truth <- gumbel_truth(grid, dimension_theta, rep(dimension_gamma, d))
  cat(sprintf("Study 4(d) %02d/%02d: d=%d, k=%d\n", id, nrow(design_d), d, row$k))
  out <- mc_lapply(n_rep, function(b) {
    dat <- simulate_pareto(row$n, d, dimension_theta, dimension_censoring,
                           gamma_x = dimension_gamma)
    fit <- try(weissman_thresholds(dat$z, dat$delta, row$k), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    y <- standardize_multiplicative(dat$z, fit$threshold)
    proposed <- try(directional_km(y, dat$delta, grid, conf_level = 0.95),
                    silent = TRUE)
    if (inherits(proposed, "try-error")) return(NULL)
    naive <- naive_tail(y, grid)
    covered <- truth >= proposed$lower & truth <= proposed$upper
    c(proposed = mean((proposed$estimate[keep] - truth[keep])^2),
      naive = mean((naive[keep] - truth[keep])^2, na.rm = TRUE),
      coverage = mean(covered[keep], na.rm = TRUE),
      effective = proposed$effective[which(!keep)[1L]])
  }, seed = 440000 + 1000 * id, n_cores = n_cores)
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  data.frame(
    d = d, n = row$n, k = row$k, censoring = dimension_censoring,
    irmse_proposed = sqrt(mean(out[, "proposed"])),
    irmse_naive = sqrt(mean(out[, "naive"])),
    coverage = mean(out[, "coverage"], na.rm = TRUE),
    effective = mean(out[, "effective"])
  )
}

path_d <- file.path(result_dir, "study4_dimension_censoring.csv")
if (reuse_results && file.exists(path_d)) {
  result_d <- utils::read.csv(path_d)
} else {
  result_d <- do.call(rbind, lapply(seq_len(nrow(design_d)), function(i) {
    run_d(design_d[i, ], i)
  }))
  result_d$scaled <- sqrt(result_d$effective) * result_d$irmse_proposed
  utils::write.csv(result_d, path_d, row.names = FALSE)
}


## -------------------------------------------------------------------------
## Figure
## -------------------------------------------------------------------------

result_a$censor_label <- factor(sprintf("%.0f%%", 100 * result_a$censoring))
p_a <- ggplot2::ggplot(
  result_a, ggplot2::aes(censor_theta, scaled, colour = censor_label,
                         shape = censor_label)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(
    unname(simulation_palette["proposed"]), unname(simulation_palette["matched"])
  )) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(
    x = expression("Censoring copula parameter " * theta[C]),
    y = expression(hat(r)[n] %*% " integrated RMSE"), title = "Dependent censoring"
  ) +
  theme_simulation(8.3) +
  ggplot2::theme(legend.position = c(0.02, 0.02),
                 legend.justification = c(0, 0),
                 legend.background = ggplot2::element_rect(fill = "white", colour = NA),
                 legend.key.size = ggplot2::unit(0.75, "lines"))

relative_b <- rbind(
  data.frame(censoring = result_b$censoring, method = "Dabrowska",
             ratio = result_b$irmse_dabrowska / result_b$irmse_proposed),
  data.frame(censoring = result_b$censoring, method = "Naive",
             ratio = result_b$irmse_naive / result_b$irmse_proposed)
)
p_b <- ggplot2::ggplot(
  relative_b, ggplot2::aes(100 * censoring, ratio, colour = method, shape = method)
) +
  ggplot2::geom_hline(yintercept = 1, linewidth = 0.35, colour = "#606060",
                      linetype = "22") +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(
    "Dabrowska" = unname(simulation_palette["dabrowska"]),
    "Naive" = unname(simulation_palette["naive"])
  )) +
  ggplot2::labs(
    x = "Tail censoring (%)", y = "Integrated RMSE relative to directional",
    title = "Heavy censoring"
  ) +
  theme_simulation(8.3) +
  ggplot2::theme(legend.position = c(0.02, 0.98),
                 legend.justification = c(0, 1),
                 legend.background = ggplot2::element_rect(fill = "white", colour = NA),
                 legend.key.size = ggplot2::unit(0.75, "lines"))

result_d_long <- rbind(
  data.frame(d = result_d$d, method = "Directional", irmse = result_d$irmse_proposed),
  data.frame(d = result_d$d, method = "Naive", irmse = result_d$irmse_naive)
)
p_d <- ggplot2::ggplot(
  result_d_long, ggplot2::aes(d, irmse, colour = method, shape = method)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(
    "Directional" = unname(simulation_palette["proposed"]),
    "Naive" = unname(simulation_palette["naive"])
  )) +
  ggplot2::scale_x_continuous(breaks = dimensions) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(
    x = "Dimension d", y = "Integrated RMSE",
    title = "Dimension at 25% censoring"
  ) +
  theme_simulation(8.3) +
  ggplot2::theme(legend.position = c(0.02, 0.98),
                 legend.justification = c(0, 1),
                 legend.background = ggplot2::element_rect(fill = "white", colour = NA),
                 legend.key.size = ggplot2::unit(0.75, "lines"))

## Asymptotic independence, ordered by decreasing tail dependence.  The estimator
## is measured against the finite-n ratio, which it targets for every n, and the
## point size carries the joint count that collapses as dependence weakens.
result_c$setting <- factor(
  sprintf("%s %s", ifelse(result_c$family == "gumbel", "Gumbel", "Gaussian"),
          as.character(result_c$parameter)),
  levels = c("Gumbel 2", "Gumbel 1.75", "Gumbel 1.5", "Gumbel 1.25",
             "Gaussian 0.75", "Gaussian 0.65", "Gaussian 0.55", "Gaussian 0.5")
)
result_c$dependence <- ifelse(result_c$family == "gumbel",
                              "Asymptotically dependent", "Asymptotically independent")
p_c <- ggplot2::ggplot(
  result_c, ggplot2::aes(setting, irmse_finite, colour = dependence,
                         shape = dependence)
) +
  ggplot2::geom_line(ggplot2::aes(group = dependence), linewidth = 0.6,
                     colour = "#9A9A9A") +
  ggplot2::geom_point(ggplot2::aes(size = effective)) +
  ggplot2::scale_colour_manual(values = c(
    "Asymptotically dependent" = unname(simulation_palette["proposed"]),
    "Asymptotically independent" = unname(simulation_palette["dabrowska"])
  )) +
  ggplot2::scale_size_continuous(range = c(1.2, 3.4), breaks = NULL, guide = "none") +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(x = NULL, y = "Integrated RMSE, finite-n target",
                title = "Asymptotic independence") +
  theme_simulation(8.3) +
  ggplot2::theme(legend.position = c(0.02, 0.02),
                 legend.justification = c(0, 0),
                 legend.background = ggplot2::element_rect(fill = "white", colour = NA),
                 legend.key.size = ggplot2::unit(0.75, "lines"),
                 axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

figure <- cowplot::plot_grid(
  p_a, p_b, p_c, p_d, labels = c("(a)", "(b)", "(c)", "(d)"), nrow = 2,
  label_size = panel_label_size, label_fontface = "plain",
  align = "hv", axis = "tblr"
)
save_simulation_figure(
  figure, file.path(figure_dir, "simulation_robustness.pdf"),
  width = 7.2, height = 5.4
)

print(result_a)
print(result_b)
print(result_c)
print(result_d)
