rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggtext)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

## Figure S1: posterior of the remote-onsite shift (TL-sGMM, SURMOUNT-OSA) ====================
# real_data/results/shift_posterior.rds holds only the delta_beta draws, spline knots and arm
# of the four TL-sGMM fits (see real_data/extract_shift_posterior.R); no patient-level data.
shift_posterior <- readRDS("real_data/results/shift_posterior.rds")
df_db <- list()
for (ii in seq_along(shift_posterior)) {
  res <- shift_posterior[[ii]]
  db <- res$delta_beta
  dlt <- "\u03B4"
  colnames(db) <- c(
    paste0(dlt, "<sub>\u03B2<sub>0</sub></sub>"),
    paste0(dlt, "<sub>\u03B2<sub>1</sub></sub>"),
    paste0("**", dlt, "**<sub>b</sub>[", seq_len(ncol(db) - 2), "]")
  )
  df_db[[ii]] <- as.data.frame(db) %>%
    pivot_longer(everything(), names_to = "coef", values_to = "delta") %>%
    mutate(
      treatment = res$treatment,
      coef = factor(coef, levels = colnames(db))
    )
}

### shift function f1(t) - f0(t)
# remote-onsite shift function delta(t) = delta_beta0 + delta_beta1 * t + B1(t) delta_b
t_grid <- seq(0, 12, by = 0.1)
df_shift <- list()
for (ii in seq_along(shift_posterior)) {
  res <- shift_posterior[[ii]]
  db <- res$delta_beta
  knots <- res$knots
  sqrt_omega <- with(svd(abs(outer(knots, knots, "-"))^3), t(v %*% (t(u) * sqrt(d))))
  B_grid <- t(solve(sqrt_omega, t(abs(outer(t_grid, knots, "-"))^3)))
  d_draw <- db[, 1] + outer(db[, 2], t_grid) + db[, -(1:2)] %*% t(B_grid)
  df_shift[[ii]] <- tibble(
    time = t_grid * 30 / 7,
    est = colMeans(d_draw),
    lo = apply(d_draw, 2, quantile, 0.025),
    hi = apply(d_draw, 2, quantile, 0.975),
    treatment = res$treatment,
    truth = if (ii > 2 && res$treatment == "TZP MTD") -2 - 2 * t_grid else 0
  )
}

plot_fig2c <- function(dat, shift, title, ylim = 100, ylim_pbo = ylim) {
  dat <- dat %>%
    mutate(treatment = factor(treatment,
      levels = c("TZP MTD", "Placebo"),
      labels = c("Tirzepatide", "Placebo")
    ))
  shift <- shift %>%
    mutate(
      treatment = factor(treatment,
        levels = c("TZP MTD", "Placebo"),
        labels = c("Tirzepatide", "Placebo")
      ),
      coef = "Week"
    )
  lv <- levels(dat$coef)
  if (is.null(lv)) lv <- unique(as.character(dat$coef))
  ylim <- rep(ylim, length.out = 2)
  ylim_pbo <- rep(ylim_pbo, length.out = 2)
  # facet panel.spacing relative to one panel width at 10 x 7 in; keeps all panels equally wide
  gap <- 0.08
  mk <- function(d, yl, yl_pbo, ttl, ystrip) {
    y_top <- data.frame(
      treatment = factor(c("Tirzepatide", "Placebo"), levels = c("Tirzepatide", "Placebo")),
      y = c(yl, yl_pbo)
    )
    ggplot(d, aes(delta, fill = treatment, color = treatment)) +
      geom_density(alpha = 0.4, show.legend = FALSE) +
      geom_blank(data = y_top, aes(y = y), inherit.aes = FALSE) +
      geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
      scale_fill_manual(values = c("Tirzepatide" = "#377EB8", "Placebo" = "#66C2A5")) +
      scale_color_manual(values = c("Tirzepatide" = "#377EB8", "Placebo" = "#66C2A5")) +
      scale_x_continuous(
        breaks = scales::breaks_extended(n = 3),
        labels = function(x) format(round(x, 3), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
      ) +
      scale_y_continuous(limits = c(0, NA)) +
      facet_grid(treatment ~ coef, scales = "free", switch = "x") +
      labs(x = NULL, y = "Density", title = ttl) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = if (ystrip) element_text(face = "bold") else element_blank(),
        strip.text.x = ggtext::element_markdown(family = "Times New Roman"),
        strip.placement = "outside",
        axis.text.x = element_text(size = 7),
        plot.title = element_text(face = "plain", hjust = 0, size = 14)
      )
  }
  mk_shift <- ggplot(shift, aes(x = time)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = treatment), alpha = 0.4) +
    geom_line(aes(y = est, color = treatment)) +
    geom_line(data = filter(shift, truth != 0), aes(y = truth), linetype = 2, linewidth = 0.3, color = "#EF6161") +
    scale_fill_manual(values = c("Tirzepatide" = "#377EB8", "Placebo" = "#66C2A5")) +
    scale_color_manual(values = c("Tirzepatide" = "#377EB8", "Placebo" = "#66C2A5")) +
    scale_x_continuous(breaks = c(seq(0, 40, 8), 52)) +
    facet_grid(treatment ~ coef, scales = "free_y", switch = "x") +
    geom_blank(data = filter(shift, truth == 0), aes(y = if_else(time < 26, -8, 8))) +
    labs(x = NULL, y = expression(italic(f)[1](italic(t)) - italic(f)[0](italic(t)))) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_text(face = "bold"),
      strip.text.x = ggtext::element_markdown(family = "Times New Roman"),
      strip.placement = "outside",
      axis.text.x = element_text(size = 7)
    )
  mk(filter(dat, coef %in% lv[1:2]), ylim[1], ylim_pbo[1], title, FALSE) +
    mk(filter(dat, coef %in% lv[-(1:2)]), ylim[2], ylim_pbo[2], NULL, FALSE) +
    mk_shift +
    plot_layout(widths = c(2 + gap, 4 + 3 * gap, 1))
}
fig2c <- (plot_fig2c(bind_rows(df_db[1:2]), bind_rows(df_shift[1:2]), "(a) Observed Data", ylim = c(55, 95)) /
  plot_fig2c(bind_rows(df_db[3:4]), bind_rows(df_shift[3:4]), "(b) Induced shift", ylim = c(10, 28), ylim_pbo = c(55, 95)))
ggsave(
  filename = "output/figures/figure_S1.pdf",
  plot = fig2c, width = 11, height = 7, dpi = 800,
  device = if (capabilities("aqua")) {
    function(filename, width, height, ...) quartz(file = filename, type = "pdf", width = width, height = height)
  } else {
    cairo_pdf
  }
)
