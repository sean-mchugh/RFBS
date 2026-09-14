rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(vioplot)
})

# ------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------

# One or a few RFBS output directories containing RFBS_log.
# Name the runs so those names appear in plot legends.
run_dirs <- c(
  example = file.path(getwd(), "RFBS_user_example_output")
)

# Fraction of the chain to discard from the start.
burnin <- 0.5

# If you have a prior-only RFBS_log object or file, you can wire it in later.
# Leave as NULL to skip prior overlays.
prior_log_file <- NULL

# Optional output directory. If NULL and only one run is supplied, results go to
# <run_dir>/publication_summary. If multiple runs are supplied, results go to a
# shared directory in the working directory.
summary_dir <- NULL

# Optional run colors. If NULL, defaults are assigned automatically.
run_colors <- NULL


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  rgb(
    rgb.val[1], rgb.val[2], rgb.val[3],
    max = 255,
    alpha = (100 - percent) * 255 / 100,
    names = name
  )
}

default_run_colors <- function(n) {
  base_cols <- c(
    "ivory", "linen", "blue4", "yellow", "yellow3", "red1",
    "red4", "green1", "green4", "purple1", "purple4",
    "orange1", "darkorange2", "orange3", "darkorange4",
    "gray70", "gray50", "gray30", "gray0"
  )
  rep(base_cols, length.out = n)
}

resolve_summary_dir <- function(run_dirs, summary_dir) {
  if (!is.null(summary_dir)) {
    return(summary_dir)
  }
  if (length(run_dirs) == 1) {
    return(file.path(unname(run_dirs[[1]]), "publication_summary"))
  }
  file.path(getwd(), "RFBS_publication_summary")
}

read_rfbs_log <- function(out_dir, burnin = 0.5) {
  log_file <- file.path(out_dir, "RFBS_log")
  if (!file.exists(log_file)) {
    stop("Could not find RFBS_log in: ", out_dir)
  }

  log_df <- read.table(log_file, header = TRUE, sep = "\t", check.names = FALSE)
  if (nrow(log_df) == 0) {
    stop("RFBS_log is empty: ", log_file)
  }

  start_row <- max(1, floor(nrow(log_df) * burnin) + 1)
  log_df[start_row:nrow(log_df), , drop = FALSE]
}

preferred_rate_order <- c(
  "_rf1_l_s",
  "_rf2_l_d",
  "_rf1_g_s",
  "_rf12_sw_s",
  "_rf2_l_s",
  "_rf2_g_d",
  "_rf2_g_s",
  "split",
  "sub",
  "equal"
)

preferred_switch_order <- c(
  "switch_rf1_l_s",
  "switch_rf2_l_d",
  "switch_rf1_g_s",
  "switch_rf12_sw_s",
  "switch_rf2_l_s",
  "switch_rf2_g_d",
  "switch_rf2_g_s",
  "switch_split",
  "switch_sub",
  "switch_equal"
)

parameter_labels <- c(
  "_rf1_l_s" = "Enabled loss event",
  "_rf2_l_d" = "Double loss event",
  "_rf1_g_s" = "Enabled gain event",
  "_rf12_sw_s" = "Established switch event",
  "_rf2_l_s" = "Established loss event",
  "_rf2_g_d" = "Double gain event",
  "_rf2_g_s" = "Established gain event",
  "split" = "Speciation split event",
  "sub" = "Speciation subset event",
  "equal" = "Speciation equal event",
  "switch_rf1_l_s" = "Enabled loss RJ",
  "switch_rf2_l_d" = "Double loss RJ",
  "switch_rf1_g_s" = "Enabled gain RJ",
  "switch_rf12_sw_s" = "Established switch RJ",
  "switch_rf2_l_s" = "Established loss RJ",
  "switch_rf2_g_d" = "Double gain RJ",
  "switch_rf2_g_s" = "Established gain RJ",
  "switch_split" = "Split RJ",
  "switch_sub" = "Subset RJ",
  "switch_equal" = "Equal RJ"
)

available_rate_columns <- function(log_df) {
  intersect(preferred_rate_order, colnames(log_df))
}

available_switch_columns <- function(log_df) {
  intersect(preferred_switch_order, colnames(log_df))
}

collect_run_data <- function(run_dirs, burnin) {
  run_names <- names(run_dirs)
  if (is.null(run_names) || any(run_names == "")) {
    run_names <- paste0("run_", seq_along(run_dirs))
    names(run_dirs) <- run_names
  }

  run_logs <- lapply(run_dirs, read_rfbs_log, burnin = burnin)
  list(run_dirs = run_dirs, run_logs = run_logs, run_names = names(run_dirs))
}

stack_parameter_draws <- function(run_logs, cols) {
  do.call(rbind, lapply(names(run_logs), function(run_name) {
    log_df <- run_logs[[run_name]]
    do.call(rbind, lapply(cols, function(col) {
      x <- as.numeric(log_df[[col]])
      x <- x[is.finite(x)]
      data.frame(
        run = run_name,
        parameter = col,
        parameter_label = ifelse(col %in% names(parameter_labels), parameter_labels[[col]], col),
        value = x,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

posterior_summary_table <- function(draw_df) {
  split_df <- split(draw_df, list(draw_df$run, draw_df$parameter), drop = TRUE)
  do.call(rbind, lapply(split_df, function(df) {
    data.frame(
      run = df$run[[1]],
      parameter = df$parameter[[1]],
      parameter_label = df$parameter_label[[1]],
      mean = mean(df$value),
      median = median(df$value),
      sd = sd(df$value),
      q05 = as.numeric(quantile(df$value, 0.05)),
      q25 = as.numeric(quantile(df$value, 0.25)),
      q75 = as.numeric(quantile(df$value, 0.75)),
      q95 = as.numeric(quantile(df$value, 0.95)),
      nonzero_fraction = mean(df$value != 0),
      stringsAsFactors = FALSE
    )
  }))
}

rj_summary_table <- function(run_logs, switch_cols) {
  do.call(rbind, lapply(names(run_logs), function(run_name) {
    log_df <- run_logs[[run_name]]
    do.call(rbind, lapply(switch_cols, function(col) {
      x <- as.numeric(log_df[[col]])
      x <- x[is.finite(x)]
      data.frame(
        run = run_name,
        switch = col,
        switch_label = ifelse(col %in% names(parameter_labels), parameter_labels[[col]], col),
        posterior_on_fraction = mean(x == 1),
        posterior_off_fraction = mean(x == 0),
        mean_value = mean(x),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

save_violin_page <- function(draw_df, outfile, run_colors) {
  params <- unique(draw_df$parameter)

  pdf(outfile, width = 10, height = 10)
  on.exit(dev.off(), add = TRUE)

  par(oma = c(6, 1, 1, 1))
  par(mar = c(5, 2, 2, 2))
  n_panels <- length(params)
  par(mfrow = c(n_panels, 1))

  run_levels <- unique(draw_df$run)
  for (param in params) {
    sub_df <- draw_df[draw_df$parameter == param, , drop = FALSE]
    split_vals <- split(sub_df$value, sub_df$run)
    split_vals <- split_vals[run_levels]
    split_vals <- unname(split_vals)

    do.call(
      vioplot,
      c(
        split_vals,
        list(
          col = run_colors[seq_along(run_levels)],
          names = rep("", length(run_levels)),
          main = unique(sub_df$parameter_label),
          cex.main = 1.4,
          cex.axis = 1.0,
          las = 2
        )
      )
    )
  }

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(
    "bottom",
    legend = run_levels,
    col = "black",
    pt.bg = run_colors[seq_along(run_levels)],
    pch = 21,
    pt.cex = 2,
    cex = 1.4,
    seg.len = 0.25,
    ncol = max(1, min(4, length(run_levels))),
    bty = "n"
  )
}

save_density_page <- function(draw_df, outfile, run_colors, prior_draw_df = NULL) {
  plot_df <- draw_df
  plot_df$run <- factor(plot_df$run, levels = unique(draw_df$run))

  p <- ggplot(plot_df, aes(x = value, color = run, fill = run)) +
    geom_density(alpha = 0.10, linewidth = 1.0) +
    facet_wrap(~ parameter_label, scales = "free", ncol = 2) +
    scale_color_manual(values = run_colors[seq_along(unique(draw_df$run))]) +
    scale_fill_manual(values = run_colors[seq_along(unique(draw_df$run))]) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = "RFBS posterior densities",
      subtitle = "Post-burnin posterior samples",
      x = "Posterior value",
      y = "Density"
    )

  if (!is.null(prior_draw_df) && nrow(prior_draw_df) > 0) {
    p <- p +
      geom_density(
        data = prior_draw_df,
        aes(x = value),
        inherit.aes = FALSE,
        color = "gray40",
        linewidth = 0.9,
        linetype = "dashed"
      )
  }

  ggsave(outfile, p, width = 10, height = 10)
}

joint_on_fractions <- function(log_df, switch_cols) {
  switch_df <- log_df[, switch_cols, drop = FALSE]
  switch_df[] <- lapply(switch_df, function(x) as.numeric(x) == 1)
  combos <- expand.grid(row = switch_cols, col = switch_cols, stringsAsFactors = FALSE)
  combos$joint_on_fraction <- mapply(function(a, b) {
    mean(switch_df[[a]] & switch_df[[b]])
  }, combos$row, combos$col)
  combos
}

save_rj_joint_heatmap <- function(log_df, switch_cols, outfile, title_prefix = "RJ") {
  if (length(switch_cols) < 2) {
    return(invisible(NULL))
  }

  combos <- joint_on_fractions(log_df, switch_cols)
  combos$row_label <- ifelse(combos$row %in% names(parameter_labels), parameter_labels[combos$row], combos$row)
  combos$col_label <- ifelse(combos$col %in% names(parameter_labels), parameter_labels[combos$col], combos$col)

  p <- ggplot(combos, aes(x = col_label, y = row_label, fill = joint_on_fraction)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", joint_on_fraction)), size = 3) +
    scale_fill_gradient(low = "white", high = "navy") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(
      title = paste(title_prefix, "joint ON frequencies"),
      subtitle = "Fraction of post-burnin samples where both switches are ON",
      fill = "Joint ON"
    )

  ggsave(outfile, p, width = 10, height = 8)
}

save_rj_pairwise_quadrants <- function(log_df, switch_cols, outfile, title_prefix = "RJ") {
  if (length(switch_cols) < 2) {
    return(invisible(NULL))
  }

  pairs <- combn(switch_cols, 2, simplify = FALSE)
  switch_df <- log_df[, switch_cols, drop = FALSE]
  switch_df[] <- lapply(switch_df, function(x) as.numeric(x) == 1)

  plots <- lapply(pairs, function(pair) {
    a <- pair[[1]]
    b <- pair[[2]]
    rj1 <- switch_df[[a]]
    rj2 <- switch_df[[b]]
    joints <- c(
      mean(rj1 & rj2),
      mean(rj1 & !rj2),
      mean(!rj1 & rj2),
      mean(!rj1 & !rj2)
    )
    dataframe <- data.frame(matrix(c(0, 1, 2, 0, 1, 2), ncol = 2))

    ggplot(dataframe, aes(x = X1, y = X2)) +
      geom_rect(xmin = 1, xmax = 2, ymin = 1, ymax = 2, color = "black", fill = "blue", alpha = joints[1]) +
      geom_rect(xmin = 1, xmax = 2, ymin = 0, ymax = 1, color = "black", fill = "blue", alpha = joints[2]) +
      geom_rect(xmin = 0, xmax = 1, ymin = 1, ymax = 2, color = "black", fill = "blue", alpha = joints[3]) +
      geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 1, color = "black", fill = "blue", alpha = joints[4]) +
      annotate(geom = "text", x = 1.5, y = 1.5, label = sprintf("%.2f", joints[1]), size = 5) +
      annotate(geom = "text", x = 1.5, y = 0.5, label = sprintf("%.2f", joints[2]), size = 5) +
      annotate(geom = "text", x = 0.5, y = 1.5, label = sprintf("%.2f", joints[3]), size = 5) +
      annotate(geom = "text", x = 0.5, y = 0.5, label = sprintf("%.2f", joints[4]), size = 5) +
      scale_x_continuous(limits = c(0, 2), breaks = c(0.5, 1.5), labels = c("Off", "On")) +
      scale_y_continuous(limits = c(0, 2), breaks = c(0.5, 1.5), labels = c("Off", "On")) +
      labs(
        title = paste0(
          ifelse(a %in% names(parameter_labels), parameter_labels[[a]], a),
          " vs.\n",
          ifelse(b %in% names(parameter_labels), parameter_labels[[b]], b)
        ),
        x = NULL,
        y = NULL
      ) +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = -0.3),
        axis.text.x = element_text(vjust = 1.5)
      )
  })

  full_plot <- wrap_plots(plots, ncol = 2) +
    plot_annotation(title = paste(title_prefix, "pairwise RJ ON/OFF plots"))

  ggsave(
    outfile,
    full_plot,
    width = 10,
    height = max(8, 3 * ceiling(length(plots) / 2)),
    limitsize = FALSE
  )
}

save_rj_big_onoff_plot <- function(log_df, switch_cols, outfile, title_prefix = "RJ") {
  if (length(switch_cols) == 0) {
    return(invisible(NULL))
  }

  switch_df <- log_df[, switch_cols, drop = FALSE]
  switch_df[] <- lapply(switch_df, function(x) as.numeric(x) == 1)

  sums <- rowSums(switch_df)
  onoff_sums <- cbind(switch_df, sums)
  onoff <- data.frame(sums = sums)
  total <- nrow(switch_df)

  total_df <- do.call(rbind, lapply(switch_cols, function(param) {
    data.frame(
      param = param,
      percent = mean(switch_df[[param]]),
      stringsAsFactors = FALSE
    )
  }))

  stacked_df <- do.call(rbind, lapply(0:length(switch_cols), function(i) {
    onoff_sums_i <- onoff_sums[onoff_sums[, ncol(onoff_sums)] == i, , drop = FALSE]
    if (nrow(onoff_sums_i) == 0) {
      return(NULL)
    }
    do.call(rbind, lapply(switch_cols, function(param) {
      data.frame(
        onoff = i,
        param = param,
        count = sum(onoff_sums_i[[param]]),
        stringsAsFactors = FALSE
      )
    }))
  }))

  levels <- switch_cols
  labels <- ifelse(levels %in% names(parameter_labels), parameter_labels[levels], levels)
  palette_vals <- rep(RColorBrewer::brewer.pal(min(8, max(3, length(levels))), "Dark2"), length.out = length(levels))

  total_plot <- ggplot(total_df, aes(x = factor(param, levels = levels), y = percent, fill = factor(param, levels = levels))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(labels = labels, values = palette_vals) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(fill = "Parameter", x = "Parameter", y = "Frequency") +
    theme_bw() +
    theme(
      aspect.ratio = 0.4,
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(8, 8, 4, 8)
    )

  onoff_plot <- ggplot(onoff, aes(sums)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5, color = "black", linewidth = 0.1, fill = "white") +
    geom_density(aes(linetype = "Observed"), adjust = 4, linewidth = 1) +
    geom_density(aes(rbinom(nrow(onoff), length(switch_cols), 0.5), linetype = "Prior"), adjust = 4, linewidth = 1) +
    scale_linetype_manual(values = c("Observed" = "solid", "Prior" = "dashed")) +
    scale_x_continuous(breaks = seq(0, length(switch_cols), 1), limits = c(-0.5, length(switch_cols) + 0.5), expand = c(0, 0)) +
    labs(y = "Density", x = NULL, linetype = "Distribution") +
    theme_bw() +
    theme(
      aspect.ratio = 0.3,
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(4, 8, 4, 8)
    )

  stacked_plot <- ggplot(stacked_df, aes(x = factor(onoff), y = count, fill = factor(param, levels = levels))) +
    geom_bar(stat = "identity", position = "fill", width = 1, color = "black", linewidth = 0.1) +
    scale_fill_manual(labels = labels, values = palette_vals) +
    labs(fill = "Parameter", x = "Number of 'ON' Parameters", y = "Representation By Bin") +
    theme_bw() +
    theme(
      aspect.ratio = 0.75,
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(4, 8, 8, 8)
    )

  big_plot <- total_plot + onoff_plot + stacked_plot + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(title = paste(title_prefix, "ON/OFF summaries"))

  ggsave(outfile, big_plot, width = 10, height = 10)
}

maybe_prior_draws <- function(prior_log_file, cols) {
  if (is.null(prior_log_file) || !file.exists(prior_log_file)) {
    return(NULL)
  }
  prior_df <- read.table(prior_log_file, header = TRUE, sep = "\t", check.names = FALSE)
  if (nrow(prior_df) == 0) {
    return(NULL)
  }
  stack_parameter_draws(list(prior = prior_df), intersect(cols, colnames(prior_df)))
}


# ------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("Package 'RColorBrewer' is required.")
}

run_data <- collect_run_data(run_dirs, burnin = burnin)
summary_dir <- resolve_summary_dir(run_data$run_dirs, summary_dir)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

if (is.null(run_colors)) {
  run_colors <- default_run_colors(length(run_data$run_names))
}
names(run_colors) <- run_data$run_names

first_log <- run_data$run_logs[[1]]
rate_cols <- available_rate_columns(first_log)
switch_cols <- available_switch_columns(first_log)

rate_draws <- stack_parameter_draws(run_data$run_logs, rate_cols)
rate_summary <- posterior_summary_table(rate_draws)
rj_summary <- rj_summary_table(run_data$run_logs, switch_cols)

write.csv(rate_summary, file.path(summary_dir, "posterior_parameter_summary.csv"), row.names = FALSE)
write.csv(rj_summary, file.path(summary_dir, "posterior_rj_summary.csv"), row.names = FALSE)

prior_draws <- maybe_prior_draws(prior_log_file, rate_cols)

save_violin_page(rate_draws, file.path(summary_dir, "posterior_density_violin.pdf"), run_colors = unname(run_colors))
save_density_page(rate_draws, file.path(summary_dir, "posterior_density_panels.pdf"), run_colors = unname(run_colors), prior_draw_df = prior_draws)

for (run_name in run_data$run_names) {
  log_df <- run_data$run_logs[[run_name]]
  if (length(switch_cols) > 0) {
    save_rj_pairwise_quadrants(
      log_df,
      switch_cols,
      file.path(summary_dir, paste0(run_name, "_RJ_pairwise_quadrants.pdf")),
      title_prefix = run_name
    )
    save_rj_joint_heatmap(
      log_df,
      switch_cols,
      file.path(summary_dir, paste0(run_name, "_RJ_joint_heatmap.pdf")),
      title_prefix = run_name
    )
    save_rj_big_onoff_plot(
      log_df,
      switch_cols,
      file.path(summary_dir, paste0(run_name, "_RJ_onoff_summary.pdf")),
      title_prefix = run_name
    )
  }
}

cat("Publication-style summaries written to:\n")
cat(summary_dir, "\n\n")
cat("Files created:\n")
cat("- posterior_parameter_summary.csv\n")
cat("- posterior_rj_summary.csv\n")
cat("- posterior_density_violin.pdf\n")
cat("- posterior_density_panels.pdf\n")
if (length(switch_cols) > 0) {
  for (run_name in run_data$run_names) {
    cat("-", paste0(run_name, "_RJ_pairwise_quadrants.pdf"), "\n")
    cat("-", paste0(run_name, "_RJ_joint_heatmap.pdf"), "\n")
    cat("-", paste0(run_name, "_RJ_onoff_summary.pdf"), "\n")
  }
} else {
  cat("- no RJ plots were created because no switch_* columns were present\n")
}
