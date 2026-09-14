rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

# User-facing posterior summary helpers for RFBS log files.
# This script is intentionally small and self-contained. It is meant to help a
# user inspect one RFBS run directory without needing the paper-era plotting scripts.

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
  post_burn <- log_df[start_row:nrow(log_df), , drop = FALSE]

  list(
    log_file = log_file,
    full = log_df,
    post_burn = post_burn,
    burnin = burnin,
    start_row = start_row
  )
}

get_rate_columns <- function(log_df) {
  cols <- colnames(log_df)
  rate_cols <- grep("^_", cols, value = TRUE)
  clado_cols <- intersect(c("split", "sub", "equal"), cols)
  c(rate_cols, clado_cols)
}

get_rj_columns <- function(log_df) {
  grep("^switch_", colnames(log_df), value = TRUE)
}

summarize_posterior_columns <- function(log_df, cols) {
  if (length(cols) == 0) {
    return(data.frame())
  }

  do.call(rbind, lapply(cols, function(col) {
    x <- as.numeric(log_df[[col]])
    x <- x[is.finite(x)]
    data.frame(
      parameter = col,
      mean = mean(x),
      median = median(x),
      sd = stats::sd(x),
      q05 = stats::quantile(x, 0.05),
      q25 = stats::quantile(x, 0.25),
      q75 = stats::quantile(x, 0.75),
      q95 = stats::quantile(x, 0.95),
      nonzero_fraction = mean(x != 0),
      stringsAsFactors = FALSE
    )
  }))
}

summarize_rj_columns <- function(log_df, switch_cols) {
  if (length(switch_cols) == 0) {
    return(data.frame())
  }

  do.call(rbind, lapply(switch_cols, function(col) {
    x <- as.numeric(log_df[[col]])
    x <- x[is.finite(x)]
    data.frame(
      switch = col,
      posterior_on_fraction = mean(x == 1),
      posterior_off_fraction = mean(x == 0),
      mean_value = mean(x),
      stringsAsFactors = FALSE
    )
  }))
}

plot_rate_densities <- function(log_df, cols, outfile, max_cols = 12) {
  cols <- cols[seq_len(min(length(cols), max_cols))]
  if (length(cols) == 0) {
    return(invisible(NULL))
  }

  long_df <- do.call(rbind, lapply(cols, function(col) {
    x <- as.numeric(log_df[[col]])
    x <- x[is.finite(x)]
    data.frame(parameter = col, value = x)
  }))

  p <- ggplot(long_df, aes(x = value)) +
    geom_density(fill = "steelblue", alpha = 0.35, linewidth = 0.6) +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    theme_bw() +
    labs(
      title = "RFBS posterior rate distributions",
      subtitle = "Post-burnin samples",
      x = "Posterior value",
      y = "Density"
    )

  ggsave(outfile, p, width = 10, height = 8)
}

plot_rj_bars <- function(rj_summary_df, outfile) {
  if (nrow(rj_summary_df) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot(rj_summary_df, aes(x = reorder(switch, posterior_on_fraction), y = posterior_on_fraction)) +
    geom_col(fill = "firebrick3") +
    coord_flip() +
    theme_bw() +
    labs(
      title = "RJ posterior ON frequencies",
      subtitle = "Fraction of post-burnin samples with switch = 1",
      x = "RJ switch",
      y = "Posterior ON fraction"
    )

  ggsave(outfile, p, width = 8, height = 5)
}

plot_rj_joint_heatmap <- function(log_df, switch_cols, outfile) {
  if (length(switch_cols) < 2) {
    return(invisible(NULL))
  }

  switch_df <- log_df[, switch_cols, drop = FALSE]
  switch_df[] <- lapply(switch_df, function(x) as.numeric(x) == 1)

  combos <- expand.grid(row = switch_cols, col = switch_cols, stringsAsFactors = FALSE)
  combos$joint_on_fraction <- mapply(function(a, b) {
    mean(switch_df[[a]] & switch_df[[b]])
  }, combos$row, combos$col)

  p <- ggplot(combos, aes(x = col, y = row, fill = joint_on_fraction)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", joint_on_fraction)), size = 3) +
    scale_fill_gradient(low = "white", high = "navy") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "RJ joint ON frequencies",
      subtitle = "Fraction of post-burnin samples where both switches are ON",
      x = NULL,
      y = NULL,
      fill = "Joint ON"
    )

  ggsave(outfile, p, width = 8, height = 6)
}

write_user_friendly_postprocess <- function(out_dir, burnin = 0.5, summary_dir = file.path(out_dir, "posterior_summary")) {
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

  log_obj <- read_rfbs_log(out_dir, burnin = burnin)
  post_burn <- log_obj$post_burn

  rate_cols <- get_rate_columns(post_burn)
  switch_cols <- get_rj_columns(post_burn)

  rate_summary <- summarize_posterior_columns(post_burn, rate_cols)
  rj_summary <- summarize_rj_columns(post_burn, switch_cols)

  write.csv(rate_summary, file.path(summary_dir, "posterior_rate_summary.csv"), row.names = FALSE)
  write.csv(rj_summary, file.path(summary_dir, "posterior_rj_summary.csv"), row.names = FALSE)

  plot_rate_densities(post_burn, rate_cols, file.path(summary_dir, "posterior_rate_densities.pdf"))
  plot_rj_bars(rj_summary, file.path(summary_dir, "posterior_rj_on_frequencies.pdf"))
  plot_rj_joint_heatmap(post_burn, switch_cols, file.path(summary_dir, "posterior_rj_joint_on_heatmap.pdf"))

  cat("Posterior summaries written to:\n")
  cat(summary_dir, "\n\n")
  cat("Files created:\n")
  cat("- posterior_rate_summary.csv\n")
  cat("- posterior_rj_summary.csv\n")
  cat("- posterior_rate_densities.pdf\n")
  cat("- posterior_rj_on_frequencies.pdf\n")
  cat("- posterior_rj_joint_on_heatmap.pdf\n")

  invisible(list(
    summary_dir = summary_dir,
    rate_summary = rate_summary,
    rj_summary = rj_summary
  ))
}


# ----------------------------
# Example usage
# ----------------------------
# Set this to any RFBS output directory containing RFBS_log.
out_dir <- file.path(getwd(), "RFBS_user_example_output")

# Burnin is the fraction of rows dropped from the start of the chain.
burnin <- 0.5

write_user_friendly_postprocess(out_dir = out_dir, burnin = burnin)
