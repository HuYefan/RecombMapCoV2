# =============================================================================
#  03_plot_events_over_time.R
#  Supporting analysis: when did the recombination events occur?  Monthly
#  counts and cumulative number of new recombination nodes, for all events and
#  for supported events (best recomb_parsimony >= 7).
# =============================================================================
#
#  Input   : data/derived/recombination_events_summary.tsv
#            (produced by 01_annotate_recombination_dates.py)
#  Output  : figures/main/FigureS2_recombination_events_over_time.pdf
#
#  Run from the repository root:
#      Rscript Figure_1/cd_recombination_breakpoints/scripts/03_plot_events_over_time.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("config.R")

events_file <- file.path(DIR_DATA_DRV, "recombination_events_summary.tsv")
if (!file.exists(events_file)) {
  # fall back to the compact node table distributed with the repository
  events_file <- file.path(DIR_DATA_RAW, "descendants_node_dates.tsv")
  events <- read.delim(events_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "")
  names(events)[1] <- "recomb_node_id"
  events$best_recomb_parsimony <- NA_integer_
} else {
  events <- read.delim(events_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
}

events$earliest_date <- as.Date(events$earliest_date)
events <- events[!is.na(events$earliest_date) &
                   events$earliest_date >= as.Date("2019-12-01") &
                   events$earliest_date <= as.Date("2025-12-31"), , drop = FALSE]

if (all(is.na(events$best_recomb_parsimony))) {
  events$supported <- TRUE
} else {
  events$supported <- events$best_recomb_parsimony >= RECOMB_PARSIMONY_MIN
}
events$month <- as.Date(format(events$earliest_date, "%Y-%m-01"))

message("events in range: ", nrow(events),
        " (supported: ", sum(events$supported), ")")

# -- monthly counts (all vs. supported) --------------------------------------
count_by <- function(df, label) {
  agg <- aggregate(list(n = df$recomb_node_id), by = list(month = df$month), FUN = length)
  agg$set <- label
  agg
}
monthly <- rbind(count_by(events, "All events"),
                 count_by(events[events$supported, , drop = FALSE], "Supported"))

# fill missing months with 0 so the two series line up
all_months <- seq(min(monthly$month), max(monthly$month), by = "month")
monthly <- merge(expand.grid(month = all_months,
                             set = c("All events", "Supported")),
                 monthly, by = c("month", "set"), all.x = TRUE)
monthly$n[is.na(monthly$n)] <- 0

p1 <- ggplot(monthly, aes(x = month, y = n, fill = set)) +
  geom_col(position = position_dodge(width = 20)) +
  scale_fill_manual(values = c("All events" = "grey60", "Supported" = "steelblue")) +
  labs(x = "Month", y = "Number of new\nrecombination nodes", fill = NULL) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.15, 0.85))

# -- cumulative counts -------------------------------------------------------
cum_by <- function(df, label) {
  d <- df[order(df$earliest_date), , drop = FALSE]
  data.frame(date = d$earliest_date,
             cum = seq_len(nrow(d)),
             set = label)
}
cumulative <- rbind(cum_by(events, "All events"),
                    cum_by(events[events$supported, , drop = FALSE], "Supported"))

p2 <- ggplot(cumulative, aes(x = date, y = cum, colour = set)) +
  geom_step(linewidth = 0.5) +
  scale_colour_manual(values = c("All events" = "grey60", "Supported" = "steelblue")) +
  labs(x = "Date", y = "Cumulative number of\nrecombination nodes", colour = NULL) +
  theme_ms() +
  theme(legend.position = c(0.2, 0.85))

p <- p1 / p2 +
  plot_annotation(
    title = sprintf("Recombination events over time (recomb_parsimony >= %d, n = %d)",
                    RECOMB_PARSIMONY_MIN, sum(events$supported))
  )

print(p)
ggsave(p, filename = file.path(DIR_FIG_SUPP, "FigureS2_recombination_events_over_time.pdf"),
       device = cairo_pdf, width = 120, height = 100, units = "mm")
