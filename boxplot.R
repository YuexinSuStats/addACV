# Draw boxplots for S1-L-n300 Monte Carlo estimates.
#
# Place this script at:
#   Rcode_med/boxplot.R
#
# Expected data location (relative to this script):
#   Scenario2/L-n300/S1-L-n300_rep1-500_summary.csv
#
# A boxplot requires the replication-level estimates. The summary file is
# located and read first; if it contains aggregate statistics only, the script
# automatically reads S1-L-n300_rep1-500.csv from the same folder.

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing package(s): ", paste(missing_packages, collapse = ", "),
    ". Install with: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

library(ggplot2)

# Find the folder containing boxplot.R. This works with Rscript, source(), and
# the RStudio Source button; no user-specific absolute path is embedded.
get_script_dir <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", all_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }

  frames <- sys.frames()
  sourced_paths <- vapply(
    frames,
    function(frame) {
      if (is.null(frame$ofile)) "" else as.character(frame$ofile)
    },
    character(1)
  )
  sourced_paths <- sourced_paths[nzchar(sourced_paths)]
  if (length(sourced_paths) > 0L) {
    script_path <- sourced_paths[[length(sourced_paths)]]
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    script_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    if (nzchar(script_path)) {
      return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
summary_name <- "S1-L-n300_rep1-500_summary.csv"

# First check the required relative location. The lower-case alternative keeps
# the script tolerant of folder-name capitalization differences.
summary_candidates <- file.path(
  script_dir,
  c("Scenario2", "scenario2"),
  "L-n300",
  summary_name
)
existing_exact <- summary_candidates[file.exists(summary_candidates)]
summary_hits <- if (length(existing_exact) > 0L) {
  existing_exact[[1L]]
} else {
  character(0)
}

# If the exact relative path is unavailable, search below Rcode_med by name.
if (length(summary_hits) == 0L) {
  summary_hits <- list.files(
    path = script_dir,
    pattern = "^S1-L-n300_rep1-500_summary[.]csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  normalized_hits <- gsub("\\\\", "/", summary_hits)
  preferred <- grepl(
    "/scenario2/l-n300/",
    tolower(normalized_hits),
    fixed = TRUE
  )
  if (any(preferred)) {
    summary_hits <- summary_hits[preferred]
  }
}

if (length(summary_hits) == 0L) {
  stop(
    "Could not find ",
    file.path("Scenario2", "L-n300", summary_name),
    " relative to boxplot.R.",
    call. = FALSE
  )
}

if (length(summary_hits) > 1L) {
  stop(
    "Multiple matching summary files were found:\n",
    paste(summary_hits, collapse = "\n"),
    call. = FALSE
  )
}

summary_csv <- normalizePath(
  summary_hits[[1L]],
  winslash = "/",
  mustWork = TRUE
)

summary_results <- read.csv(
  summary_csv,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

method_order <- c("WEE", "ICV1", "ICV2", "ICV4")
parameter_order <- c("beta1", "beta2")
estimate_columns <- as.vector(
  outer(parameter_order, method_order, paste, sep = ".")
)

# Standard summary output contains Bias/SD/SE rather than the 500 individual
# estimates, so it cannot itself define boxplot quartiles. Use the companion
# replication-level file when those estimate columns are absent.
if (all(c("rep", estimate_columns) %in% names(summary_results))) {
  input_csv <- summary_csv
  results <- summary_results
} else {
  raw_name <- sub("_summary[.]csv$", ".csv", basename(summary_csv), ignore.case = TRUE)
  raw_csv <- file.path(dirname(summary_csv), raw_name)

  if (!file.exists(raw_csv)) {
    stop(
      "The summary CSV contains aggregate statistics only. A boxplot needs ",
      "replication-level estimates. Place ", raw_name,
      " in the same folder as ", basename(summary_csv), ".",
      call. = FALSE
    )
  }

  input_csv <- normalizePath(raw_csv, winslash = "/", mustWork = TRUE)
  results <- read.csv(
    input_csv,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

missing_columns <- setdiff(c("rep", estimate_columns), names(results))
if (length(missing_columns) > 0L) {
  stop(
    "Required column(s) missing from the CSV: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# Convert the wide result table to one row per replication, method, and
# parameter. If an ICV validity flag exists, invalid fits are excluded.
long_parts <- lapply(method_order, function(method) {
  valid_column <- paste0("valid.", method)
  valid_fit <- if (valid_column %in% names(results)) {
    !is.na(results[[valid_column]]) & as.logical(results[[valid_column]])
  } else {
    rep(TRUE, nrow(results))
  }

  do.call(
    rbind,
    lapply(parameter_order, function(parameter) {
      estimate_column <- paste(parameter, method, sep = ".")
      data.frame(
        rep = results$rep,
        Method = method,
        Parameter = parameter,
        Estimate = suppressWarnings(as.numeric(results[[estimate_column]])),
        Valid = valid_fit,
        stringsAsFactors = FALSE
      )
    })
  )
})

plot_data <- do.call(rbind, long_parts)
plot_data <- plot_data[
  plot_data$Valid & is.finite(plot_data$Estimate),
  c("rep", "Method", "Parameter", "Estimate"),
  drop = FALSE
]

if (nrow(plot_data) == 0L) {
  stop("No finite, valid coefficient estimates are available to plot.", call. = FALSE)
}

plot_data$Method <- factor(plot_data$Method, levels = method_order)
plot_data$Parameter <- factor(
  plot_data$Parameter,
  levels = parameter_order,
  labels = c("hat(beta)[1]", "hat(beta)[2]")
)

true_beta <- 1

p <- ggplot(plot_data, aes(x = Method, y = Estimate)) +
  geom_hline(
    yintercept = true_beta,
    linetype = "dashed",
    color = "#B14C2E",
    linewidth = 0.65
  ) +
  geom_boxplot(
    width = 0.62,
    fill = "#6BAED6",
    color = "#252525",
    linewidth = 0.45,
    outlier.shape = 16,
    outlier.size = 0.8,
    outlier.alpha = 0.45,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ Parameter,
    nrow = 1,
    scales = "free_y",
    labeller = label_parsed
  ) +
  labs(
    x = "Method",
    y = "Estimate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#252525"),
    panel.spacing = grid::unit(1.1, "lines"),
    plot.margin = margin(6, 12, 8, 10)
  )

output_stem <- paste0(tools::file_path_sans_ext(basename(input_csv)), "_boxplot")
png_file <- file.path(script_dir, paste0(output_stem, ".png"))
pdf_file <- file.path(script_dir, paste0(output_stem, ".pdf"))

ggsave(
  filename = png_file,
  plot = p,
  width = 8.5,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 8.5,
  height = 4.8,
  units = "in",
  device = "pdf"
)

print(p)

cat("Valid estimates per panel and method:\n")
print(with(plot_data, table(Parameter, Method)))
cat("\nSaved:\n", png_file, "\n", pdf_file, "\n", sep = "")
