data_path <- "data_HF.csv"
out_csv <- "HF_anova.csv"
out_txt <- "HF_anova.txt"
title_text <- "HF monthly total tick density ANOVA analysis"

raw <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
raw$sampling_date <- as.Date(raw[["sampling time"]], format = "%Y/%m/%d")
raw$month_num <- as.integer(format(raw$sampling_date, "%m"))
observed_months <- sort(unique(raw$month_num[!is.na(raw$month_num)]))
raw$month <- factor(month.name[raw$month_num], levels = month.name[observed_months])
raw$density <- as.numeric(trimws(raw[["tick_density_total"]]))

dat <- raw[!is.na(raw$month) & !is.na(raw$density), ]
dat$month <- droplevels(dat$month)
dat$log_density <- if (any(dat$density <= 0)) log1p(dat$density) else log(dat$density)
log_label <- if (any(dat$density <= 0)) "Log1p-transformed total tick density" else "Log-transformed total tick density"

se <- function(x) {
  n <- length(x)
  if (n <= 1) return(NA_real_)
  sd(x) / sqrt(n)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.4f", p)
}

levene_median <- function(y, group) {
  abs_dev <- ave(y, group, FUN = function(x) abs(x - median(x)))
  fit <- aov(abs_dev ~ group)
  tab <- summary(fit)[[1]]
  data.frame(
    df1 = tab["group", "Df"],
    df2 = tab["Residuals", "Df"],
    F = tab["group", "F value"],
    p_value = tab["group", "Pr(>F)"],
    row.names = NULL
  )
}

anova_result <- function(y, group, label) {
  fit <- aov(y ~ group)
  tab <- summary(fit)[[1]]
  shapiro <- shapiro.test(residuals(fit))
  lev <- levene_median(y, group)
  tukey <- as.data.frame(TukeyHSD(fit, "group")$group)
  tukey$comparison <- rownames(tukey)
  rownames(tukey) <- NULL
  names(tukey) <- c("mean_difference", "lower_95_ci", "upper_95_ci", "p_adjusted", "comparison")
  tukey <- tukey[, c("comparison", "mean_difference", "lower_95_ci", "upper_95_ci", "p_adjusted")]
  list(
    label = label,
    fit = fit,
    overall = data.frame(
      analysis = label,
      test = "One-way ANOVA",
      df1 = tab["group", "Df"],
      df2 = tab["Residuals", "Df"],
      F = tab["group", "F value"],
      p_value = tab["group", "Pr(>F)"],
      row.names = NULL
    ),
    shapiro = data.frame(
      analysis = label,
      test = "Shapiro-Wilk test on ANOVA residuals",
      W = unname(shapiro$statistic),
      p_value = shapiro$p.value,
      row.names = NULL
    ),
    levene = data.frame(
      analysis = label,
      test = "Median-centered Levene test",
      df1 = lev$df1,
      df2 = lev$df2,
      F = lev$F,
      p_value = lev$p_value,
      row.names = NULL
    ),
    tukey = cbind(analysis = label, tukey)
  )
}

bind_fill <- function(parts) {
  all_names <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(df) {
    missing <- setdiff(all_names, names(df))
    for (m in missing) df[[m]] <- NA
    df[, all_names]
  })
  do.call(rbind, parts)
}

make_section <- function(section, df) {
  df$section <- section
  df[, c("section", setdiff(names(df), "section"))]
}

desc <- aggregate(
  density ~ month,
  dat,
  function(x) c(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    se = se(x),
    median = median(x),
    min = min(x),
    max = max(x)
  )
)
desc <- do.call(data.frame, desc)
names(desc) <- c("month", "n", "mean", "sd", "se", "median", "min", "max")
desc$section <- "descriptive_statistics"
desc <- desc[, c("section", "month", "n", "mean", "sd", "se", "median", "min", "max")]

raw_analysis <- anova_result(dat$density, dat$month, "Raw total tick density")
log_analysis <- anova_result(dat$log_density, dat$month, log_label)

welch_raw <- tryCatch(
  oneway.test(density ~ month, data = dat, var.equal = FALSE),
  error = function(e) e
)
if (inherits(welch_raw, "error")) {
  welch_note <- welch_raw$message
  welch <- data.frame(
    section = "overall_tests",
    analysis = "Raw total tick density",
    test = "Welch one-way ANOVA",
    df1 = NA_real_,
    df2 = NA_real_,
    F = NA_real_,
    p_value = NA_real_,
    note = welch_note
  )
} else {
  welch_note <- "calculated"
  welch <- data.frame(
    section = "overall_tests",
    analysis = "Raw total tick density",
    test = "Welch one-way ANOVA",
    df1 = unname(welch_raw$parameter[1]),
    df2 = unname(welch_raw$parameter[2]),
    F = unname(welch_raw$statistic),
    p_value = welch_raw$p.value,
    note = welch_note
  )
}

csv_parts <- list(
  data_summary = data.frame(
    section = "data_summary",
    metric = c("source_file", "n_records", "months_compared"),
    value = c(data_path, nrow(dat), paste(levels(dat$month), collapse = "; "))
  ),
  descriptive_statistics = desc,
  assumption_tests = bind_fill(list(
    make_section("assumption_tests", raw_analysis$shapiro),
    make_section("assumption_tests", raw_analysis$levene),
    make_section("assumption_tests", log_analysis$shapiro),
    make_section("assumption_tests", log_analysis$levene)
  )),
  overall_tests = bind_fill(list(
    make_section("overall_tests", raw_analysis$overall),
    welch,
    make_section("overall_tests", log_analysis$overall)
  )),
  tukey_hsd_posthoc = bind_fill(list(
    make_section("tukey_hsd_posthoc", raw_analysis$tukey),
    make_section("tukey_hsd_posthoc", log_analysis$tukey)
  ))
)

csv_out <- bind_fill(csv_parts)
write.csv(csv_out, out_csv, row.names = FALSE)

peak_idx <- which.max(desc$mean)
peak_month <- desc$month[peak_idx]
peak_mean <- desc$mean[peak_idx]
peak_se <- desc$se[peak_idx]

raw_sig <- raw_analysis$overall$p_value < 0.05
log_sig <- log_analysis$overall$p_value < 0.05
raw_tukey_sig <- raw_analysis$tukey[raw_analysis$tukey$p_adjusted < 0.05, ]
log_tukey_sig <- log_analysis$tukey[log_analysis$tukey$p_adjusted < 0.05, ]
welch_text <- if (is.na(welch$p_value)) {
  paste0("Raw-density Welch ANOVA was not calculated: ", welch_note, ".")
} else {
  paste0("Raw-density Welch ANOVA: F(",
         sprintf("%.4f", welch$df1), ", ", sprintf("%.4f", welch$df2), ") = ",
         sprintf("%.4f", welch$F), ", P = ", fmt_p(welch$p_value), ".")
}
transform_note <- if (log_label == "Log1p-transformed total tick density") {
  "Because zero density values were present, log1p transformation [log(1 + density)] was used for the sensitivity analysis."
} else {
  "Because all density values were positive, natural log transformation was used for the sensitivity analysis."
}

txt <- c(
  title_text,
  paste(rep("=", nchar(title_text)), collapse = ""),
  "",
  paste0("Data source: ", data_path),
  "Outcome: tick_density_total (ticks per 100 m)",
  paste0("Months compared: ", paste(levels(dat$month), collapse = ", ")),
  paste0("Number of observations: ", nrow(dat)),
  "",
  "Workflow",
  "--------",
  paste(
    "Each original sampling record was treated as one observation.",
    "Monthly differences in total tick density were first assessed using one-way ANOVA.",
    "Normality was checked using the Shapiro-Wilk test applied to ANOVA residuals, and homogeneity of variances was checked using a median-centered Levene test.",
    "Tukey's HSD test was used for post-hoc multiple comparisons.",
    "A transformed-density ANOVA was additionally performed as a sensitivity analysis.",
    "Welch's ANOVA was also calculated for transparency.",
    transform_note
  ),
  "",
  "Descriptive result",
  "------------------",
  paste0("The highest mean total tick density was observed in ", peak_month,
         " (mean = ", sprintf("%.4f", peak_mean),
         ", SE = ", sprintf("%.4f", peak_se), " ticks per 100 m)."),
  "",
  "Assumption checks",
  "-----------------",
  paste0("Raw-density ANOVA residual normality: Shapiro-Wilk W = ",
         sprintf("%.4f", raw_analysis$shapiro$W), ", P = ", fmt_p(raw_analysis$shapiro$p_value), "."),
  paste0("Raw-density homogeneity of variances: Levene F(",
         raw_analysis$levene$df1, ", ", raw_analysis$levene$df2, ") = ",
         sprintf("%.4f", raw_analysis$levene$F), ", P = ", fmt_p(raw_analysis$levene$p_value), "."),
  paste0(log_label, " residual normality: Shapiro-Wilk W = ",
         sprintf("%.4f", log_analysis$shapiro$W), ", P = ", fmt_p(log_analysis$shapiro$p_value), "."),
  paste0(log_label, " homogeneity of variances: Levene F(",
         log_analysis$levene$df1, ", ", log_analysis$levene$df2, ") = ",
         sprintf("%.4f", log_analysis$levene$F), ", P = ", fmt_p(log_analysis$levene$p_value), "."),
  "",
  "Overall tests",
  "-------------",
  paste0("Raw-density one-way ANOVA: F(",
         raw_analysis$overall$df1, ", ", raw_analysis$overall$df2, ") = ",
         sprintf("%.4f", raw_analysis$overall$F), ", P = ", fmt_p(raw_analysis$overall$p_value), "."),
  welch_text,
  paste0(log_label, " one-way ANOVA: F(",
         log_analysis$overall$df1, ", ", log_analysis$overall$df2, ") = ",
         sprintf("%.4f", log_analysis$overall$F), ", P = ", fmt_p(log_analysis$overall$p_value), "."),
  "",
  "Interpretation",
  "--------------",
  if (raw_sig) {
    paste0("Using the pre-specified raw-density one-way ANOVA, total tick density differed significantly among months (P = ",
           fmt_p(raw_analysis$overall$p_value), ").")
  } else {
    paste0("Using the pre-specified raw-density one-way ANOVA, no statistically significant difference in total tick density was detected among months (P = ",
           fmt_p(raw_analysis$overall$p_value), ").")
  },
  if (log_sig) {
    paste0("The transformed-density sensitivity analysis also indicated significant monthly differences (P = ",
           fmt_p(log_analysis$overall$p_value), ").")
  } else {
    paste0("The transformed-density sensitivity analysis did not indicate significant monthly differences (P = ",
           fmt_p(log_analysis$overall$p_value), ").")
  },
  paste0("Number of significant Tukey HSD pairwise comparisons on raw density: ", nrow(raw_tukey_sig), "."),
  paste0("Number of significant Tukey HSD pairwise comparisons on transformed density: ", nrow(log_tukey_sig), "."),
  "",
  "Suggested Methods text",
  "----------------------",
  paste(
    "Differences in monthly total tick density were assessed using one-way analysis of variance (ANOVA), with each sampling record treated as an independent observation.",
    "Prior to ANOVA, normality was evaluated using the Shapiro-Wilk test applied to model residuals, and homogeneity of variances was assessed using a median-centered Levene test.",
    "Tukey's honestly significant difference (HSD) test was used for post-hoc multiple comparisons.",
    "A transformed-density ANOVA was additionally conducted as a sensitivity analysis, and Welch's ANOVA was calculated to assess robustness to unequal variances."
  ),
  "",
  "Suggested Results text",
  "----------------------",
  paste0("For raw total tick density, the Shapiro-Wilk test yielded W = ",
         sprintf("%.4f", raw_analysis$shapiro$W), ", P = ", fmt_p(raw_analysis$shapiro$p_value),
         ", and Levene's test yielded F(",
         raw_analysis$levene$df1, ", ", raw_analysis$levene$df2, ") = ",
         sprintf("%.4f", raw_analysis$levene$F), ", P = ", fmt_p(raw_analysis$levene$p_value), ". ",
         if (raw_sig) "One-way ANOVA detected a significant difference in monthly total tick density " else "One-way ANOVA did not detect a significant difference in monthly total tick density ",
         "[F(", raw_analysis$overall$df1, ", ", raw_analysis$overall$df2, ") = ",
         sprintf("%.4f", raw_analysis$overall$F), ", P = ", fmt_p(raw_analysis$overall$p_value),
         "]. The highest mean density occurred in ", peak_month,
         " (mean = ", sprintf("%.4f", peak_mean), " ticks per 100 m)."),
  "",
  paste0("Full numerical outputs are provided in ", out_csv, ".")
)

writeLines(txt, out_txt, useBytes = TRUE)
