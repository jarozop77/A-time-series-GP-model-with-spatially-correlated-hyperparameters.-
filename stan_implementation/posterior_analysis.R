# Purpose: summarize an existing Stan fit without rerunning the model.
# Inputs: fit.rds and analysis_data.rds. Outputs: summaries and diagnostics under outputs/.

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args || identical(Sys.getenv("REPRO_TEST"), "1")
mode_label <- if (test_mode) "test" else "full"
output_dir <- file.path("outputs", "standard_implementation", mode_label)
fit_file <- file.path(output_dir, "fit.rds")
if (!file.exists(fit_file)) stop("Required fit not found: ", fit_file, ". Run standard_implementation/run_model.R", if (test_mode) " --test" else "", " first.")
required_packages <- c("posterior", "bayesplot", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

fit <- readRDS(fit_file)
draws <- fit$draws()
summary_table <- posterior::summarise_draws(draws, "mean", "sd", ~posterior::quantile2(.x, probs = c(.025, .5, .975)), "rhat", "ess_bulk", "ess_tail")
utils::write.csv(as.data.frame(summary_table), file.path(output_dir, "posterior_summary.csv"), row.names = FALSE)
diagnostics <- fit$diagnostic_summary()
utils::write.csv(as.data.frame(diagnostics), file.path(output_dir, "sampler_diagnostics.csv"), row.names = FALSE)
available <- posterior::variables(draws)
trace_vars <- intersect(c("mu_log_phi", "gamma_phi_sigma", "sd_log_phi", "sigma_icar"), available)
if (length(trace_vars)) {
  p <- bayesplot::mcmc_trace(draws, pars = trace_vars)
  ggplot2::ggsave(file.path(output_dir, "trace_global_parameters.png"), p, width = 8, height = 5, dpi = 150)
}
message("Stan posterior summaries saved in ", output_dir)
