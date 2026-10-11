## ------------------------------------------------------------------------
## Analysis of hetero_ddf_sim.R results: Type-I error rates (alpha = 0.05)
## with binom.test(), and a 6-panel ECDF figure (one panel per method),
## mirroring the figure at
## https://lgraz.com/posts/lmm-heteroskedastic/#type-i-error-rate-on-heteroscedastic-data
##
## Uses the "anova" test column (the omnibus test for the trt effect --
## Satterthwaite/Kenward-Roger F-test, LRT, or Wald chisq depending on
## method), since that's the effect the ddf corrections in hetero_ddf_sim.R
## are meant to fix. Under the simulated null (no trt effect), well-calibrated
## p-values should be ~Uniform(0,1), i.e. their ECDF should track the diagonal;
## methods with anti-conservative (inflated Type-I error) p-values bow above
## the diagonal near small p.
## ------------------------------------------------------------------------

library(ggplot2)

alpha <- 0.05

results_long <- readRDS("hetero_ddf_sim_results.rds")

## ---- Type-I error table (both "anova" and "contrast" tests) ----
tab <- aggregate(p_value ~ method + test, data = results_long,
                  FUN = function(p) mean(p < alpha, na.rm = TRUE))
names(tab)[3] <- "type1_error"
tab$n_na <- aggregate(p_value ~ method + test, data = results_long,
                       FUN = function(p) sum(is.na(p)))$p_value
tab$binom_p <- mapply(function(m, t) {
    p <- results_long$p_value[results_long$method == m & results_long$test == t]
    p <- p[!is.na(p)]
    binom.test(sum(p < alpha), length(p), alpha)$p.value
}, tab$method, tab$test)
tab <- tab[order(tab$test, tab$method), ]

cat(sprintf("---- Type-I error rate at alpha = %.2f (binom.test against nominal) ----\n", alpha))
print(tab, row.names = FALSE)

## ---- 6-panel ECDF figure (one panel per method, "anova" test) ----
method_order <- c("lmerTest", "lme_varIdent", "glmmTMB_none", "glmmTMB_wald",
                   "glmmTMB_satterthwaite", "glmmTMB_kenward-roger")

panel_data <- results_long[results_long$test == "anova", ]
panel_data <- panel_data[!is.na(panel_data$p_value), ]
panel_data$method <- factor(panel_data$method, levels = method_order)

rate_labs <- aggregate(p_value ~ method, data = panel_data,
                        FUN = function(p) mean(p < alpha))
names(rate_labs)[2] <- "rate"
rate_labs$label <- sprintf("Type I error = %.3f", rate_labs$rate)

gg <- ggplot(panel_data, aes(p_value)) +
    stat_ecdf(geom = "step", colour = "steelblue", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = alpha, linetype = "dashed", colour = "red") +
    geom_text(data = rate_labs, aes(x = 1, y = 0, label = label),
              hjust = 1, vjust = -0.5, size = 3) +
    facet_wrap(~ method, nrow = 2) +
    scale_x_sqrt(limits = c(0, 1), breaks = c(0, alpha, 0.25, 0.5, 0.75, 1)) +
    scale_y_sqrt(limits = c(0, 1)) +
    labs(x = "p-value", y = "ECDF") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

ggsave("hetero_ddf_sim_ecdf.png", gg, width = 9, height = 6, dpi = 150)

cat("\nFigure written to hetero_ddf_sim_ecdf.png\n")
