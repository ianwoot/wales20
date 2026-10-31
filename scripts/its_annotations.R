## ITS quantities for the annotated event-study figure.
## For each of the six road-user/severity outcomes (restricted roads, pop
## offset), fit the ITS model (same spec as its_decomposition.R) and extract
## with LA-clustered CIs:
##   pre slope  : gamma (wales_t)                          - baseline trend
##   step       : beta1 (tp)                              - level shift at Sep 2023
##   slope chg  : beta2 (tpx)                             - bend after Sep 2023
##   post slope : gamma + beta2 (with cov-based CI)        - the trend actually
##                                                          followed post-policy
## Output: output/its_annotations.csv
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
r <- panel[panel$grp == "restricted", ]
r$wales_moy <- factor(ifelse(r$wales, paste0("wm", r$moy), "base"))
r$wales_t <- as.numeric(r$wales) * r$month_idx
r$tpx <- as.numeric(r$wales) * pmax(r$month_idx - 81, 0)

fit_nb <- function(fml, dat) {
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  list(fit = fit, V = V, theta = th)
}

out <- do.call(rbind, lapply(c("n_all", "n_ksi", "n_ped_ksi", "n_cyc_ksi",
                               "n_ped", "n_cyc"), function(oc) {
  cat("ITS:", oc, "...\n")
  fml <- as.formula(paste(oc,
    "~ factor(la) + factor(month_idx) + wales_moy + wales_t + tp + tpx + offset(log(pop))"))
  z <- fit_nb(fml, r)
  cf <- coef(z$fit); V <- z$V
  b <- cf[c("wales_t", "tp", "tpx")]
  se <- sqrt(diag(V))[c("wales_t", "tp", "tpx")]
  # post slope = gamma + beta2, CI from the 2x2 covariance block
  k <- c("wales_t", "tpx")
  post <- sum(b[k])
  post_se <- sqrt(t(c(1, 1)) %*% V[k, k] %*% c(1, 1))
  data.frame(
    outcome = oc,
    pre_slope = b["wales_t"], pre_lo = b["wales_t"] - 1.96 * se["wales_t"],
    pre_hi = b["wales_t"] + 1.96 * se["wales_t"],
    step = b["tp"], step_lo = b["tp"] - 1.96 * se["tp"],
    step_hi = b["tp"] + 1.96 * se["tp"],
    slope_chg = b["tpx"], chg_lo = b["tpx"] - 1.96 * se["tpx"],
    chg_hi = b["tpx"] + 1.96 * se["tpx"],
    post_slope = post, post_lo = post - 1.96 * post_se,
    post_hi = post + 1.96 * post_se)
}))
write.csv(out, "output/its_annotations.csv", row.names = FALSE)
print(out, digits = 4, row.names = FALSE)
cat("done\n")
