## Controlled interrupted time-series (ITS) decomposition of the DiD.
##
## Replaces the event-study dummies with a parametric Wales time profile:
##   wales_t : Wales x linear trend (whole window)      -> pre-period slope
##   tp      : Wales x post                              -> LEVEL shift at Sep 2023
##   tpx     : Wales x (month_idx - 81)+                 -> SLOPE change after
## with LA FE, calendar-month FE, Wales-specific seasonality, population
## offset, NB (two-step) errors and LA-clustered SEs.
##
## Interpretation: tp answers "did the lines shift at the intervention?";
## tpx answers "did the slope change after it?"; wales_t answers "was there
## a differential pre-trend at all?" (null expected from flat annual ratios).
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
r <- panel[panel$grp == "restricted", ]
r$wales_moy <- factor(ifelse(r$wales, paste0("wm", r$moy), "base"))
r$wales_t  <- as.numeric(r$wales) * r$month_idx
r$tpx      <- as.numeric(r$wales) * pmax(r$month_idx - 81, 0)

fit_nb <- function(fml, dat) {
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  list(fit = fit, se = sqrt(diag(V)), theta = th)
}

out <- do.call(rbind, lapply(c("n_all", "n_ksi", "adj_ksi", "n_ped_ksi", "n_cyc_ksi",
                              "n_ped", "n_cyc"),
                             function(oc) {
  cat("ITS:", oc, "... ")
  fml <- as.formula(paste(oc,
    "~ factor(la) + factor(month_idx) + wales_moy + wales_t + tp + tpx + offset(log(pop))"))
  z <- fit_nb(fml, r)
  cf <- coef(z$fit); se <- z$se
  mk <- function(v) {
    b <- cf[v]; s <- se[v]
    data.frame(term = v, beta = b, se = s, z = b / s,
               irr = exp(b), ci_lo = exp(b - 1.96 * s), ci_hi = exp(b + 1.96 * s))
  }
  res <- do.call(rbind, lapply(c("wales_t", "tp", "tpx"), mk))
  res$outcome <- oc; res$theta <- z$theta
  cat("level IRR", round(exp(cf["tp"]), 3), " slope-change z", round((cf["tpx"]/se["tpx"]), 2), "\n")
  res
}))
write.csv(out, "output/its_decomposition.csv", row.names = FALSE)
cat("\n")
print(out[, c("outcome", "term", "irr", "ci_lo", "ci_hi", "z")], digits = 3, row.names = FALSE)
cat("done\n")
