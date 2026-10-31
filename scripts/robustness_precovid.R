## Robustness: pre-COVID window (2017-2019) vs 2024, restricted roads.

# robustness_precovid.R is a robustness check that re-runs the headline DiD while excluding 
# the COVID-era months from the pre-period. The concern it addresses: the main model's baseline 
# (2017–2025) includes 2020–2022, when lockdowns suppressed casualty counts in both countries. 
# If Wales and England didn't recover from the pandemic in the same way, some of the measured 
# "treatment effect" could actually be differential COVID recovery — so the check tests whether 
# the result survives when the pandemic years are dropped entirely.
#
# How it works:
#  
# • Pre-period: Jan 2017 – Dec 2019 only (month_idx ≤ 36) — a clean, pre-COVID baseline.
# • Post-period: full-year 2024 (month_idx 85–96) — the first complete calendar year after 
# the intervention, also avoiding any transition noise in late 2023.
# • Everything else is the same model as the headline analysis: LA fixed effects, calendar-month
# fixed effects, tp (Wales × post), Wales-specific seasonality, population offset, negative binomial, 
# LA-clustered sandwich SEs. The only thing that changes is which months define "before" and "after."

suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
r <- panel[panel$grp == "restricted", ]
r19 <- r[r$month_idx <= 36 | (r$month_idx >= 85 & r$month_idx <= 96), ]  # 2017-2019 vs full-year 2024 (months 85-96)

add_wales_moy <- function(dat) {
  dat$wales_moy <- factor(ifelse(dat$wales, paste0("wm", dat$moy), "base")); dat }

fit_nb <- function(fml, dat) {
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  list(fit = fit, se = sqrt(diag(V)))
}

out <- do.call(rbind, lapply(c("n_all", "n_ksi", "adj_ksi", "n_ped_ksi", "n_cyc_ksi"),
                             function(oc) {
  d <- add_wales_moy(r19)
  fml <- as.formula(paste(oc, "~ factor(la) + factor(month_idx) + tp + wales_moy + offset(log(pop))"))
  z <- fit_nb(fml, d)
  cf <- coef(z$fit); b <- cf["tp"]; s <- z$se["tp"]
  data.frame(outcome = oc, beta = b, se = s, irr = exp(b),
             ci_lo = exp(b - 1.96 * s), ci_hi = exp(b + 1.96 * s))
}))
write.csv(out, "output/did_precovid_2024.csv", row.names = FALSE)
print(out, digits = 3)
cat("done\n")
