## Event studies with the FULL pre-intervention period shown monthly.
## Same specification as analysis.R fit_es, but leads are no longer binned at
## <= -25: every pre-treatment month from et = -80 (Jan 2017) to et = -2 is
## its own dummy (reference: et = -1, August 2023). Writes
## output/event_study_full.csv for the six road-user/severity series used in
## Figures 1, 1a and the combined grid.
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

panel <- readRDS("data/panel.rds")
r <- panel[panel$grp == "restricted", ]

fit_es_full <- function(dat, outcome, offset = "log(pop)") {
  et <- dat$event_time
  lab <- ifelse(!dat$wales | et == -1, "eng", paste0("t", et))
  dat$ev <- factor(lab, levels = c("eng",
                 paste0("t", sort(unique(et[et >= 0]))),
                 paste0("t", sort(unique(et[et < -1])))))
  fml <- as.formula(paste(outcome, "~ factor(la) + factor(month_idx) +",
                          "ev + offset(", offset, ")"))
  cat("  poisson... ")
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[outcome]], fitted(p), limit = 100)
  cat("theta", round(th, 1), " nb... ")
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
             control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  se <- sqrt(diag(V)); cf <- coef(fit)
  idx <- grep("^evt", names(cf))
  es <- data.frame(term = sub("^ev", "", names(cf))[idx],
                   beta = cf[idx], se = se[idx], outcome = outcome,
                   stringsAsFactors = FALSE)
  es$et <- as.integer(sub("^t", "", es$term))
  es$lo <- es$beta - 1.96 * es$se; es$hi <- es$beta + 1.96 * es$se
  cat("vcov done\n")
  es
}

out <- list()
for (oc in c("n_all", "n_ksi", "n_ped_ksi", "n_cyc_ksi", "n_ped", "n_cyc")) {
  cat("event study (full pre-period):", oc, "\n")
  out[[oc]] <- fit_es_full(r, oc)
}
es <- do.call(rbind, out)
write.csv(es, "output/event_study_full.csv", row.names = FALSE)
cat("saved output/event_study_full.csv:", nrow(es), "rows;",
    "et range:", range(es$et), "\n")
