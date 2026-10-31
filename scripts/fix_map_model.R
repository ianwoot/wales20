## Corrected LA-level (map) model.
##
## Two fixes vs the original section 5 of analysis.R:
##
## 1. Outcome: severity-ADJUSTED KSI. The Dyfed-Powys force switched to
##    injury-based reporting (CRASH/IBRS) in Q2 2023, just before the policy.
##    Carmarthenshire, Pembrokeshire, Ceredigion and Powys sit in that force,
##    so their REPORTED KSI series jump at the switch and inflate the LA-level
##    post/pre contrast. The adjusted series removes the artefact.
##
## 2. Standard errors: model-based (Hessian) for the LA-specific coefficients.
##    Cluster-robust (LA) variance is degenerate here: each LA's interaction
##    coefficient is identified only within its own cluster, whose score is
##    exactly zero at the MLE, so the sandwich meat collapses and the
##    clustered SEs are an order of magnitude too small (see
##    scripts/diagnose_map_se.R). The common month effects are estimated from
##    296 English authorities and are effectively exact, so the Hessian SE
##    is the right within-authority uncertainty measure.
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich)
}))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
rr <- panel[panel$grp == "restricted", ]

welsh_las <- sort(unique(rr$la[rr$wales]))
D <- sapply(welsh_las, function(x) as.numeric(rr$la == x & rr$post == 1))
colnames(D) <- welsh_las

## two-step NB on adjusted KSI
fml_adj <- adj_ksi ~ factor(la) + factor(month_idx) + offset(log(pop)) + D
p <- glm(fml_adj, data = rr, family = poisson(), control = glm.control(maxit = 300))
th <- MASS::theta.ml(rr$adj_ksi, fitted(p), limit = 100)
fit <- glm(fml_adj, data = rr, family = MASS::negative.binomial(theta = th),
           control = glm.control(maxit = 300))
V <- vcov(fit)                      # model-based (see header)
cf <- coef(fit); se <- sqrt(diag(V))
idx <- grep("^D", names(cf))

map_tab <- data.frame(la = sub("^D", "", names(cf))[idx],
                      beta = cf[idx], se = se[idx],
                      stringsAsFactors = FALSE)
map_tab$outcome <- "adj_ksi"
map_tab$theta <- th
map_tab$irr <- exp(map_tab$beta)
map_tab$ci_lo <- exp(map_tab$beta - 1.96 * map_tab$se)
map_tab$ci_hi <- exp(map_tab$beta + 1.96 * map_tab$se)

## also record the reported-KSI betas for comparison (same spec)
fml_rep <- n_ksi ~ factor(la) + factor(month_idx) + offset(log(pop)) + D
p2 <- glm(fml_rep, data = rr, family = poisson(), control = glm.control(maxit = 300))
th2 <- MASS::theta.ml(rr$n_ksi, fitted(p2), limit = 100)
fit2 <- glm(fml_rep, data = rr, family = MASS::negative.binomial(theta = th2),
            control = glm.control(maxit = 300))
map_tab$beta_reported <- coef(fit2)[idx]

map_tab <- map_tab[order(map_tab$irr), ]
write.csv(map_tab, "output/la_effects.csv", row.names = FALSE)
print(map_tab[, c("la", "irr", "ci_lo", "ci_hi", "beta_reported")], digits = 3,
      row.names = FALSE)

## refresh the map entry in models.rds
if (file.exists("output/models.rds")) {
  mo <- readRDS("output/models.rds")
  mo$map <- map_tab
  saveRDS(mo, "output/models.rds")
}
cat("done\n")
