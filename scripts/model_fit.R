## Overall model-fit metrics for the main DiD analyses.
## For each headline outcome (restricted roads, population offset) and the
## 40+ mph placebo contrast, fit:
##   full : outcome ~ factor(la) + factor(month_idx) + wales_moy + tp + offset(log pop)   [NB]
##   null : same without tp                                                                  [NB]
##   pois : full with Poisson                                                               [for AIC]
## and report: logLik, McFadden pseudo-R2 vs null, LR test of tp (df=1),
## Cameron-Trivedi count R2 = cor(y, fitted)^2, Pearson dispersion vs the
## NB variance (should be ~1 if the NB fits), AIC for nb/pois/null.
suppressWarnings(suppressMessages(library(MASS)))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious

fit2 <- function(fml, dat, family) {
  if (family == "nb") {
    p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
    th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
    glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
        control = glm.control(maxit = 300))
  } else glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
}

run_one <- function(oc, dat) {
  dat$wales_moy <- factor(ifelse(dat$wales, paste0("wm", dat$moy), "base"))
  ffull <- as.formula(paste(oc, "~ factor(la) + factor(month_idx) + wales_moy + tp + offset(log(pop))"))
  fnull <- as.formula(paste(oc, "~ factor(la) + factor(month_idx) + wales_moy + offset(log(pop))"))
  cat("fitting", oc, "...\n")
  p_null <- glm(fnull, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th_null <- MASS::theta.ml(dat[[oc]], fitted(p_null), limit = 100)
  m_null <- glm(fnull, data = dat, family = MASS::negative.binomial(theta = th_null),
                control = glm.control(maxit = 300))
  m_pois <- glm(ffull, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th1 <- MASS::theta.ml(dat[[oc]], fitted(m_pois), limit = 100)
  m_full <- glm(ffull, data = dat, family = MASS::negative.binomial(theta = th1),
                control = glm.control(maxit = 300))
  y <- dat[[oc]]; mu <- fitted(m_full)
  v <- mu + mu^2 / th1
  data.frame(
    outcome = oc,
    ll_full = as.numeric(logLik(m_full)), ll_null = as.numeric(logLik(m_null)),
    mcfadden_r2 = 1 - as.numeric(logLik(m_full)) / as.numeric(logLik(m_null)),
    LR = 2 * (as.numeric(logLik(m_full)) - as.numeric(logLik(m_null))),
    LR_p = pchisq(2 * (as.numeric(logLik(m_full)) - as.numeric(logLik(m_null))), 1,
                  lower.tail = FALSE),
    ct_r2 = cor(y, mu)^2,
    pearson_disp = sum((y - mu)^2 / v) / m_full$df.residual,
    theta = th1,
    aic_nb = AIC(m_full), aic_pois = AIC(m_pois), aic_null = AIC(m_null))
}

r <- panel[panel$grp == "restricted", ]
out <- do.call(rbind, lapply(c("n_all", "n_ksi", "adj_ksi", "n_fatal", "n_slight",
                              "n_ped_ksi", "n_ped", "n_cyc_ksi", "n_cyc",
                              "n_child_ksi", "n_older_ksi"), run_one, dat = r))
out$spec <- "headline (restricted roads, pop offset)"
h <- panel[panel$grp == "higher", ]
pl <- run_one("n_all", h)
pl$spec <- "placebo (40+ mph roads)"
out <- rbind(out, pl)
write.csv(out, "output/model_fit.csv", row.names = FALSE)
print(out, digits = 4, row.names = FALSE)
cat("done\n")
