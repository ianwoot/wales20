## Difference-in-differences analysis: Wales default 20mph (Sep 2023).
## Base R + MASS::negative.binomial + sandwich/lmtest for LA-clustered SEs.

suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious

## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
## Differential seasonality: Wales x month-of-year dummies only
add_wales_moy <- function(dat) {
  dat$wales_moy <- factor(ifelse(dat$wales, paste0("wm", dat$moy), "base"))
  dat
}

## Two-step negative binomial: (1) Poisson fit, (2) dispersion theta via
## theta.ml on fitted means, (3) one NB fit with fixed theta.
## Cluster-robust (LA) SEs; sandwich handles any residual overdispersion.
fit_nb <- function(fml, dat) {
  tt <- Sys.time()
  step <- function(s) cat(s, round(as.numeric(difftime(Sys.time(), tt, units = "secs")), 1), "s | ")
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  step("pois")
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  step("theta")
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  step("nb")
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  step("vcov\n")
  list(fit = fit, se = sqrt(diag(V)), theta = th)
}

## Static DiD: LA FE + calendar-month FE + Wales seasonality + offset
fit_did <- function(dat, outcome, offset = "log(pop)", extra = "tp") {
  dat <- add_wales_moy(dat)
  fml <- as.formula(paste(outcome, "~ factor(la) + factor(month_idx) +",
                          extra, "+ wales_moy + offset(", offset, ")"))
  z <- fit_nb(fml, dat)
  cf <- coef(z$fit); idx <- match("tp", names(cf))
  data.frame(outcome = outcome, beta = cf[idx], se = z$se[idx],
             irr = exp(cf[idx]), ci_lo = exp(cf[idx] - 1.96 * z$se[idx]),
             ci_hi = exp(cf[idx] + 1.96 * z$se[idx]), theta = z$theta,
             stringsAsFactors = FALSE)
}

## Event study: Wales x event-month dummies, reference = Aug 2023 (et = -1)
fit_es <- function(dat, outcome, offset = "log(pop)") {
  et <- dat$event_time
  lab <- ifelse(!dat$wales | et == -1, "eng",
         ifelse(et <= -25, "tle-25", paste0("t", et)))
  dat$ev <- factor(lab, levels = c("eng", "tle-25",
                 paste0("t", sort(unique(et[et >= 0]))),
                 paste0("t", sort(unique(et[et < -1 & et > -25])))))
  fml <- as.formula(paste(outcome, "~ factor(la) + factor(month_idx) +",
                          "ev + offset(", offset, ")"))
  z <- fit_nb(fml, dat)
  cf <- coef(z$fit); idx <- grep("^ev(t|le-25)", names(cf))
  data.frame(term = sub("^ev", "", names(cf))[idx],
             beta = cf[idx], se = z$se[idx], outcome = outcome,
             stringsAsFactors = FALSE)
}

r <- panel[panel$grp == "restricted", ]
h <- panel[panel$grp == "higher", ]
cat("restricted panel:", nrow(r), "rows; total KSI:", sum(r$n_ksi), "\n")

## ------------------------------------------------------------------
## 1. Headline DiD models (restricted roads: 20/30 mph network)
## ------------------------------------------------------------------
specs <- list(
  c("n_all",       "All casualties"),
  c("n_ksi",      "KSI (reported)"),
  c("adj_ksi",    "KSI (severity-adjusted)"),
  c("n_fatal",    "Fatal"),
  c("n_slight",   "Slight"),
  c("n_ped_ksi",  "Pedestrian KSI"),
  c("n_ped",      "Pedestrian casualties"),
  c("n_cyc_ksi",  "Cyclist KSI"),
  c("n_cyc",      "Cyclist casualties"),
  c("n_child_ksi","Child (0-15) KSI"),
  c("n_older_ksi","Older (65+) KSI"))

t0 <- Sys.time()
tab_main <- do.call(rbind, lapply(seq_along(specs), function(i) {
  s <- specs[[i]]
  cat("fitting", s[1], "-", s[2], "... ")
  out <- fit_did(r, s[1]); out$label <- s[2]; out$offset <- "population"; out
  cat("IRR", round(out$irr, 3), "\n"); out }))
tab_main$spec <- "pop offset"
cat("main table done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

specs_vkm <- list(c("n_all","All casualties"), c("n_ksi","KSI (reported)"),
                  c("adj_ksi","KSI (severity-adjusted)"),
                  c("n_ped_ksi","Pedestrian KSI"))
# Sensitivity check using an alternative offset vehicle km rather than population. 
# Didn't use in the paper - very similar to population offset - but left here for interest
tab_vkm <- do.call(rbind, lapply(specs_vkm, function(s) {
  out <- fit_did(r, s[1], "log(vkm_m_km * 1e6)")
  out$label <- s[2]; out$offset <- "vehicle-km"; out }))
tab_vkm$spec <- "vkm offset"

## Poisson with robust SEs as estimator check (headline)
rr <- add_wales_moy(r)
pf <- glm(n_ksi ~ factor(la) + factor(month_idx) + tp + wales_moy + offset(log(pop)),
          data = rr, family = poisson(), control = glm.control(maxit = 300))
pse <- sqrt(diag(vcovCL(pf, cluster = rr$la, type = "HC1")))["tp"]
pois_check <- data.frame(spec = "Poisson robust", label = "KSI (reported)",
                         beta = coef(pf)["tp"], se = pse,
                         irr = exp(coef(pf)["tp"]),
                         ci_lo = exp(coef(pf)["tp"] - 1.96 * pse),
                         ci_hi = exp(coef(pf)["tp"] + 1.96 * pse),
                         theta = NA_real_, offset = "population",
                         outcome = "n_ksi")

## ------------------------------------------------------------------
## 2. Placebo: 40+ mph roads (no treatment expected)
## ------------------------------------------------------------------
tab_plac <- do.call(rbind, lapply(list(
  c("n_all", "All casualties"), c("n_ksi", "KSI (reported)"),
  c("adj_ksi", "KSI (severity-adjusted)"),
  c("n_ped_ksi", "Pedestrian KSI")), function(s) {
  out <- fit_did(h, s[1]); out$label <- s[2]; out$offset <- "population"; out }))
tab_plac$spec <- "placebo 40+ mph"

tab <- rbind(tab_main, tab_vkm, pois_check, tab_plac)
write.csv(tab, "output/did_main.csv", row.names = FALSE)
cat("\n=== Headline DiD (IRR, 95% CI) ===\n")
print(tab[, c("label","spec","irr","ci_lo","ci_hi")], digits = 3)

# Not unpacked in the paper owing to space, but if you examine this table you'll find signs of a
# pedestrian KSI reduction on 40+ roads. This is based on only around 350 pedestrians across the UK and 
# it's not clear whether or not it is a solid effect.

## ------------------------------------------------------------------
## 3. Event study
## ------------------------------------------------------------------
t1 <- Sys.time()
es_list <- list()
for (oc in c("n_ksi", "adj_ksi", "n_ped_ksi", "n_cyc_ksi", "n_all")) {
  es_list[[oc]] <- fit_es(r, oc)
  cat("event study done:", oc, " (", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min)\n", sep = "")
}
es <- do.call(rbind, es_list)
es$et <- ifelse(es$term == "tle-25", -25, as.integer(sub("^t", "", es$term)))
es$lo <- es$beta - 1.96 * es$se; es$hi <- es$beta + 1.96 * es$se
write.csv(es, "output/event_study.csv", row.names = FALSE)

## ------------------------------------------------------------------
## 4. Permutation placebo: 1,000 random pseudo-Wales groups of 22
##    English LAs; aggregate two-group Poisson with month FE.
##
## NB this is KSIs only. There's also a separate version for all the casualties
## in permutation_all.R
## ------------------------------------------------------------------
set.seed(20250917)
eng_la <- sort(unique(r$la[!r$wales]))
month_tab <- sort(unique(r$month_idx))

## monthly population denominator per LA
pop_m <- r[, c("la", "month_idx", "pop")]
pop_m <- pop_m[!duplicated(pop_m), ]
y_m <- r[, c("la", "month_idx", "n_ksi")]

## precompute LA-month matrices for speed
la_vec <- sort(unique(r$la))
Y <- tapply(y_m$n_ksi, list(y_m$la, y_m$month_idx), sum)
P <- tapply(pop_m$pop, list(pop_m$la, pop_m$month_idx), function(z) z[1])
Y[is.na(Y)] <- 0
eng_idx <- match(eng_la, la_vec)
wls_idx <- grep("^W06", la_vec)

post_vec <- as.integer(month_tab >= 81)  # Sep 2023 = month_idx 81
X_post <- cbind(1, post_vec)
agg_fit <- function(treated_idx, control_idx) {
  yt <- colSums(Y[treated_idx, , drop = FALSE])
  yc <- colSums(Y[control_idx, , drop = FALSE])
  pt <- colSums(P[treated_idx, , drop = FALSE])
  pc <- colSums(P[control_idx, , drop = FALSE])
  nm <- length(month_tab)
  g <- rep(c(1, 0), each = nm)
  post2 <- rep(post_vec, 2)
  y <- c(yt, yc)
  off <- log(c(pt, pc))
  mm <- cbind(1, g, post2, g * post2,
              sapply(month_tab[-1],
                     function(m) rep(as.integer(month_tab == m), 2)))
  dimnames(mm) <- NULL
  storage.mode(mm) <- "double"
  fit <- glm.fit(x = mm, y = y, family = poisson(), offset = off)
  fit$coefficients[4]
}
perm <- numeric(1000)
t2 <- Sys.time()
for (i in 1:1000) {
  tr <- sample(eng_idx, 22)
  perm[i] <- agg_fit(tr, setdiff(eng_idx, tr))
}
cat("permutation done in", round(as.numeric(difftime(Sys.time(), t2, units = "secs")), 1), "s\n")

wales_beta_agg <- agg_fit(wls_idx, eng_idx)
perm_out <- data.frame(
  wales_agg_beta = wales_beta_agg,
  perm_mean = mean(perm),
  perm_p = mean(abs(perm) >= abs(wales_beta_agg)))
write.csv(data.frame(draw = 1:1000, beta = perm), "output/perm_draws.csv", row.names = FALSE)
cat("\nPermutation: Wales agg beta =", round(wales_beta_agg, 3),
    " mean placebo =", round(mean(perm), 3),
    " two-sided perm p =", perm_out$perm_p, "\n")

## ------------------------------------------------------------------
## 5. LA-level treatment effects (for the map)
## ------------------------------------------------------------------
welsh_las <- sort(unique(rr$la[rr$wales]))
D <- sapply(welsh_las, function(x) as.numeric(rr$la == x & rr$post == 1))
colnames(D) <- welsh_las
fml_la <- n_ksi ~ factor(la) + factor(month_idx) + offset(log(pop)) + D
z <- fit_nb(fml_la, rr)
cf_la <- coef(z$fit)
idx <- grep("^D", names(cf_la))
map_tab <- data.frame(
  la = sub("^D", "", names(cf_la))[idx],
  beta = cf_la[idx], se = z$se[idx], stringsAsFactors = FALSE)
map_tab$irr <- exp(map_tab$beta)
map_tab$ci_lo <- exp(map_tab$beta - 1.96 * map_tab$se)
map_tab$ci_hi <- exp(map_tab$beta + 1.96 * map_tab$se)
write.csv(map_tab, "output/la_effects.csv", row.names = FALSE)

saveRDS(list(tab = tab, es = es, perm = perm, wales_agg = wales_beta_agg,
             map = map_tab, pois_check = pois_check),
        "output/models.rds")
cat("\nanalysis complete\n")
