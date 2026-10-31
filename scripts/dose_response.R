## Dose-response analysis: does the effect scale with treatment intensity?
##
## Intensity measure: for each Welsh authority, the pre-period share of its
## restricted-road collisions that occurred on 30-mph-recorded roads. Roads
## recorded at 30 mph before September 2023 and 20 mph after ARE the
## reclassified roads, so this share approximates the fraction of the
## restricted network that the policy touched (casualty-weighted). A Welsh
## authority whose restricted casualties were already mostly on 20 mph roads
## (dense pre-existing zones) received a smaller dose; an English authority
## received none by construction (its 30s stayed 30).
##
## Tests (same DiD machinery: NB two-step, LA + month FE, Wales seasonality,
## LA-clustered SEs, population offset):
##   M1: tp only                     - pooled effect (headline form)
##   M2: wint = tp x intensity only   - fully dose-scaled effect
##   M3: tp + wint                   - unconstrained: residual pooled effect
##                                     at zero intensity plus dose-response
## Plus a two-stage check: WLS of the LA-level adjusted-KSI effects
## (output/la_effects.csv) on intensity across the 22 Welsh authorities.
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

## ------------------------------------------------------------------
## 1. Intensity from the collision file (pre-period 30-share)
## ------------------------------------------------------------------
raw <- "data/raw"
xwalk <- c(
  "E06000028" = "E06000058", "E06000029" = "E06000058", "E07000048" = "E06000058",
  "E07000049" = "E06000059", "E07000050" = "E06000059", "E07000051" = "E06000059",
  "E07000052" = "E06000059", "E07000053" = "E06000059",
  "E07000190" = "E06000066", "E07000191" = "E06000066", "E07000246" = "E06000066",
  "E07000004" = "E06000060", "E07000005" = "E06000060",
  "E07000006" = "E06000060", "E07000007" = "E06000060",
  "E07000201" = "E07000245", "E07000204" = "E07000245",
  "E07000205" = "E07000244", "E07000206" = "E07000244",
  "E07000150" = "E06000061", "E07000152" = "E06000061",
  "E07000153" = "E06000061", "E07000156" = "E06000061",
  "E07000151" = "E06000062", "E07000154" = "E06000062", "E07000155" = "E06000062",
  "E07000026" = "E06000063", "E07000028" = "E06000063", "E07000029" = "E06000063",
  "E07000027" = "E06000064", "E07000030" = "E06000064", "E07000031" = "E06000064",
  "E07000187" = "E06000066", "E07000188" = "E06000066", "E07000189" = "E06000066",
  "E07000163" = "E06000065", "E07000164" = "E06000065", "E07000165" = "E06000065",
  "E07000166" = "E06000065", "E07000167" = "E06000065", "E07000168" = "E06000065",
  "E07000169" = "E06000065", "EHEATHROW" = "E09000017")

col <- read.csv(file.path(raw, "stats19/collisions_2017_2025.csv"),
                colClasses = c("character", "integer", "character", "character",
                               "integer", "integer", "integer", "integer"))
col$la <- unname(xwalk[col$la_ons]); col$la[is.na(col$la)] <- col$la_ons[is.na(col$la)]
col <- col[col$la != "-1", ]
col <- col[grepl("^E|^W", col$la), ]
col$month_idx <- (col$collision_year - 2017) * 12 +
                 as.integer(substr(col$date, 4, 5))

pre <- col[col$month_idx < 81 & col$speed_limit %in% c(20, 30), ]
n30 <- tapply(as.integer(pre$speed_limit == 30), pre$la, sum)
nr  <- tapply(pre$speed_limit, pre$la, length)
int_tab <- data.frame(la = names(nr), n_rest = as.integer(nr),
                      intensity = as.numeric(n30[names(nr)]) / as.integer(nr))
int_tab <- int_tab[order(int_tab$la), ]
cat("Intensity (pre-period 30-share of restricted-road collisions):\n")
cat("Wales range:", round(range(int_tab$intensity[grepl("^W", int_tab$la)]), 3),
    "England range:", round(range(int_tab$intensity[grepl("^E", int_tab$la)]), 3),
    "\n")
write.csv(int_tab, "output/treatment_intensity.csv", row.names = FALSE)

## ------------------------------------------------------------------
## 2. Merge into the panel (English dose = 0: nothing was reclassified)
## ------------------------------------------------------------------
panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
r <- panel[panel$grp == "restricted", ]
r$intensity <- ifelse(r$wales, int_tab$intensity[match(r$la, int_tab$la)], 0)
stopifnot(!any(is.na(r$intensity)))
r$wales_moy <- factor(ifelse(r$wales, paste0("wm", r$moy), "base"))
r$wint <- as.numeric(r$wales & r$post == 1) * r$intensity

## ------------------------------------------------------------------
## 3. Models
## ------------------------------------------------------------------
fit_nb <- function(fml, dat) {
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  list(fit = fit, se = sqrt(diag(V)), theta = th)
}

dose_fit <- function(oc, terms, label, dat) {
  cat("fitting", label, "-", oc, "... ")
  fml <- as.formula(paste(oc,
    "~ factor(la) + factor(month_idx) +", terms, "+ wales_moy + offset(log(pop))"))
  z <- fit_nb(fml, dat)
  cf <- coef(z$fit); se <- z$se
  idx <- intersect(c("tp", "wint"), names(cf))
  data.frame(outcome = oc, model = label, term = idx,
             beta = cf[idx], se = se[idx], z = cf[idx] / se[idx],
             irr = exp(cf[idx]),
             ci_lo = exp(cf[idx] - 1.96 * se[idx]),
             ci_hi = exp(cf[idx] + 1.96 * se[idx]))
}

out <- list()
for (oc in c("n_all", "adj_ksi")) {
  out[[paste(oc, "M1", sep = ".")]] <- dose_fit(oc, "tp", "M1 pooled", r)
  out[[paste(oc, "M2", sep = ".")]] <- dose_fit(oc, "wint", "M2 dose-scaled", r)
  out[[paste(oc, "M3", sep = ".")]] <- dose_fit(oc, "tp + wint", "M3 unconstrained", r)
}
res <- do.call(rbind, out)
cat("\n=== Dose-response models ===\n")
print(res[, c("outcome", "model", "term", "irr", "ci_lo", "ci_hi", "z")],
      digits = 3, row.names = FALSE)

## implied pooled effect under M2 at the Welsh mean intensity
wmean <- mean(unique(r[r$wales, c("la", "intensity")])$intensity)
cat("\nWelsh mean intensity:", round(wmean, 3), "\n")
for (oc in c("n_all", "adj_ksi")) {
  b <- res$beta[res$outcome == oc & res$model == "M2 dose-scaled"]
  s <- res$se[res$outcome == oc & res$model == "M2 dose-scaled"]
  cat(oc, ": M2 implied pooled IRR at mean intensity =",
      round(exp(b * wmean), 3), "(", round(exp(b * wmean - 1.96 * s * wmean), 3),
      "-", round(exp(b * wmean + 1.96 * s * wmean), 3), ")\n")
}

## ------------------------------------------------------------------
## 4. Two-stage check: LA-level effects vs intensity
## ------------------------------------------------------------------
map <- read.csv("output/la_effects.csv")            # adjusted KSI, model-based SEs
map$intensity <- int_tab$intensity[match(map$la, int_tab$la)]
w2 <- lm(beta ~ intensity, data = map, weights = 1 / map$se^2)
sm <- summary(w2)$coefficients
cat("\n=== Two-stage: LA-level adjusted-KSI effect vs intensity ===\n")
cat("slope:", round(sm[2, 1], 3), " (SE", round(sm[2, 2], 3),
    ", z", round(sm[2, 3], 2), ")\n")
cat("implied effect at intensity = 1:", round(exp(sm[2, 1] + sm[1, 1]), 3),
    "; at intensity = 0:", round(exp(sm[1, 1]), 3), "\n")
cor2 <- cor(map$beta, map$intensity)
cat("correlation(LA effect, intensity):", round(cor2, 2), "\n")

write.csv(res, "output/dose_response.csv", row.names = FALSE)
saveRDS(list(res = res, two_stage = sm, intensity = int_tab, map = map),
        "output/dose_response.rds")
cat("done\n")
