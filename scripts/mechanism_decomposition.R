## Mechanism decomposition: which link in the causal chain did the policy act on?
##
## The all-casualty effect (per vehicle-km) decomposes into:
##   (a) collisions per vehicle-km        -- did crashes become less frequent?
##   (b) casualties per collision          -- given a crash, fewer people hurt?
## and the KSI effect adds:
##   (c) KSI per casualty                 -- given an injury, less severe?
##
## Each link is estimated with the same DiD specification as the headline
## models (NB two-step, LA FE, calendar-month FE, Wales seasonality,
## LA-clustered SEs); the offsets chain the links together:
##   (a) offset = log(vkm); (b) offset = log(collisions); (c) offset = log(casualties)
## The product of the link IRRs should approximately recover the headline
## all-casualty (a x b) and adjusted-KSI (a x b x c) vehicle-km estimates.
##
## Zero cells: rates are only defined where the denominator is positive, so
## (b) is fitted on cells with >=1 collision and (c) on cells with >=1
## casualty; y is necessarily 0 in the dropped cells, so they carry no rate
## information. This makes the product an approximate, not exact, identity.
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

## ------------------------------------------------------------------
## 1. Panel skeleton (population, vkm, casualties) from the built panel
## ------------------------------------------------------------------
panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
skel <- panel[, c("la", "month_idx", "grp", "year", "moy", "pop", "vkm_m_km",
                  "wales", "post", "tp", "n_all", "n_ksi", "adj_ksi")]
skel <- skel[!duplicated(skel), ]

## ------------------------------------------------------------------
## 2. Collision counts per LA x month x road group (from extracted file)
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
sl <- col$speed_limit
col$grp <- ifelse(sl %in% c(20, 30), "restricted", ifelse(sl >= 40, "higher", NA))
col <- col[!is.na(col$grp), ]
coll_agg <- as.data.frame(table(la = col$la, month_idx = col$month_idx,
                                grp = col$grp),
                          responseName = "n_coll")
coll_agg$month_idx <- as.integer(as.character(coll_agg$month_idx))
coll_agg$la <- as.character(coll_agg$la); coll_agg$grp <- as.character(coll_agg$grp)

## ------------------------------------------------------------------
## 3. Merge and build the mechanism panel (restricted roads)
## ------------------------------------------------------------------
m <- merge(skel[skel$grp == "restricted", ], coll_agg, all.x = TRUE)
m$n_coll[is.na(m$n_coll)] <- 0
m$wales_moy <- factor(ifelse(m$wales, paste0("wm", m$moy), "base"))
cat("mechanism panel:", nrow(m), "rows; total collisions:",
    sum(m$n_coll), "; total casualties:", sum(m$n_all), "\n")

## ------------------------------------------------------------------
## 4. Estimation helpers (same machinery as the headline models)
## ------------------------------------------------------------------
fit_nb <- function(fml, dat) {
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  list(fit = fit, se = sqrt(diag(V)), theta = th)
}
est <- function(label, outcome, offset, dat) {
  cat("fitting", label, "... ")
  fml <- as.formula(paste(outcome,
    "~ factor(la) + factor(month_idx) + tp + wales_moy + offset(", offset, ")"))
  z <- fit_nb(fml, dat)
  cf <- coef(z$fit); b <- cf["tp"]; s <- z$se["tp"]
  cat("IRR", round(exp(b), 3), "\n")
  data.frame(link = label, outcome = outcome, offset = offset,
             beta = b, se = s, irr = exp(b),
             ci_lo = exp(b - 1.96 * s), ci_hi = exp(b + 1.96 * s),
             theta = z$theta, n = nrow(dat))
}

## ------------------------------------------------------------------
## 5. The three links
## ------------------------------------------------------------------
res <- list()

## (a) collisions per vehicle-km
res[["a"]] <- est("a: collisions per vehicle-km", "n_coll", "log(vkm_m_km * 1e6)", m)

## (b) casualties per collision (cells with >=1 collision)
mb <- m[m$n_coll >= 1, ]
res[["b"]] <- est("b: casualties per collision", "n_all", "log(n_coll)", mb)

## (c) KSI per casualty, reported and adjusted (cells with >=1 casualty)
mc1 <- m[m$n_all >= 1, ]
res[["c1"]] <- est("c: KSI per casualty (reported)", "n_ksi", "log(n_all)", mc1)
res[["c2"]] <- est("c: KSI per casualty (adjusted)", "adj_ksi", "log(n_all)", mc1)

out <- do.call(rbind, res)
write.csv(out, "output/mechanism_decomposition.csv", row.names = FALSE)

## ------------------------------------------------------------------
## 6. Product check against the headline vehicle-km estimates
## ------------------------------------------------------------------
tab <- read.csv("output/did_main.csv")
head_all <- tab[tab$spec == "vkm offset" & tab$label == "All casualties", ]
head_ksi <- tab[tab$spec == "vkm offset" & tab$label == "KSI (severity-adjusted)", ]
cat("\n=== Product check ===\n")
irr_a <- out$irr[grepl("^a:", out$link)]
irr_b <- out$irr[grepl("^b:", out$link)]
irr_c2 <- out$irr[grepl("^c: KSI per casualty \\(adjusted)", out$link)]
prod_all <- irr_a * irr_b
prod_ksi <- prod_all * irr_c2
cat("all casualties: a x b =", round(prod_all, 3),
    " vs headline vkm estimate", round(head_all$irr, 3), "\n")
cat("adjusted KSI:  a x b x c_adj =", round(prod_ksi, 3),
    " vs headline vkm estimate", round(head_ksi$irr, 3), "\n")
cat("note: product is approximate - links (b) and (c) drop zero-denominator",
    "cells and each link has its own dispersion parameter.\n")
cat("\n=== Summary ===\n")
print(out[, c("link", "irr", "ci_lo", "ci_hi")], digits = 3, row.names = FALSE)
cat("done\n")
