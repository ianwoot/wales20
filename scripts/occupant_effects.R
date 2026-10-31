## Vehicle-occupant effects: did the 20 mph default protect people INSIDE cars?
##
## The all-casualty estimate bundles occupants with VRUs; this script separates
## them. Occupants = driver or passenger (casualty_class 1 or 2) in a car, taxi
## or van (casualty_type 8, 9, 19). Same DiD machinery as the headline models
## (NB two-step, LA + month FE, Wales seasonality, population offset,
## LA-clustered SEs), on restricted roads.
suppressWarnings(suppressMessages({
  library(MASS); library(sandwich); library(lmtest)
}))

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

## collision attributes (la, month, speed-limit group)
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
ck <- col[, c("collision_index", "la", "month_idx", "grp")]

## casualties with occupant flags
cas <- read.csv(file.path(raw, "stats19/casualties_2017_2025.csv"),
                colClasses = c("character", "integer", "integer", "integer",
                               "integer", "integer", "integer", "numeric", "numeric"))
cas <- merge(cas, ck, by = "collision_index")
is_occ <- cas$casualty_class %in% c(1, 2) & cas$casualty_type %in% c(8, 9, 19)
cat("occupant casualties:", sum(is_occ), "of", nrow(cas), "in England and Wales\n")

occ <- cas[is_occ, ]
X <- cbind(
  occ_all  = 1,
  occ_ksi  = as.numeric(occ$casualty_severity <= 2),
  occ_adjksi = as.numeric(occ$casualty_severity == 1) +
              ifelse(occ$adj_serious < 0, 0, occ$adj_serious))
agg <- rowsum(X, group = factor(paste(occ$la, occ$month_idx, occ$grp, sep = "_")))
keys <- do.call(rbind, strsplit(rownames(agg), "_", fixed = TRUE))
occ_counts <- data.frame(la = keys[, 1], month_idx = as.integer(keys[, 2]),
                          grp = keys[, 3], agg, row.names = NULL)

## merge onto the panel skeleton and fit
panel <- readRDS("data/panel.rds")
skel <- panel[panel$grp == "restricted",
              c("la", "month_idx", "grp", "moy", "pop", "wales", "post", "tp")]
skel <- skel[!duplicated(skel), ]
m <- merge(skel, occ_counts, by = c("la", "month_idx", "grp"), all.x = TRUE)
for (v in c("occ_all", "occ_ksi", "occ_adjksi")) m[[v]][is.na(m[[v]])] <- 0
m$wales_moy <- factor(ifelse(m$wales, paste0("wm", m$moy), "base"))

fit_nb <- function(fml, dat) {
  p <- glm(fml, data = dat, family = poisson(), control = glm.control(maxit = 300))
  th <- MASS::theta.ml(dat[[all.vars(fml)[1]]], fitted(p), limit = 100)
  fit <- glm(fml, data = dat, family = MASS::negative.binomial(theta = th),
            control = glm.control(maxit = 300))
  V <- vcovCL(fit, cluster = dat$la, type = "HC1")
  list(fit = fit, se = sqrt(diag(V)))
}
res <- do.call(rbind, lapply(c("occ_all", "occ_ksi", "occ_adjksi"), function(oc) {
  cat("fitting", oc, "... ")
  fml <- as.formula(paste(oc,
    "~ factor(la) + factor(month_idx) + tp + wales_moy + offset(log(pop))"))
  z <- fit_nb(fml, m)
  b <- coef(z$fit)["tp"]; s <- z$se["tp"]
  cat("IRR", round(exp(b), 3), "\n")
  data.frame(outcome = oc, beta = b, se = s, irr = exp(b),
             ci_lo = exp(b - 1.96 * s), ci_hi = exp(b + 1.96 * s))
}))
write.csv(res, "output/occupant_effects.csv", row.names = FALSE)
cat("\n=== Vehicle-occupant DiD (restricted roads, pop offset) ===\n")
print(res[, c("outcome", "irr", "ci_lo", "ci_hi")], digits = 3, row.names = FALSE)
cat("done\n")
