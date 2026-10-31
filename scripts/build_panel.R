## Build LA x month casualty panel for the Wales 20mph difference-in-differences study.
## Base R throughout. Outputs data/panel.rds.

suppressWarnings(suppressMessages(library(readxl)))

raw <- "data/raw"
dir.create("data", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## ------------------------------------------------------------------
## 1. Read extracted STATS19 files (2017-2024)
## ------------------------------------------------------------------
col <- read.csv(file.path(raw, "stats19/collisions_2017_2025.csv"),
                colClasses = c("character", "integer", "character", "character",
                               "integer", "integer", "integer", "integer"))
cas <- read.csv(file.path(raw, "stats19/casualties_2017_2025.csv"),
                colClasses = c("character", "integer", "integer", "integer",
                               "integer", "integer", "integer", "numeric", "numeric"))
cat("collisions:", nrow(col), " casualties:", nrow(cas), "\n")

## ------------------------------------------------------------------
## 2. LA code crosswalk: defunct codes -> current LAD codes
##    (all 2019-2023 English reorganisations, verified empirically)
## ------------------------------------------------------------------
xwalk <- c(
  ## BCP 2019
  "E06000028" = "E06000058", "E06000029" = "E06000058", "E07000048" = "E06000058",
  ## Dorset 2019
  "E07000049" = "E06000059", "E07000050" = "E06000059", "E07000051" = "E06000059",
  "E07000052" = "E06000059", "E07000053" = "E06000059",
  ## Somerset West and Taunton 2019 (then Somerset 2023)
  "E07000190" = "E06000066", "E07000191" = "E06000066", "E07000246" = "E06000066",
  ## Buckinghamshire 2020
  "E07000004" = "E06000060", "E07000005" = "E06000060",
  "E07000006" = "E06000060", "E07000007" = "E06000060",
  ## Suffolk 2019/2020
  "E07000201" = "E07000245", "E07000204" = "E07000245",
  "E07000205" = "E07000244", "E07000206" = "E07000244",
  ## North Northamptonshire 2021
  "E07000150" = "E06000061", "E07000152" = "E06000061",
  "E07000153" = "E06000061", "E07000156" = "E06000061",
  ## West Northamptonshire 2021
  "E07000151" = "E06000062", "E07000154" = "E06000062", "E07000155" = "E06000062",
  ## Cumberland 2023
  "E07000026" = "E06000063", "E07000028" = "E06000063", "E07000029" = "E06000063",
  ## Westmorland and Furness 2023
  "E07000027" = "E06000064", "E07000030" = "E06000064", "E07000031" = "E06000064",
  ## Somerset 2023
  "E07000187" = "E06000066", "E07000188" = "E06000066", "E07000189" = "E06000066",
  ## North Yorkshire 2023
  "E07000163" = "E06000065", "E07000164" = "E06000065", "E07000165" = "E06000065",
  "E07000166" = "E06000065", "E07000167" = "E06000065", "E07000168" = "E06000065",
  "E07000169" = "E06000065",
  ## Heathrow airport pseudo-code -> Hillingdon
  "EHEATHROW" = "E09000017")

## ------------------------------------------------------------------
## 3. Merge casualties with collision attributes, clean
## ------------------------------------------------------------------
col$la <- unname(xwalk[col$la_ons]); col$la[is.na(col$la)] <- col$la_ons[is.na(col$la)]
if (any(grepl("^(S|N)", col$la))) warning("unexpected non-E/W codes retained")
if (any(!grepl("^(E|W|S)", col$la))) warning("unexpected codes:", col$la[!grepl("^(E|W|S)", col$la)])

col <- col[col$la != "-1", ]
keep_country <- grepl("^E|^W", col$la)          # drop Scotland
col <- col[keep_country, ]

## month index 1..96 over Jan 2017 - Dec 2024
col$y <- as.integer(substr(col$date, 7, 10))
col$m <- as.integer(substr(col$date, 4, 5))
stopifnot(!any(is.na(col$y) | is.na(col$m)), all(col$y >= 2017 & col$y <= 2025))
col$month_idx <- (col$y - 2017) * 12 + col$m

cols_keep <- c("collision_index", "la", "y", "m", "month_idx", "speed_limit",
              "road_type", "first_road_class", "urban_or_rural")
cas <- merge(cas, col[, cols_keep], by = "collision_index", all.x = FALSE)
cat("merged casualties:", nrow(cas), "\n")

## speed limit sanity
print(table(cas$speed_limit, useNA = "ifany"))

## road group: restricted (20/30) vs higher (>=40); NA for missing/other
sl <- cas$speed_limit
cas$grp <- ifelse(sl %in% c(20, 30), "restricted",
           ifelse(sl >= 40, "higher", NA_character_))
cat("casualties with unusable speed limit:", sum(is.na(cas$grp)), "\n")

## ------------------------------------------------------------------
## 4. Outcome weights per casualty
## ------------------------------------------------------------------
sev <- cas$casualty_severity           # 1 fatal, 2 serious, 3 slight
adj_s <- ifelse(cas$adj_serious < 0, 0, cas$adj_serious)   # -1 -> 0 (unknown)
X <- cbind(
  n_all      = 1,
  n_fatal    = as.numeric(sev == 1),
  n_slight   = as.numeric(sev == 3),
  n_ksi      = as.numeric(sev <= 2),
  adj_fatal  = as.numeric(sev == 1),
  adj_serious = adj_s,
  n_ped      = as.numeric(cas$casualty_class == 3),
  n_ped_ksi  = as.numeric(cas$casualty_class == 3 & sev <= 2),
  n_cyc      = as.numeric(cas$casualty_type == 1),
  n_cyc_ksi  = as.numeric(cas$casualty_type == 1 & sev <= 2),
  n_child    = as.numeric(cas$age_of_casualty >= 0 & cas$age_of_casualty <= 15),
  n_child_ksi = as.numeric(cas$age_of_casualty >= 0 & cas$age_of_casualty <= 15 & sev <= 2),
  n_older    = as.numeric(cas$age_of_casualty >= 65),
  n_older_ksi = as.numeric(cas$age_of_casualty >= 65 & sev <= 2))

grp_fac <- factor(paste(cas$la, cas$month_idx, cas$grp, sep = "_"))
agg <- rowsum(X, group = grp_fac, reorder = TRUE)
keys <- do.call(rbind, strsplit(rownames(agg), "_", fixed = TRUE))
panel_counts <- data.frame(la = keys[, 1], month_idx = as.integer(keys[, 2]),
                           grp = keys[, 3], agg, row.names = NULL)

## ------------------------------------------------------------------
## 5. Balanced panel frame: LA x month x group
## ------------------------------------------------------------------
myeb <- read_excel(file.path(raw, "pop/myebtablesuk20112025.xlsx"), sheet = "MYEB2", skip = 1)
pop <- myeb[grepl("^E0[6-9]|^W06", myeb$ladcode23), ]
pop_la <- aggregate(pop[, paste0("population_", 2017:2025)],
                    by = list(la = pop$ladcode23), FUN = sum)
pop_la$la <- as.character(pop_la$la)
names(pop_la) <- c("la", paste0("pop_", 2017:2025))
pop_long <- data.frame(la = rep(pop_la$la, each = 9),
                       year = rep(2017:2025, times = nrow(pop_la)),
                       pop = as.vector(t(as.matrix(pop_la[, -1]))))
pop_long$la <- as.character(pop_long$la)

las <- sort(pop_la$la)
cat("panel LAs:", length(las), " Wales:", sum(grepl("^W", las)), "\n")

months <- 1:108
frame <- expand.grid(la = las, month_idx = months, grp = c("restricted", "higher"))
frame$year <- 2017 + (frame$month_idx - 1) %/% 12
frame$moy  <- (frame$month_idx - 1) %% 12 + 1
panel <- merge(frame, panel_counts, by = c("la", "month_idx", "grp"), all.x = TRUE)
outc <- c("n_all","n_fatal","n_slight","n_ksi","adj_fatal","adj_serious",
          "n_ped","n_ped_ksi","n_cyc","n_cyc_ksi","n_child","n_child_ksi",
          "n_older","n_older_ksi")
panel[outc] <- lapply(panel[outc], function(z) { z[is.na(z)] <- 0; z })
panel <- panel[order(panel$la, panel$grp, panel$month_idx), ]

## population offset
panel <- merge(panel, pop_long, by = c("la", "year"), all.x = TRUE)
panel$wales <- grepl("^W", panel$la)

## ------------------------------------------------------------------
## 6. Traffic (TRA8904): annual vehicle-km by traffic authority
## ------------------------------------------------------------------
tra <- read.csv(file.path(raw, "traffic/tra8904_full.csv"), header = FALSE,
                colClasses = "character", skip = 1)   # row 1 is the sheet title
names(tra)[1:6] <- c("code", "region", "ita", "name", "notes", "units")
year_cols <- paste0("V", 7:39)          # up to 33 year columns
hrow <- which(tra[[7]] == "1993" | grepl("^1993", tra[[7]]))[1]
if (is.na(hrow)) stop("could not locate TRA8904 year header row")
yrs <- suppressWarnings(as.integer(gsub("[^0-9].*$", "", tra[hrow, year_cols])))
sel <- year_cols[!is.na(yrs) & yrs %in% 2017:2025]
tra_yrs <- yrs[!is.na(yrs) & yrs %in% 2017:2025]
vkm <- tra[grepl("^E0[6-9]|^E100|^W06", tra$code), c("code", sel)]
names(vkm) <- c("auth", paste0("v", tra_yrs))
vkm_mat <- suppressWarnings(apply(as.matrix(vkm[, -1]), 2, as.numeric))
vkm_long <- data.frame(auth = rep(vkm$auth, each = 9),
                       year = rep(tra_yrs, times = nrow(vkm)),
                       vkm_m_km = as.vector(t(vkm_mat)))

## district -> traffic authority: shire districts inherit county traffic;
## UAs, London boroughs and metropolitan boroughs have their own TRA8904 rows.
lu <- read.csv(file.path(raw, "lad_to_county_apr2023.csv"))
names(lu)[1:4] <- c("la", "ladnm", "cty", "ctynm")
lu <- lu[grepl("^E07", lu$la), c("la", "cty")]
panel <- merge(panel, lu, by = "la", all.x = TRUE)
la_chr <- as.character(panel$la)
panel$auth <- ifelse(is.na(panel$cty), la_chr, as.character(panel$cty))

## Authorities created by mergers (BCP/Dorset, Northamptonshire pair,
## Cumbria pair) lack pre-merger TRA8904 rows: impute from the parent
## county aggregate, split by population share of the successor areas.
parent_split <- list(
  list(children = c("E06000058", "E06000059"), parent = "E10000009", upto = 2018),
  list(children = c("E06000061", "E06000062"), parent = "E10000021", upto = 2020),
  list(children = c("E06000063", "E06000064"), parent = "E10000006", upto = 2022))
for (ps in parent_split) {
  pw <- aggregate(pop ~ year, data = pop_long[pop_long$la %in% ps$children, ], FUN = sum)
  for (ch in ps$children) {
    for (yy in 2017:ps$upto) {
      idx <- which(vkm_long$auth == ch & vkm_long$year == yy)
      if (length(idx) == 1 && !is.na(vkm_long$vkm_m_km[idx])) next
      pv <- vkm_long$vkm_m_km[vkm_long$auth == ps$parent & vkm_long$year == yy]
      share <- pop_long$pop[pop_long$la == ch & pop_long$year == yy] /
                pw$pop[pw$year == yy]
      val <- pv * share
      if (length(idx) == 1) {
        vkm_long$vkm_m_km[idx] <- val
      } else {
        vkm_long <- rbind(vkm_long, data.frame(auth = ch, year = yy, vkm_m_km = val))
      }
    }
  }
}
panel <- merge(panel, vkm_long, by.x = c("auth", "year"), by.y = c("auth", "year"), all.x = TRUE)
cat("LAs missing vkm:", paste(unique(panel$la[is.na(panel$vkm_m_km)]), collapse = " "), "\n")

## ------------------------------------------------------------------
## 7. Treatment calendar
## ------------------------------------------------------------------
## Default 20mph limit became law in Wales on 17 September 2023.
sep23 <- (2023 - 2017) * 12 + 9   # month_idx of Sept 2023
panel$post <- as.integer(panel$month_idx >= sep23)
panel$event_time <- panel$month_idx - sep23
panel$tp <- panel$wales * panel$post # important: dummy indicator for post-intervention Wales

saveRDS(panel, "data/panel.rds")

## ------------------------------------------------------------------
## 8. QA printout
## ------------------------------------------------------------------
w <- panel[panel$grp == "restricted", ]
cat("\n--- QA: annual casualties on 20/30 roads, Wales ---\n")
print(aggregate(cbind(n_all, n_ksi, adj_serious) ~ year,
                data = w[w$wales, ], FUN = sum))
cat("\n--- QA: annual casualties on 20/30 roads, England ---\n")
print(aggregate(cbind(n_all, n_ksi, adj_serious) ~ year,
                data = w[!w$wales, ], FUN = sum))
cat("\n--- QA: speed-limit mix of casualties, Wales, pre vs post ---\n")
casw <- cas[grepl("^W", cas$la), ]
print(table(casw$speed_limit, post = casw$month_idx >= sep23))
cat("\npanel rows:", nrow(panel), " LAs:", length(las), "\n")

