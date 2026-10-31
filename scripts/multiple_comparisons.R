## Multiple-comparisons corrections for the 22 LA-level map estimates.
## Holm-Bonferroni (familywise error) and Benjamini-Hochberg (FDR) adjusted p-values.
m <- read.csv("output/la_effects.csv")
m$p <- 2 * pnorm(-abs(m$beta / m$se))
m$p_holm <- p.adjust(m$p, method = "holm")
m$p_bh   <- p.adjust(m$p, method = "BH")
m <- m[order(m$p), ]
suppressWarnings(suppressMessages(library(sf)))
shp <- st_read("data/raw/geo/lad_boundaries.geojson", quiet = TRUE)
nm <- setNames(shp$LAD24NM, shp$LAD24CD)
m$name <- nm[m$la]
cat("Ordered by raw p:\n")
print(m[, c("name", "irr", "p", "p_holm", "p_bh")], row.names = FALSE, digits = 4)
cat("\nHolm significant at 0.05:", sum(m$p_holm < 0.05), "\n")
cat("BH significant at 0.05 (FDR):", sum(m$p_bh <= 0.05), "\n")
write.csv(m, "output/la_effects_mc.csv", row.names = FALSE)
