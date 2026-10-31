## Permutation test for ALL CASUALTIES (n_all) on restricted roads,
## mirroring the KSI permutation in analysis.R section 4:
## 1,000 random pseudo-Wales groups of 22 English authorities; aggregate
## two-group Poisson with month fixed effects and population offset.
set.seed(20250917)

panel <- readRDS("data/panel.rds")
r <- panel[panel$grp == "restricted", ]
la_vec <- sort(unique(r$la))
Y <- tapply(r$n_all, list(r$la, r$month_idx), sum); Y[is.na(Y)] <- 0
pm <- r[, c("la", "month_idx", "pop")]; pm <- pm[!duplicated(pm), ]
P <- tapply(pm$pop, list(pm$la, pm$month_idx), function(z) z[1])
eng_idx <- which(grepl("^E", la_vec))
wls_idx <- which(grepl("^W", la_vec))
post_vec <- as.integer(1:108 >= 81)
month_tab <- 1:108

agg_fit <- function(treated_idx, control_idx) {
  yt <- colSums(Y[treated_idx, , drop = FALSE])
  yc <- colSums(Y[control_idx, , drop = FALSE])
  pt <- colSums(P[treated_idx, , drop = FALSE])
  pc <- colSums(P[control_idx, , drop = FALSE])
  g <- rep(c(1, 0), each = 108); post2 <- rep(post_vec, 2)
  mm <- cbind(1, g, post2, g * post2,
              sapply(month_tab[-1],
                     function(m) rep(as.integer(month_tab == m), 2)))
  storage.mode(mm) <- "double"
  fit <- glm.fit(x = mm, y = c(yt, yc), family = poisson(),
                 offset = log(c(pt, pc)))
  fit$coefficients[4]
}

t0 <- Sys.time()
perm <- numeric(1000)
for (i in 1:1000) {
  tr <- sample(eng_idx, 22)
  perm[i] <- agg_fit(tr, setdiff(eng_idx, tr))
}
wales_agg <- agg_fit(wls_idx, eng_idx)
pval <- mean(abs(perm) >= abs(wales_agg))
cat("permutation (all casualties):", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
cat("Wales agg beta:", round(wales_agg, 4), " IRR:", round(exp(wales_agg), 3),
    " two-sided perm p:", pval, "\n")
cat("placebo quantiles:", round(quantile(perm, c(.025, .5, .975)), 3), "\n")

write.csv(data.frame(draw = 1:1000, beta = perm),
          "output/perm_draws_all.csv", row.names = FALSE)
saveRDS(list(perm = perm, wales_agg = wales_agg, p = pval),
        "output/perm_all.rds")
cat("done\n")
