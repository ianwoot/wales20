## Figures for the Wales 20mph DiD study. ggplot2 for plotting only.

# There are several plots here that were explored but not used, as well as interim versions 
# of eventual plots. They're left in this script for completeness.

suppressWarnings(suppressMessages({ library(ggplot2); library(sf) }))

panel <- readRDS("data/panel.rds")
panel$adj_ksi <- panel$adj_fatal + panel$adj_serious
es <- read.csv("output/event_study_full.csv")
tab <- read.csv("output/did_main.csv")
perm <- read.csv("output/perm_draws.csv")
map_tab <- read.csv("output/la_effects.csv")

theme_set(theme_bw(base_size = 11) +
           theme(panel.grid.minor = element_blank(),
                 strip.background = element_rect(fill = "grey92", linewidth = 0),
                 legend.position = "bottom"))
INTERVENTION <- as.Date("2023-09-17")

## ==================================================================
## 1. Event-study coefficient plot (centrepiece)
## ==================================================================
es$out_lab <- c(n_ksi = "All road users: KSI",
                n_ped_ksi = "Pedestrians: KSI",
                n_cyc_ksi = "Cyclists: KSI")[es$outcome]
es <- es[es$outcome %in% c("n_ksi", "n_ped_ksi", "n_cyc_ksi"), ]
es$lo_irr <- exp(es$lo); es$hi_irr <- exp(es$hi)

## WLS trend fits through the coefficients (weights 1/se^2), pre and post.
## The pre-slope tests parallel trends (the reference-month level shift
## cancels in a slope); the post-slope tests fade-out or deepening.
trend_lines <- do.call(rbind, lapply(split(es, es$outcome), function(d) {
  mk <- function(rng, phase) {
    dd <- d[d$et >= rng[1] & d$et <= rng[2], ]
    if (phase == "pre") dd <- dd[!(dd$et >= -42 & dd$et <= -18), ]  # exclude COVID window from the pre-trend line only
    f <- lm(beta ~ et, data = dd, weights = 1 / dd$se^2)
    nd <- data.frame(et = seq(rng[1], rng[2], length.out = 30))
    pr <- predict(f, newdata = nd, se.fit = TRUE)
    data.frame(et = nd$et, y = pr$fit,
               lo = pr$fit - 1.96 * pr$se.fit,
               hi = pr$fit + 1.96 * pr$se.fit,
               phase = phase, slope = coef(f)[2],
               z = summary(f)$coefficients[2, 3], outcome = d$outcome[1])
  }
  rbind(mk(c(min(d$et), -2), "pre"), mk(c(0, max(es$et)), "post"))
}))
trend_lines$out_lab <- es$out_lab[match(trend_lines$outcome, es$outcome)]
pre_t <- trend_lines[trend_lines$phase == "pre", ]
post_t <- trend_lines[trend_lines$phase == "post", ]

slope_ann <- do.call(rbind, lapply(split(trend_lines, trend_lines$outcome), function(d) {
  p <- d[d$phase == "pre", ][1, ]; q <- d[d$phase == "post", ][1, ]
  data.frame(outcome = p$outcome,
             out_lab = p$out_lab,
             lab = paste0("pre-trend ", sprintf("%+.1f%%/mo (z %+.1f)", 100 * p$slope, p$z),
                          "\npost-trend ", sprintf("%+.1f%%/mo (z %+.1f)", 100 * q$slope, q$z)))
}))

p1 <- ggplot(es, aes(x = et, y = beta)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "red3", linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(data = pre_t, aes(x = et, ymin = lo, ymax = hi), alpha = 0.10,
              fill = "grey40", inherit.aes = FALSE) +
  geom_line(data = pre_t, aes(x = et, y = y), colour = "grey30", linetype = "22",
            linewidth = 0.55, inherit.aes = FALSE) +
  geom_ribbon(data = post_t, aes(x = et, ymin = lo, ymax = hi), alpha = 0.10,
              fill = "navy", inherit.aes = FALSE) +
  geom_line(data = post_t, aes(x = et, y = y), colour = "navy", linewidth = 0.55,
            inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "steelblue") +
  geom_point(size = 1.1, colour = "steelblue") +
  geom_text(data = slope_ann, aes(x = -24, y = Inf, label = lab),
            hjust = 0, vjust = 1.4, size = 2.6, colour = "grey20",
            inherit.aes = FALSE) +
  facet_wrap(~ out_lab, ncol = 3) +
  scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 12, 27),
                     labels = c(-80, -60, -40, -20, 0, 12, 27)) +
  labs(x = "Months since intervention (Sep 2023)",
       y = "Coefficient (log IRR vs England)",
       title = "Event-study estimates with fitted pre- and post-intervention trends",
       subtitle = "Negative binomial, LA and month fixed effects, population offset, LA-clustered 95% CIs. Reference: August 2023.\nDashed grey: WLS pre-trend fit; navy: WLS post-trend fit (weights 1/SE\u00b2); bands are approximate \u00b11.96 SE.",
       caption = "Pre-slopes are statistically null in all panels (no differential pre-trend once the reference-month level shift cancels in the slope); post-slopes are null (no fade-out or deepening within the window). Levels are normalised to Aug 2023 and are not directly comparable to the static DiD IRRs in Section 3.") +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "IRR"))

ggsave("figures/fig1_event_study.png", p1 + coord_cartesian(ylim = c(-3, NA)), width = 11.5, height = 4.6, dpi = 200)

## ==================================================================
## 1a. Event study for ALL CASUALTIES by road user group:
##     all road users, pedestrians, cyclists. No titles or captions.
## ==================================================================
es_all <- read.csv("output/event_study_full.csv")
es_all <- es_all[es_all$outcome == "n_all", ]
es_ru <- es_all
es_ru$lo <- es_ru$beta - 1.96 * es_ru$se
es_ru$hi <- es_ru$beta + 1.96 * es_ru$se
es_ru$out_lab <- c(n_all = "All road users", n_ped = "Pedestrians",
                   n_cyc = "Cyclists")[es_ru$outcome]
mk_tl <- function(d, rng, phase) {
  dd <- d[d$et >= rng[1] & d$et <= rng[2], ]
    if (phase == "pre") dd <- dd[!(dd$et >= -42 & dd$et <= -18), ]  # exclude COVID window from the pre-trend line only
  f <- lm(beta ~ et, data = dd, weights = 1 / dd$se^2)
  nd <- data.frame(et = seq(rng[1], rng[2], length.out = 30))
  pr <- predict(f, newdata = nd, se.fit = TRUE)
  data.frame(et = nd$et, y = pr$fit, lo = pr$fit - 1.96 * pr$se.fit,
             hi = pr$fit + 1.96 * pr$se.fit, phase = phase,
             out_lab = d$out_lab[1])
}
tl_ru <- do.call(rbind, lapply(split(es_ru, es_ru$outcome), function(d)
  rbind(mk_tl(d, c(min(d$et), -2), "pre"), mk_tl(d, c(0, max(d$et)), "post"))))

p1a <- ggplot(es_ru, aes(x = et, y = beta)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "red3", linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(data = tl_ru[tl_ru$phase == "pre", ],
              aes(x = et, ymin = lo, ymax = hi), alpha = 0.10,
              fill = "grey40", inherit.aes = FALSE) +
  geom_line(data = tl_ru[tl_ru$phase == "pre", ],
            aes(x = et, y = y), colour = "grey30", linetype = "22",
            linewidth = 0.55, inherit.aes = FALSE) +
  geom_ribbon(data = tl_ru[tl_ru$phase == "post", ],
              aes(x = et, ymin = lo, ymax = hi), alpha = 0.10,
              fill = "navy", inherit.aes = FALSE) +
  geom_line(data = tl_ru[tl_ru$phase == "post", ],
            aes(x = et, y = y), colour = "navy", linewidth = 0.55,
            inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "steelblue") +
  geom_point(size = 1.1, colour = "steelblue") +
  facet_wrap(~ out_lab, ncol = 3) +
  scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 12, 27),
                     labels = c(-80, -60, -40, -20, 0, 12, 27)) +
  labs(x = "Months since intervention (Sep 2023)",
       y = "Coefficient (log IRR vs England)") +
  scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = "IRR"))

ggsave("figures/fig1a_event_study_all.png", p1a + coord_cartesian(ylim = c(-3, NA)), width = 11.5, height = 4.3, dpi = 200)

## ==================================================================
## 1b. Combined 2 x 3 event-study grid: rows = severity (all severities,
##     KSI), columns = road user (all road users, cyclists, pedestrians).
## ==================================================================
es_a <- read.csv("output/event_study_full.csv")
es_a <- es_a[es_a$outcome %in% c("n_all", "n_ksi", "n_ped_ksi", "n_cyc_ksi", "n_ped", "n_cyc"), ]
es_b <- NULL
es6 <- rbind(es_a, es_b)
es6$lo <- es6$beta - 1.96 * es6$se
es6$hi <- es6$beta + 1.96 * es6$se
es6$group <- c(n_all = "All road users", n_cyc = "Cyclists",
               n_ped = "Pedestrians", n_ksi = "All road users",
               n_cyc_ksi = "Cyclists", n_ped_ksi = "Pedestrians")[es6$outcome]
es6$severity <- c(n_all = "All severities", n_cyc = "All severities",
                  n_ped = "All severities", n_ksi = "KSI",
                  n_cyc_ksi = "KSI", n_ped_ksi = "KSI")[es6$outcome]
es6$group <- factor(es6$group, levels = c("All road users", "Cyclists", "Pedestrians"))
es6$severity <- factor(es6$severity, levels = c("All severities", "KSI"))
mk_tl2 <- function(d, rng, phase) {
  dd <- d[d$et >= rng[1] & d$et <= rng[2], ]
    if (phase == "pre") dd <- dd[!(dd$et >= -42 & dd$et <= -18), ]  # exclude COVID window from the pre-trend line only
  f <- lm(beta ~ et, data = dd, weights = 1 / dd$se^2)
  nd <- data.frame(et = seq(rng[1], rng[2], length.out = 30))
  pr <- predict(f, newdata = nd, se.fit = TRUE)
  data.frame(et = nd$et, y = pr$fit, lo = pr$fit - 1.96 * pr$se.fit,
             hi = pr$fit + 1.96 * pr$se.fit, phase = phase,
             group = d$group[1], severity = d$severity[1])
}
tl6 <- do.call(rbind, lapply(split(es6, paste(es6$group, es6$severity)), function(d)
  rbind(mk_tl2(d, c(min(d$et), -2), "pre"), mk_tl2(d, c(0, max(d$et)), "post"))))

p1b <- ggplot(es6, aes(x = et, y = beta)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "red3", linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(data = tl6[tl6$phase == "pre", ],
              aes(x = et, ymin = lo, ymax = hi), alpha = 0.10,
              fill = "grey40", inherit.aes = FALSE) +
  geom_line(data = tl6[tl6$phase == "pre", ],
            aes(x = et, y = y), colour = "grey30", linetype = "22",
            linewidth = 0.45, inherit.aes = FALSE) +
  geom_ribbon(data = tl6[tl6$phase == "post", ],
              aes(x = et, ymin = lo, ymax = hi), alpha = 0.10,
              fill = "navy", inherit.aes = FALSE) +
  geom_line(data = tl6[tl6$phase == "post", ],
            aes(x = et, y = y), colour = "navy", linewidth = 0.45,
            inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "steelblue") +
  geom_point(size = 0.9, colour = "steelblue") +
  facet_grid(severity ~ group) +
  scale_x_continuous(breaks = c(-80, -40, 0, 27),
                     labels = c(-80, -40, 0, 27)) +
  labs(x = "Months since intervention (Sep 2023)",
       y = "Coefficient (log IRR vs England)") +
  theme(panel.spacing = unit(0.6, "lines"),
        strip.text = element_text(size = 8, face = "bold"))

ggsave("figures/fig1_combined.png", p1b + coord_cartesian(ylim = c(-3, NA)), width = 11.5, height = 7, dpi = 200)

## ==================================================================
## 1c. Annotated combined grid: ITS statistics per panel.
##     Pre-slope / step / post-slope labels (with 95% CIs) in the top
##     band, curved arrows to the fitted trend lines and the step.
## ==================================================================
ann <- read.csv("output/its_annotations.csv")
get_tl_y <- function(g, sev, phase, target) {
  sub <- tl6[as.character(tl6$group) == g & as.character(tl6$severity) == sev &
             tl6$phase == phase, ]
  sub[which.min(abs(sub$et - target)), "y"]
}
lab_df <- do.call(rbind, lapply(split(es6, es6$outcome), function(d) {
  o <- d$outcome[1]; a <- ann[ann$outcome == o, ]
  g <- as.character(d$group[1]); sev <- as.character(d$severity[1])
  ytxt <- max(d$hi) + 0.55
  pre_y  <- get_tl_y(g, sev, "pre", -50) + 0.10
  post_y <- get_tl_y(g, sev, "post", 15) + 0.10
  beta0  <- d$beta[d$et == 0]
  data.frame(
    group = g, severity = sev,
    label = c(sprintf("Pre %+.3f/mo\n[%+.3f, %+.3f]",
                      a$pre_slope, a$pre_lo, a$pre_hi),
              sprintf("Step %+.2f\n[%+.2f, %+.2f]",
                      a$step, a$step_lo, a$step_hi),
              sprintf("Bend %+.3f/mo\n[%+.3f, %+.3f]",
                      a$slope_chg, a$chg_lo, a$chg_hi),
              sprintf("Post %+.3f/mo\n[%+.3f, %+.3f]",
                      a$post_slope, a$post_lo, a$post_hi)),
    x = c(-79, -52, -25, 27), y = rep(ytxt, 4), hjust = c(0, 0, 0, 1),
    ax = c(-66, -40, -13, 13), ay = rep(ytxt - 0.10, 4),
    bx = c(-50, 0, 6, 22), by = c(pre_y, beta0 + 0.15,
              get_tl_y(g, sev, "post", 6) + 0.10, get_tl_y(g, sev, "post", 22) + 0.10))
}))
lab_df$group <- factor(lab_df$group, levels = levels(es6$group))
lab_df$severity <- factor(lab_df$severity, levels = levels(es6$severity))

p1c <- p1b +
  geom_curve(data = lab_df,
             aes(x = ax, y = ay, xend = bx, yend = by),
             curvature = 0.18, linewidth = 0.3, colour = "grey35",
             arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
             inherit.aes = FALSE) +
  geom_text(data = lab_df,
            aes(x = x, y = y, label = label, hjust = hjust),
            size = 2.4, lineheight = 0.9, colour = "grey15", vjust = 0.5,
            inherit.aes = FALSE)

ggsave("figures/fig1_its.png", p1c + coord_cartesian(ylim = c(-3, NA)),
       width = 12, height = 7.5, dpi = 200)

## ==================================================================
## 2. Counterfactual trend chart
## ==================================================================
r <- panel[panel$grp == "restricted", ]
agg <- function(d) aggregate(cbind(y = d$n_ksi) ~ month_idx + wales, data = d, FUN = sum)
m_series <- function(is_wales) {
  s <- agg(r[r$wales == is_wales, ])
  s$t <- s$month_idx; s$moy <- (s$month_idx - 1) %% 12 + 1; s$pre <- s$month_idx < 81
  s
}
w <- m_series(TRUE); e <- m_series(FALSE)
fit_pre <- function(s) glm(y ~ factor(moy) + t, data = s[s$pre, ], family = poisson())
wf <- fit_pre(w); ef <- fit_pre(e)
w$pred <- predict(wf, newdata = w, type = "response")
e$pred <- predict(ef, newdata = e, type = "response")
e$dev <- e$y / e$pred                          # England post deviation from own pre-trend
cf <- w$pred * e$dev                           # England-adjusted Welsh counterfactual
w$counterfactual <- cf
w$date <- as.Date(sprintf("%d-%02d-01", 2017 + (w$month_idx - 1) %/% 12,
                          w$moy))
avoided <- sum(w$y[w$month_idx >= 81] - w$counterfactual[w$month_idx >= 81])

p2 <- ggplot(w, aes(x = date)) +
  geom_line(aes(y = counterfactual), colour = "grey30", linetype = "dashed", linewidth = 0.7) +
  geom_point(aes(y = y), size = 0.7, colour = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = INTERVENTION, colour = "red3", linewidth = 0.5) +
  annotate("text", x = INTERVENTION, y = max(w$y) * 0.98,
           label = " 17 Sep 2023\ndefault 20 mph", hjust = 0, size = 3, colour = "red3") +
  annotate("text", x = as.Date("2020-06-01"), y = max(w$y) * 0.25,
           label = paste0("dashed: counterfactual\n(England-adjusted pre-trend)\nKSI avoided, Sep 2023-Dec 2025: ", round(avoided, 0)),
           size = 3, colour = "grey30", hjust = 0) +
  labs(x = NULL, y = "KSI casualties per month, Wales",
       title = "Observed vs counterfactual: killed or seriously injured on 20/30 mph roads",
       caption = "Counterfactual = Wales pre-trend (month-of-year + linear trend) scaled by England's post-period deviation from its own pre-trend.",
       subtitle = "Points: observed monthly counts; the 2020 dip reflects COVID lockdowns. Dashed: projected without treatment.") +
  coord_cartesian(ylim = c(0, max(w$y) * 1.05))

ggsave("figures/fig2_counterfactual.png", p2, width = 8.5, height = 5, dpi = 200)

## ==================================================================
## 3. DiD bar chart: IRRs by outcome, population offset, with 40+ placebo row
## ==================================================================
fit <- read.csv("output/model_fit.csv")
r2_of <- function(oc) sprintf("%.2f", fit$ct_r2[fit$outcome == oc &
                                                fit$spec != "placebo (40+ mph roads)"])
keep <- tab[tab$spec == "pop offset", ]
plac <- tab[tab$spec == "placebo 40+ mph" & tab$label == "All casualties", ]
plac$label <- "All casualties (40+ mph roads: placebo)"
keep <- rbind(keep, plac)
ord <- c("All casualties", "KSI (reported)", "KSI (severity-adjusted)",
         "Fatal", "Slight", "Pedestrian casualties", "Pedestrian KSI",
         "Cyclist casualties", "Cyclist KSI", "Child (0-15) KSI", "Older (65+) KSI",
         "All casualties (40+ mph roads: placebo)")
r2map <- c("All casualties" = r2_of("n_all"), "KSI (reported)" = r2_of("n_ksi"),
           "KSI (severity-adjusted)" = r2_of("adj_ksi"), "Fatal" = r2_of("n_fatal"),
           "Slight" = r2_of("n_slight"),
           "Pedestrian casualties" = r2_of("n_ped"), "Pedestrian KSI" = r2_of("n_ped_ksi"),
           "Cyclist casualties" = r2_of("n_cyc"), "Cyclist KSI" = r2_of("n_cyc_ksi"),
           "Child (0-15) KSI" = r2_of("n_child_ksi"), "Older (65+) KSI" = r2_of("n_older_ksi"),
           "All casualties (40+ mph roads: placebo)" =
             sprintf("%.2f", fit$ct_r2[fit$outcome == "n_all" &
                                       fit$spec == "placebo (40+ mph roads)"]))
keep$label <- paste0(keep$label, " (R\u00b2 = ", r2map[keep$label], ")")
ord <- paste0(ord, " (R\u00b2 = ", r2map[ord], ")")
keep$label <- factor(keep$label, levels = rev(ord))

p3 <- ggplot(keep[keep$label %in% ord, ], aes(x = irr, y = label)) +
  geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3,
                 colour = "steelblue", linewidth = 0.6) +
  geom_point(size = 2, colour = "steelblue") +
  scale_x_log10() +
  labs(x = "Incidence-rate ratio, Wales x post vs England (log scale)",
       y = NULL,
       title = NULL,
       subtitle = NULL,
       caption = NULL)

ggsave("figures/fig3_did_bars.png", p3, width = 8.5, height = 5.5, dpi = 200)

## ==================================================================
## 4. Placebo permutation distribution
## ==================================================================
wales_agg <- as.numeric(readRDS("output/models.rds")$wales_agg)
p4 <- ggplot(perm, aes(x = beta)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(xintercept = wales_agg, colour = "red3", linewidth = 0.8) +
  annotate("text", x = wales_agg, y = Inf, vjust = 2,
           label = "Wales", colour = "red3", hjust = -0.15, size = 3.5) +
  labs(x = "Placebo DiD coefficient (log IRR)",
       y = "Random pseudo-Wales groups (n = 1000)",
       title = "Permutation test: 22-authority pseudo-treatment groups drawn from England",
       subtitle = paste0("Wales estimate = ", round(wales_agg, 3),
                         " (IRR ", round(exp(wales_agg), 2), "); two-sided permutation p = ",
                         mean(abs(perm$beta) >= abs(wales_agg)), "."),
       caption = "Aggregate Poisson: group and month fixed effects, population offset; England-only placebo pool.")

ggsave("figures/fig4_permutation.png", p4, width = 7.5, height = 4.8, dpi = 200)

## ==================================================================
## 4a. Permutation test for ALL CASUALTIES (the stronger effect)
## ==================================================================
perm_all <- read.csv("output/perm_draws_all.csv")
pal <- readRDS("output/perm_all.rds")
p4a <- ggplot(perm_all, aes(x = beta)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(xintercept = pal$wales_agg, colour = "red3", linewidth = 0.8) +
  annotate("text", x = pal$wales_agg, y = Inf, vjust = 2,
           label = "Wales", colour = "red3", hjust = -0.15, size = 3.5) +
  labs(x = "Placebo DiD coefficient (log IRR)",
       y = "Random pseudo-Wales groups (n = 1000)",
       title = "Permutation test: all casualties on 20/30 mph roads",
       subtitle = paste0("Wales estimate = ", round(pal$wales_agg, 3),
                         " (IRR ", round(exp(pal$wales_agg), 2),
                         "); two-sided permutation p < 0.001.\nNo pseudo-treatment group comes within a quarter of the Welsh divergence."),
       caption = "Aggregate Poisson: group and month fixed effects, population offset; England-only placebo pool. Compare Figure 4 (KSI, p = 0.065).")

ggsave("figures/fig4a_permutation_all.png", p4a, width = 7.5, height = 4.8, dpi = 200)

## ==================================================================
## 5. Map of LA-level effects
## ==================================================================
shp <- st_read("data/raw/geo/lad_boundaries.geojson", quiet = TRUE)
names(shp)[1:2] <- c("la", "name")
shp <- shp[grepl("^E|^W", shp$la), ]
shp$irr <- map_tab$irr[match(shp$la, map_tab$la)]
shp$ci_lo <- map_tab$ci_lo[match(shp$la, map_tab$la)]
shp$ci_hi <- map_tab$ci_hi[match(shp$la, map_tab$la)]
shp$fill_col <- cut(shp$irr, breaks = c(0, 0.7, 0.85, 1, 1.15, 1.3, Inf),
                    labels = c("<0.70", "0.70-0.85", "0.85-1.00",
                               "1.00-1.15", "1.15-1.30", ">1.30"))

p5 <- ggplot() +
  geom_sf(data = shp[is.na(shp$irr), ], fill = "grey88", colour = "white", linewidth = 0.1) +
  geom_sf(data = shp[!is.na(shp$irr), ],
          aes(fill = fill_col), colour = "grey20", linewidth = 0.15) +
  scale_fill_brewer(palette = "RdBu", direction = 1, drop = FALSE,
                    name = "Adjusted-KSI IRR\n(Wales x post)") +
  labs(title = "Authority-level effects on 20/30 mph roads, Sep 2023 - Dec 2025",
       subtitle = "Wales: LA-specific DiD IRRs for severity-adjusted KSI (negative binomial, model-based 95% CIs in output/la_effects.csv). England in grey.",
       caption = "Population offset; LA and calendar-month fixed effects. Adjusted KSI removes the Dyfed-Powys IBRS recording change (Q2 2023) that otherwise inflates the Carmarthenshire/Pembrokeshire/Ceredigion/Powys series.") +
  theme(legend.position = "right")

ggsave("figures/fig5_map.png", p5, width = 7.2, height = 9, dpi = 200)

## ==================================================================
## 6. Behavioural panel: pedal cycle traffic, England vs Wales
## ==================================================================
cyc <- read.csv("data/raw/traffic/tra0413_raw.csv", header = FALSE, fill = TRUE,
                colClasses = "character")
cyc <- cyc[cyc[[1]] == "TRA0413", ]
hrow <- which(grepl("^1993", cyc[[6]]))
yrs <- suppressWarnings(as.integer(gsub("[^0-9].*$", "", cyc[hrow, 6:ncol(cyc)])))
mk <- function(nm) {
  row <- cyc[cyc[[2]] == nm, ]
  v <- suppressWarnings(as.numeric(unlist(row[1, 6:ncol(cyc)])))
  data.frame(country = nm, year = yrs, bn_km = v)
}
cyc_df <- rbind(mk("England"), mk("Wales"))
cyc_df <- cyc_df[!is.na(cyc_df$year) & cyc_df$year >= 2010, ]
idx <- function(d, base = 2019) d$bn_km / d$bn_km[d$year == base] * 100
cyc_df$index <- unlist(lapply(split(cyc_df, cyc_df$country), idx))

p6 <- ggplot(cyc_df, aes(x = year, y = index, colour = country)) +
  geom_hline(yintercept = 100, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 2023.7, colour = "red3", linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(England = "grey30", Wales = "red3")) +
  labs(x = NULL, y = "Pedal cycle traffic index (2019 = 100)",
       colour = NULL,
       title = "Objective behavioural trace: pedal cycle vehicle-km",
       subtitle = "DfT road traffic estimates (TRA0413), England vs Wales",
       caption = "No post-2023 divergence in cycling volume between Wales and England. ALS/NTS survey measures are England-only; see report.")

ggsave("figures/fig6_cycle_traffic.png", p6, width = 7.5, height = 4.5, dpi = 200)

## ==================================================================
## 5a. Discovery map (exploratory): colour only Welsh LAs that are
##     Benjamini-Hochberg discoveries at FDR 5%; other Welsh LAs white,
##     England grey. BH-adjusted p-values from output/la_effects_mc.csv.
## ==================================================================
mc <- read.csv("output/la_effects_mc.csv")
shp$p_bh <- mc$p_bh[match(shp$la, mc$la)]
shp$sig <- ifelse(is.na(shp$p_bh), NA,
           ifelse(shp$p_bh <= 0.05 & shp$irr < 1, "sig_below",
           ifelse(shp$p_bh <= 0.05 & shp$irr > 1, "sig_above", "null")))
eng_shp <- shp[is.na(shp$p_bh), ]
wls_null <- shp[!is.na(shp$p_bh) & shp$sig == "null", ]
wls_sig <- shp[shp$sig %in% c("sig_below", "sig_above"), ]
wls_sig$sig_lab <- ifelse(wls_sig$sig == "sig_below",
                          "Significant reduction (IRR < 1)",
                          "Significant increase (IRR > 1)")

p5a <- ggplot() +
  geom_sf(data = eng_shp, fill = "grey88", colour = "white", linewidth = 0.1) +
  geom_sf(data = wls_null, fill = "white", colour = "grey55", linewidth = 0.15) +
  geom_sf(data = wls_sig, aes(fill = sig_lab), colour = "grey20", linewidth = 0.15) +
  scale_fill_manual(values = c("Significant reduction (IRR < 1)" = "#2166ac",
                               "Significant increase (IRR > 1)" = "#b2182b"),
                    name = NULL) +
  theme(legend.position = "right")

ggsave("figures/fig5a_map_significant.png", p5a, width = 7.2, height = 9, dpi = 200)

cat("figures written\n")
