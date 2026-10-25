library(shiny)
library(shinyWidgets)
library(DT)
library(ggplot2)
suppressMessages(library(mvtnorm))

crukdark  <- "#00007E"
cruklight <- "#009CEE"
crukpink  <- "#FF0087"

## ---------------------------------------------------------------------------
## Adjusted p-values that p.adjust does not provide
## ---------------------------------------------------------------------------
## -expm1(k*log1p(-p)) instead of 1-(1-p)^k: the direct form rounds 1-p to 1
## and returns an adjusted p of exactly 0 for p below about 1e-16
p_sidak      <- function(p) pmin(-expm1(length(p) * log1p(-p)), 1)
p_holm_sidak <- function(p) {
  n <- length(p); o <- order(p); ro <- order(o)
  pmin(cummax(-expm1((n - seq_len(n) + 1) * log1p(-p[o])))[ro], 1)
}

## ---------------------------------------------------------------------------
## One simulated experiment: k groups of n_arm observations.
## Group 1 is the control. `n_eff` of the remaining groups carry a real effect
## of size d (in units of the control SD). Variances may differ between groups,
## and the data may be normal or right-skewed.
## ---------------------------------------------------------------------------
make_data <- function(k, n_arm, d, n_eff, sd_ratio, skew) {
  mu  <- rep(0, k)
  ## the groups WITHOUT an effect keep the control's SD of 1; the groups WITH an
  ## effect share the SD set on the slider, which may be above or below 1
  sds <- rep(1, k)
  if (n_eff > 0) {
    idx <- seq(k - n_eff + 1, k)
    mu[idx]  <- d
    sds[idx] <- sd_ratio
  }
  g   <- factor(rep(seq_len(k), each = n_arm), labels = paste0("G", seq_len(k)))
  e   <- if (skew) {
    ## standardised, right-skewed: exponential shifted and scaled to mean 0, sd 1
    rexp(k * n_arm) - 1
  } else rnorm(k * n_arm)
  y <- rep(mu, each = n_arm) + rep(sds, each = n_arm) * e
  list(y = y, g = g, mu = mu, truth = mu != 0)
}

## which pairs each family compares, and which of them are true nulls
family_pairs <- function(k, family, mu) {
  pairs <- if (family == "dunnett") cbind(1, 2:k) else t(combn(k, 2))
  ## the comparison is a true null when the two group means coincide
  list(pairs = pairs,
       is_null = apply(pairs, 1, function(ij) mu[ij[1]] == mu[ij[2]]))
}

## correlation between two contrasts in a balanced design: +1/2 when they share
## a group on the same side, -1/2 on opposite sides, 0 when they are disjoint
contrast_corr <- function(pairs) {
  m <- nrow(pairs); R <- diag(m)
  for (a in seq_len(m)) for (b in seq_len(m)) if (a != b) {
    sh <- intersect(pairs[a, ], pairs[b, ])
    if (length(sh) == 1)
      R[a, b] <- if ((pairs[a, 1] == sh) == (pairs[b, 1] == sh)) 0.5 else -0.5
  }
  R
}

## Critical values are fixed by the design, so they are computed once per
## configuration rather than once per simulated experiment. For a balanced
## design these are exact:
##   Dunnett  : equicoordinate quantile of the multivariate t with corr 1/2
##   Tukey    : studentized range, q(alpha, k, nu) / sqrt(2)
##   Steel    : same idea with the asymptotic normal rank distribution
## The step-down versions reuse the quantile for the number of comparisons
## still in play. For Dunnett that quantile is exact, because any subset of the
## control contrasts is again equicorrelated at 1/2. For Tukey it is not: the
## max-t quantile of a subset depends on WHICH pairs are left, not only on how
## many, and the range over the smallest number of groups that could produce
## them can fall below it (e.g. 6 groups whose true nulls are two disjoint
## triples). So the Tukey step-down uses the smaller of two bounds that are
## valid for every subset of size j: the single-step range value, and the Sidak
## value for j two-sided statistics (Sidak's inequality holds for t statistics
## sharing one variance estimate). Both are monotone in j, so the step-down is
## valid, and it is at least as powerful as single-step Tukey and as Holm-Sidak.
critical_values <- function(k, n_arm, family, alpha, sd_ratio = 1, n_eff = 0) {
  nu    <- k * (n_arm - 1)
  pairs <- family_pairs(k, family, rep(0, k))$pairs
  m     <- nrow(pairs)
  R     <- contrast_corr(pairs)

  if (family == "dunnett") {
    single <- qmvt(1 - alpha, tail = "both.tails", df = nu,
                   corr = R)$quantile
    ## exchangeable within any subset, so the quantile depends only on the size
    stepdown <- vapply(seq_len(m), function(j) {
      if (j == 1) qt(1 - alpha / 2, nu)
      else if (j == m) single   # same quantile; qmvt is random, so reuse it
      else qmvt(1 - alpha, tail = "both.tails", df = nu,
                corr = matrix(0.5, j, j) + diag(0.5, j))$quantile
    }, numeric(1))
  } else {
    single <- qtukey(1 - alpha, nmeans = k, df = nu) / sqrt(2)
    ## step-down: the smaller of the single-step range value and the Sidak
    ## value for the j comparisons still in play (see the note above)
    stepdown <- vapply(seq_len(m), function(j)
      min(single, qt(1 - (1 - (1 - alpha)^(1 / j)) / 2, nu)), numeric(1))
  }
  rank_single <- qmvnorm(1 - alpha, tail = "both.tails", corr = R)$quantile
  ## Hasler-Hothorn: max-t on Welch statistics. The correlation matrix and the
  ## df depend on the variances, so they are set from the DESIGN (the true SDs
  ## and n) once per configuration rather than re-estimated in every run.
  sds <- rep(1, k)
  if (n_eff > 0) sds[seq(k - n_eff + 1, k)] <- sd_ratio
  vg  <- sds^2 / n_arm
  cvp <- vg[pairs[, 1]] + vg[pairs[, 2]]
  Rh  <- diag(m)
  for (a in seq_len(m)) for (b in seq_len(m)) if (a != b) {
    sh <- intersect(pairs[a, ], pairs[b, ])
    if (length(sh) == 1) {
      sg <- if ((pairs[a, 1] == sh) == (pairs[b, 1] == sh)) 1 else -1
      Rh[a, b] <- sg * vg[sh] / sqrt(cvp[a] * cvp[b])
    }
  }
  nu_hh <- min((n_arm - 1) * cvp^2 /
                 (vg[pairs[, 1]]^2 + vg[pairs[, 2]]^2))
  hh <- qmvt(1 - alpha, tail = "both.tails", df = max(2, round(nu_hh)),
             corr = Rh, algorithm = GenzBretz(maxpts = 25000,
                                              abseps = 1e-4))$quantile

  list(nu = nu, pairs = pairs, single = single, stepdown = stepdown,
       rank_single = rank_single, m = m, smm = smm_spline(m, alpha), hh = hh)
}

## step-down max-t: test the largest statistic against the quantile for all m
## comparisons, then the next against the quantile for m-1, and stop at the
## first failure
stepdown_reject <- function(absT, crit) {
  m <- length(absT); o <- order(absT, decreasing = TRUE)
  rej <- logical(m)
  for (s in seq_len(m)) {
    if (absT[o[s]] > crit[m - s + 1]) rej[o[s]] <- TRUE else break
  }
  rej
}

## Welch statistic and Satterthwaite df for one pair
welch_pair <- function(y, g, i, j) {
  a <- y[g == levels(g)[i]]; b <- y[g == levels(g)[j]]
  na <- length(a); nb <- length(b)
  va <- var(a) / na; vb <- var(b) / nb
  c(t = (mean(b) - mean(a)) / sqrt(va + vb),
    df = (va + vb)^2 / (va^2 / (na - 1) + vb^2 / (nb - 1)))
}

## standardised Wilcoxon rank sum computed within one pair (Steel / Steel-Dwass)
rank_pair <- function(y, g, i, j) {
  a <- y[g == levels(g)[i]]; b <- y[g == levels(g)[j]]
  na <- length(a); nb <- length(b)
  r <- rank(c(a, b)); W <- sum(r[seq_len(na)])
  (W - na * (na + nb + 1) / 2) / sqrt(na * nb * (na + nb + 1) / 12)
}

method_sets <- list(
  dunnett = c("Bonferroni", "Holm", "Šidák", "Holm–Šidák",
              "Dunnett", "Holm–Dunnett",
              "Hasler–Hothorn", "Hasler–Hothorn Šidák",
              "Steel", "Steel Šidák",
              "Dunn + Bonferroni", "Dunn + Holm"),
  tukey   = c("Bonferroni", "Holm", "Šidák", "Holm–Šidák",
              "Tukey", "Holm–Tukey",
              "Hasler–Hothorn", "Tamhane T2", "Dunnett T3", "Games–Howell",
              "Steel–Dwass", "Steel–Dwass Šidák",
              "Dunn + Bonferroni", "Dunn + Holm"))

## which block each method belongs to
method_group <- c(
  "Bonferroni" = "Corrections to raw p-values", "Holm" = "Corrections to raw p-values",
  "Šidák" = "Corrections to raw p-values", "Holm–Šidák" = "Corrections to raw p-values",
  "Dunnett" = "Equal variances, normal", "Holm–Dunnett" = "Equal variances, normal",
  "Tukey" = "Equal variances, normal", "Holm–Tukey" = "Equal variances, normal",
  "Hasler–Hothorn" = "Unequal variances", "Hasler–Hothorn Šidák" = "Unequal variances",
  "Tamhane T2" = "Unequal variances", "Dunnett T3" = "Unequal variances",
  "Games–Howell" = "Unequal variances",
  "Steel" = "Non-parametric", "Steel Šidák" = "Non-parametric",
  "Steel–Dwass" = "Non-parametric", "Steel–Dwass Šidák" = "Non-parametric",
  "Dunn + Bonferroni" = "Dunn (not a max-t procedure)",
  "Dunn + Holm" = "Dunn (not a max-t procedure)")

method_cols <- c(
  "Bonferroni" = crukdark, "Holm" = cruklight,
  "Šidák" = "#5B8C00", "Holm–Šidák" = "#7B2D8E",
  "Dunnett" = crukpink, "Holm–Dunnett" = "#B3006B",
  "Tukey" = crukpink, "Holm–Tukey" = "#B3006B",
  "Hasler–Hothorn" = "#E08A00", "Hasler–Hothorn Šidák" = "#B36E00",
  "Tamhane T2" = "#B36E00", "Dunnett T3" = "#8A5400",
  "Games–Howell" = "#F2B43C",
  "Steel" = "#00897B", "Steel Šidák" = "#00635A",
  "Steel–Dwass" = "#00897B", "Steel–Dwass Šidák" = "#00635A",
  "Dunn + Bonferroni" = "#6D4C41", "Dunn + Holm" = "#A1887F")

## one point shape per block, so the legend groups the methods visually
group_shapes <- c("Corrections to raw p-values" = 16, "Equal variances, normal" = 17,
                  "Unequal variances" = 15, "Non-parametric" = 18,
                  "Dunn (not a max-t procedure)" = 8)
method_shapes <- group_shapes[method_group]
names(method_shapes) <- names(method_group)

## Studentised maximum modulus: P(max|T_i| <= t) for m independent normals
## sharing one chi denominator with df nu, i.e. E_S[(2*Phi(t*S) - 1)^m].
## qmvt is far too slow to call per experiment, so this 1-D integral is
## splined over df once per design.
smm_cdf <- function(t, m, nu)
  integrate(function(u) (2 * pnorm(t * u) - 1)^m *
              2 * u * dchisq(u^2 * nu, df = nu) * nu,
            0, Inf, rel.tol = 1e-8)$value

smm_spline <- function(m, alpha) {
  grid <- c(2, 3, 4, 5, 6, 8, 10, 14, 20, 30, 50, 100, 300, 2000)
  q <- vapply(grid, function(nu)
    uniroot(function(t) smm_cdf(t, m, nu) - (1 - alpha),
            c(1, 60), tol = 1e-6)$root, numeric(1))
  f <- splinefun(log(grid), q)
  function(nu) f(log(pmin(pmax(nu, 2), 2000)))
}

## Dunn: mean ranks over ALL N observations, standardised under the global null
dunn_p <- function(y, g, pairs) {
  r  <- rank(y); N <- length(y)
  Rb <- tapply(r, g, mean); nn <- tapply(r, g, length)
  z  <- apply(pairs, 1, function(ij)
    (Rb[ij[2]] - Rb[ij[1]]) /
      sqrt(N * (N + 1) / 12 * (1 / nn[ij[1]] + 1 / nn[ij[2]])))
  2 * pnorm(abs(as.numeric(z)), lower.tail = FALSE)
}

## ---------------------------------------------------------------------------
## What each method is: what it is for, what the
## statistic is, and whether the critical value comes from each statistic on
## its own (marginal) or from the joint distribution D of all of them (joint).
## ---------------------------------------------------------------------------
method_info <- list(
  "Bonferroni" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "Yes: the threshold is a bound on the max-t value that holds whatever the correlations, so none is assumed",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Marginally",
    crit  = "each p-value must beat alpha divided by the number of comparisons",
    side  = "No. Bonferroni and Holm assume nothing about dependence, so one-sided tests are equally valid.",
    note  = "Always valid, but the least powerful of all: every other test on the same statistics (Holm, Šidák, Holm–Šidák, Dunnett, Tukey and their step-down versions) rejects at least as much."),
  "Holm" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "Yes, step-down: at each step the threshold is the Bonferroni bound on the max-t value for the hypotheses still in play",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Marginally",
    crit  = "the smallest p-value must beat alpha divided by the number of comparisons, the next one alpha divided by one fewer, and so on until one fails",
    side  = "No. Bonferroni and Holm assume nothing about dependence, so one-sided tests are equally valid.",
    note  = "Always at least as powerful as Bonferroni, with the same (no) assumptions."),
  "Šidák" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value computed as if the correlations were zero",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Marginally",
    crit  = "each p-value must beat the threshold that would give exactly alpha if the comparisons were independent",
    side  = "Yes, for the all-pairs family. Two-sided it is valid whatever the correlation; one-sided it would need every correlation positive.",
    note  = "The comparisons here are correlated, so it is conservative rather than exact."),
  "Holm–Šidák" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "Yes, step-down: at each step the threshold is the Šidák max-t value, as if the correlations were zero, for the hypotheses still in play",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Marginally",
    crit  = "the Šidák threshold for the whole family, then the one for a family of one fewer, and so on",
    side  = "Yes, as for Šidák.",
    note  = "The most powerful of the four corrections that ignore the correlations."),
  "Dunnett" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "Each arm against the control.",
    maxt  = "Yes: the threshold is the max-t value read from D, using its actual correlations",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Jointly",
    crit  = "the max-t threshold read from D, a multivariate Student t in which every pair of comparisons is correlated by one half because they share the control",
    side  = "No. A one-sided test simply reads the one-sided threshold from the same D.",
    note  = "Exact under normality and equal spreads, for any group sizes."),
  "Holm–Dunnett" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "Each arm against the control.",
    maxt  = "Yes, step-down: at each step the threshold is the max-t value read from D for the hypotheses still in play",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Jointly",
    crit  = "the largest statistic faces the threshold from D for the whole family, the next one the threshold for a family of one fewer",
    side  = "No, as for Dunnett.",
    note  = "Valid because the pooled spread does not depend on the group means."),
  "Tukey" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value read from D, using its actual correlations",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Jointly",
    crit  = "the studentised range value, which with equal group sizes is exactly the max-t threshold read from D",
    side  = "No. A one-sided version takes its threshold from D instead of the range.",
    note  = "Exact under normality and equal spreads with equal group sizes; with unequal sizes (the Tukey&ndash;Kramer version) it is slightly conservative."),
  "Holm–Tukey" = list(
    model = "Normal data, with the same spread in every group.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes, step-down: at each step the threshold is the smaller of Tukey&#39;s value and the Šidák value for the comparisons still in play",
    stat  = "The gap between the two group means, divided by its standard error, with one spread pooled over all the groups. Taken together these statistics follow a joint law D.",
    used  = "Jointly",
    crit  = "the smaller of Tukey&#39;s critical value, read from D, and the Šidák value for the number of comparisons still in play",
    side  = "No, as for Tukey.",
    note  = "At least as powerful as both Tukey and Holm&ndash;Šidák, at no cost in assumptions. Unlike Dunnett&#39;s, the exact max-t threshold for a subset of the pairs depends on which pairs are left, not only on how many, so this uses two bounds that hold for every subset. With many groups the Šidák bound is the larger of the two at most steps, and the gain over single-step Tukey is then small."),
  "Hasler–Hothorn" = list(
    model = "Normal data, but each group may have its own spread.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value read from D, with the correlations taken from the group spreads",
    stat  = "The gap between the two group means, judged against the spread of those two groups only. Their joint law D is approximately a multivariate Student t, with the correlations determined by the group spreads.",
    used  = "Jointly",
    crit  = "the max-t threshold read from D, with the correlations taken from the group spreads",
    side  = "No. A one-sided test simply reads the one-sided threshold from the same D.",
    note  = "Approximate on four counts. (1) Each comparison&#39;s degrees of freedom are themselves an approximation. (2) Each comparison has its own, but D allows only one, so a single compromise value (here the smallest) serves all. (3) Each comparison has its own denominator, while D treats them all as sharing one; this is why it can sit slightly <b>above</b> α (about 0.055 for five groups of ten) even when the spreads are equal. (4) In practice the correlations and the degrees of freedom are estimated from the data; to keep the simulation fast this app takes them from the true spreads, so the table shows a best case."),
  "Hasler–Hothorn Šidák" = list(
    model = "Normal data, but each group may have its own spread.",
    hyp   = "Each arm against the control.",
    maxt  = "Yes: the threshold is the max-t value computed as if the correlations were zero",
    stat  = "The gap between the two group means, judged against the spread of those two groups only. Their joint law D is approximately a multivariate Student t, with the correlations determined by the group spreads.",
    used  = "Marginally",
    crit  = "the Šidák threshold applied to each comparison&#39;s own p-value",
    side  = "Two-sided is always safe. One-sided is fine here too, provided every arm is tested in the same direction, since each comparison then shares the control with the same sign.",
    note  = "Conservative: it ignores that the comparisons share the control."),
  "Tamhane T2" = list(
    model = "Normal data, but each group may have its own spread.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value computed as if the correlations were zero",
    stat  = "The gap between the two group means, judged against the spread of those two groups only. Their joint law D is approximately a multivariate Student t, with the correlations determined by the group spreads.",
    used  = "Marginally",
    crit  = "the Šidák threshold applied to each comparison&#39;s own p-value",
    side  = "Yes. One-sided it would need every correlation positive.",
    note  = "Reliably conservative: the safe choice when the spreads are badly unequal and the groups are small."),
  "Dunnett T3" = list(
    model = "Normal data, but each group may have its own spread.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value for independent statistics (the studentised maximum modulus), i.e. as if the correlations were zero",
    stat  = "The gap between the two group means, judged against the spread of those two groups only. Their joint law D is approximately a multivariate Student t, with the correlations determined by the group spreads.",
    used  = "Marginally",
    crit  = "the threshold for the largest of several comparisons whose means are independent but which share one variance estimate (the studentised maximum modulus); that sharing makes it a little lower than T2&#39;s, especially when the groups are small",
    side  = "Yes, as for T2.",
    note  = "Slightly less conservative than T2; the two give the same answer once the groups are large."),
  "Games–Howell" = list(
    model = "Normal data, but each group may have its own spread.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes, approximately: the threshold is Tukey&#39;s range value, the max-t value for equal spreads",
    stat  = "The gap between the two group means, judged against the spread of those two groups only. Their joint law D is approximately a multivariate Student t, with the correlations determined by the group spreads.",
    used  = "Jointly",
    crit  = "Tukey&#39;s studentised range value, which was derived for equal spreads, so here it only approximates the right threshold",
    side  = "The range is two-sided by construction: there is no standard one-sided version.",
    note  = "More powerful than T2 and T3, but mildly liberal when the groups are small."),
  "Steel" = list(
    model = "No normality needed. Under each pair&#39;s own null the two compared groups share one distribution.",
    hyp   = "Each arm against the control.",
    maxt  = "Yes: the threshold is the max-t value read from D, using its actual correlations",
    stat  = "A standardised Wilcoxon rank sum, ranked within the compared pair only. Their joint law D is close to multivariate normal once the groups are large enough.",
    used  = "Jointly",
    crit  = "the max-t threshold read from D",
    side  = "No. A one-sided test simply reads the one-sided threshold from the same D.",
    note  = "Ranking inside the pair ignores the other groups, so each statistic stays valid whatever they do."),
  "Steel Šidák" = list(
    model = "No normality needed. Under each pair&#39;s own null the two compared groups share one distribution.",
    hyp   = "Each arm against the control.",
    maxt  = "Yes: the threshold is the max-t value computed as if the correlations were zero",
    stat  = "A standardised Wilcoxon rank sum, ranked within the compared pair only. Their joint law D is close to multivariate normal once the groups are large enough.",
    used  = "Marginally",
    crit  = "the Šidák threshold applied to each comparison&#39;s own p-value",
    side  = "Two-sided is always safe. One-sided is fine here too, provided every arm is tested in the same direction.",
    note  = "Conservative."),
  "Steel–Dwass" = list(
    model = "No normality needed. Under each pair&#39;s own null the two compared groups share one distribution.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value read from D, using its actual correlations",
    stat  = "A standardised Wilcoxon rank sum, ranked within the compared pair only. Their joint law D is close to multivariate normal once the groups are large enough.",
    used  = "Jointly",
    crit  = "the max-t threshold read from D",
    side  = "No. A one-sided test reads the one-sided threshold from the same D, whatever the signs of the correlations.",
    note  = "Ranking inside the pair ignores the other groups, so each statistic stays valid whatever they do."),
  "Steel–Dwass Šidák" = list(
    model = "No normality needed. Under each pair&#39;s own null the two compared groups share one distribution.",
    hyp   = "Every pair of groups.",
    maxt  = "Yes: the threshold is the max-t value computed as if the correlations were zero",
    stat  = "A standardised Wilcoxon rank sum, ranked within the compared pair only. Their joint law D is close to multivariate normal once the groups are large enough.",
    used  = "Marginally",
    crit  = "the Šidák threshold applied to each comparison&#39;s own p-value",
    side  = "Yes. One-sided it would need every correlation positive.",
    note  = "Conservative, and valid despite the negative correlations because the test is two-sided."),
  "Dunn + Bonferroni" = list(
    model = "No normality needed, but all the groups must share one distribution.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "No: its standard errors are right only when all the groups are alike, so D is known only under the global null and no max-t threshold can be read from it",
    stat  = "The difference in mean ranks over all N observations, standardised as if all the groups were alike. Their joint law D could be derived, but only its margins are used.",
    used  = "Marginally",
    crit  = "none is read from D: each p-value is corrected on its own",
    side  = "No. Bonferroni and Holm assume nothing about dependence, so one-sided tests are equally valid.",
    note  = "The pooled ranking makes each pair depend on the other groups, so the standard error is right only if all the groups are alike. <b>Equal distributions are what FWER control requires</b>, and Bonferroni cannot repair this. Steel (arms against a control) and Steel–Dwass (every pair) rank within each pair and avoid it.",
    corr  = "Bonferroni"),
  "Dunn + Holm" = list(
    model = "No normality needed, but all the groups must share one distribution.",
    hyp   = "The comparisons chosen on the left: each arm against the control, or every pair of groups.",
    maxt  = "No: its standard errors are right only when all the groups are alike, so D is known only under the global null and no max-t threshold can be read from it",
    stat  = "The difference in mean ranks over all N observations, standardised as if all the groups were alike. Their joint law D could be derived, but only its margins are used.",
    used  = "Marginally",
    crit  = "none is read from D; the correction is applied to the raw p-values",
    side  = "No. Bonferroni and Holm assume nothing about dependence, so one-sided tests are equally valid.",
    note  = "More powerful than Dunn + Bonferroni, but built on the same standard error. <b>Equal distributions are what FWER control requires</b>, and Holm cannot repair this.",
    corr  = "Holm"))

## rejections of every method for one simulated experiment
one_experiment <- function(k, n_arm, d, n_eff, sd_ratio, skew, family, alpha, cv) {
  dat   <- make_data(k, n_arm, d, n_eff, sd_ratio, skew)
  fp    <- family_pairs(k, family, dat$mu)
  pairs <- fp$pairs
  y <- dat$y; g <- dat$g

  ## pooled-variance contrasts
  means <- tapply(y, g, mean)
  s2    <- sum(tapply(y, g, function(v) sum((v - mean(v))^2))) / cv$nu
  se    <- sqrt(s2 * 2 / n_arm)
  tt    <- (means[pairs[, 2]] - means[pairs[, 1]]) / se
  absT  <- abs(as.numeric(tt))
  p_t   <- 2 * pt(absT, df = cv$nu, lower.tail = FALSE)

  type <- if (family == "dunnett") "Dunnett" else "Tukey"
  out <- list(
    "Bonferroni"                = p.adjust(p_t, "bonferroni") <= alpha,
    "Holm"                      = p.adjust(p_t, "holm") <= alpha,
    "Šidák"            = p_sidak(p_t) <= alpha,
    "Holm–Šidák"   = p_holm_sidak(p_t) <= alpha)
  out[[type]] <- absT > cv$single
  out[[paste0("Holm–", type)]] <- stepdown_reject(absT, cv$stepdown)

  ## Welch-based: each comparison keeps its own two variances and its own df
  w    <- t(apply(pairs, 1, function(ij) welch_pair(y, g, ij[1], ij[2])))
  absW <- abs(w[, "t"]); dfW <- w[, "df"]
  p_w  <- 2 * pt(absW, df = dfW, lower.tail = FALSE)

  out[["Hasler–Hothorn"]] <- absW > cv$hh

  if (family == "dunnett") {
    out[["Hasler–Hothorn Šidák"]] <- p_sidak(p_w) <= alpha
  } else {
    out[["Tamhane T2"]]   <- p_sidak(p_w) <= alpha
    out[["Dunnett T3"]]   <- absW > cv$smm(dfW)
    out[["Games–Howell"]] <-
      ptukey(sqrt(2) * absW, nmeans = k, df = dfW, lower.tail = FALSE) <= alpha
  }

  ## rank-based, ranked within each pair (Steel / Steel-Dwass)
  z    <- apply(pairs, 1, function(ij) rank_pair(y, g, ij[1], ij[2]))
  p_z  <- 2 * pnorm(abs(z), lower.tail = FALSE)
  rnm  <- if (family == "dunnett") "Steel" else "Steel–Dwass"
  out[[rnm]] <- abs(z) > cv$rank_single
  out[[paste(rnm, "Šidák")]] <- p_sidak(p_z) <= alpha

  ## Dunn: pooled ranks over ALL N observations, then Bonferroni or Holm.
  ## The ranks are pooled over every group, including those the family does
  ## not compare -- that is exactly what separates Dunn from Steel.
  p_d <- dunn_p(y, g, pairs)
  out[["Dunn + Bonferroni"]] <- p.adjust(p_d, "bonferroni") <= alpha
  out[["Dunn + Holm"]]       <- p.adjust(p_d, "holm") <= alpha

  list(rej = out, is_null = fp$is_null)
}

simulate_mc <- function(k, n_arm, d, n_eff, sd_ratio, skew, family, alpha,
                        n_runs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cv <- critical_values(k, n_arm, family, alpha, sd_ratio, n_eff)

  runs <- lapply(seq_len(n_runs), function(i)
    one_experiment(k, n_arm, d, n_eff, sd_ratio, skew, family, alpha, cv))

  methods <- names(runs[[1]]$rej)
  is_null <- runs[[1]]$is_null
  m1 <- sum(!is_null)

  summ <- do.call(rbind, lapply(methods, function(mm) {
    V <- vapply(runs, function(r) sum(r$rej[[mm]][r$is_null]), numeric(1))
    S <- vapply(runs, function(r) sum(r$rej[[mm]][!r$is_null]), numeric(1))
    R <- V + S
    ## Q = V/R is the share of the reported comparisons that are false, with the
    ## usual convention Q = 0 when nothing is reported; the FDR is its mean
    Q <- ifelse(R > 0, V / R, 0)
    data.frame(Method = mm, Group = unname(method_group[mm]),
               FWER = mean(V > 0),
               Power = if (m1 > 0) mean(S / m1) else NA_real_,
               AnyPower = if (m1 > 0) mean(S > 0) else NA_real_,
               FalsePositives = mean(V), TruePositives = mean(S),
               FDR = mean(Q))
  }))
  summ$Method <- factor(summ$Method, levels = methods)
  summ <- summ[, c("Group", "Method", "FWER", "Power", "AnyPower",
                   "FalsePositives", "TruePositives", "FDR")]

  list(summary = summ,
       settings = data.frame(k = k, n_arm = n_arm, comparisons = cv$m,
                             effect = d, n_eff = n_eff, sd_ratio = sd_ratio,
                             skew = skew, family = family, alpha = alpha,
                             runs = n_runs, true_nulls = sum(is_null)))
}

## ---------------------------------------------------------------------------
steps <- list(
  list(id = "s1", title = "1. Using the correlation between comparisons",
       settings = list(k = 5, n_eff = 0, sd_ratio = 1, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("No group differs: every comparison is a true null."),
         tags$li("Compare the ", tags$b("FWER"), " of Bonferroni with that of Dunnett or Tukey.")),
       answer = tags$div(
         tags$p(tags$b("Bonferroni, Dunnett and Tukey all control the FWER, but not equally
                       tightly.")),
         tags$p("Bonferroni and Holm treat the comparisons as if nothing were known
                about them. Dunnett and Tukey use the fact that the comparisons share
                groups, so their statistics are correlated — +1/2 between any two
                Dunnett contrasts in a balanced design — and the critical value is
                read off that joint distribution."),
         tags$p("Because they use the real joint distribution, Dunnett and Tukey spend
                the whole α budget (exactly α in theory), while Bonferroni sits below
                it: about 0.046 against 0.038 against a control, and 0.047 against 0.035
                with every pair. Each simulated FWER wobbles by a few thousandths, but
                the gap does not: Bonferroni rejects a subset of what Dunnett or Tukey
                rejects in the very same experiment. That unused margin is wasted ",
                tags$b("statistical power", .noWS = "outside"), ", which the next step
                makes visible."),
         tags$p(tags$b("One method is not below α here:"), " Hasler–Hothorn, with every
                pair, sits at about 0.056. It gives each comparison its own variance
                estimate, but reads its threshold from a joint distribution that assumes
                they all share one. That approximation is mildly liberal with small
                groups — the price of not assuming equal spreads."))),

  list(id = "s2", title = "2. The statistical power that correlation buys",
       settings = list(k = 5, n_eff = 2, sd_ratio = 1, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("Two of the groups now carry a real effect."),
         tags$li("Compare the ", tags$b("Statistical power"), " of Dunnett (or Tukey)
                  with Bonferroni and Šidák: all three are single-step."),
         tags$li("Then compare Holm–Dunnett (or Holm–Tukey) with Holm and Holm–Šidák:
                  all step-down.")),
       answer = tags$div(
         tags$p(tags$b("Two separate gains, and they combine.")),
         tags$p("The first gain is the joint distribution: Dunnett and Tukey are more
                powerful than Bonferroni and Šidák because they stop correcting
                for multiplicity that the correlation has already removed."),
         tags$p("The second is the step-down. Holm–Dunnett and Holm–Tukey apply the
                same max-t idea sequentially, relaxing the critical value as
                hypotheses are rejected. They reject everything the single-step
                version rejects and sometimes more — again for free. Against a control
                the two gains add up almost exactly (+0.021 each, +0.042 together); with
                every pair they overlap, so together they give less than their sum."),
         tags$p("Compare like with like. Against a control, Dunnett reaches about 0.79
                against 0.77 for Bonferroni, and Holm–Dunnett about 0.81 against 0.79 for
                Holm; with every pair, Tukey about 0.69 against 0.65, and Holm–Tukey
                about 0.71 against 0.70. Holm already gains from stepping down, so it can
                match or even beat single-step Dunnett or Tukey — here it matches Dunnett
                and beats Tukey. It never beats Holm–Dunnett or Holm–Tukey, which use a
                threshold at least as low as Holm's at every step. With five groups the
                gains are a few hundredths; step 5 shows what matters more when there are
                many groups: the choice of the family of hypotheses."))),

  list(id = "s3", title = "3. When the variances differ",
       ## SD 0.25 for the groups with an effect: the TRUE-null comparisons then
       ## involve only SD-1 groups, the pooled variance is too small for them, and
       ## the pooled methods break above alpha in both families. (With SD 4 against
       ## a control they only become very conservative, which hides the lesson.)
       settings = list(k = 5, n_eff = 2, sd_ratio = 0.25, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("The two groups carrying an effect now have SD 0.25; the other three,
                  the control included, keep SD 1."),
         tags$li("Watch the ", tags$b("FWER"), " of Dunnett or Tukey, and of Bonferroni
                  and Holm too."),
         tags$li("Compare them with the Welch-based methods: Hasler–Hothorn, and with
                  every pair of groups also T2, T3 and Games–Howell.")),
       answer = tags$div(
         tags$p(tags$b("Equal variances is a real assumption, not a formality.")),
         tags$p("Dunnett and Tukey pool all the groups into one variance estimate. Here
                the two quiet groups pull it down, to 0.625 against a true variance of 1.
                Against a control, every TRUE-null comparison is between two SD-1 groups;
                with every pair, three of the four are (the fourth, the two quiet groups
                against each other, gets a standard error three times too large and is
                almost never rejected). For those SD-1 comparisons the standard error is
                too small, the t statistics too large, and the FWER rises clearly above α."),
         tags$p(tags$b("Bonferroni and Holm fail too."), " They assume nothing about the
                dependence between the comparisons, but they correct p-values computed
                from the same pooled t statistics: when the raw p-values are wrong, no
                correction can repair them. The ", tags$b("step-down"), " methods suffer
                most — Holm, Holm–Šidák and Holm–Dunnett or Holm–Tukey all reach 0.10 to
                0.14 — because once they have rejected the real effects they relax the
                threshold for the remaining TRUE-null comparisons, whose t statistics are
                inflated."),
         tags$p("The Welch-based procedures give every comparison its own variance
                estimate and its own Satterthwaite degrees of freedom, and stay below α.
                They cost ", tags$b("statistical power"), " when the variances really
                are equal (in step 2, Hasler–Hothorn about 0.75 against Dunnett 0.79,
                and T2 about 0.56 against Tukey 0.69), and they are the safe choice when
                you do not know."),
         tags$p("Now move the SD slider above 1, say to 4, so that the groups with an
                effect become the noisy ones. Against a control, the TRUE-null
                comparisons then involve only quiet groups, the pooled variance is far
                too large for them, and every pooled method becomes very conservative:
                FWER 0.000. Yet the same pooled variance is too ", tags$i("small"), " for
                the comparisons that do contain a noisy group, so the pooled methods end
                up with more ", tags$b("statistical power"), " than the Welch-based ones
                (Dunnett about 0.14 against Hasler–Hothorn 0.07). A wrong standard error
                does not simply make a test cautious: it distorts both error rates, in
                directions you cannot predict. With every pair of groups the pair of the
                two noisy groups is a TRUE null judged with too small a variance, and the
                FWER again exceeds α."))),

  list(id = "s4", title = "4. When the data are not normal",
       ## 20 per group and d = 0.7: with 10 per group the rank tests' critical values
       ## are conservative and their gain on skewed data, with every pair, is hidden
       settings = list(k = 5, n_eff = 2, sd_ratio = 1, skew = TRUE, d = 0.7, n_arm = 20),
       look = tags$ul(
         tags$li("The data are now strongly right-skewed rather than normal; 20
                  observations per group, and the groups with an effect have a mean
                  0.7 above the control's."),
         tags$li("Compare the rank-based methods (Steel or Steel–Dwass, and Dunn) with
                  the t-based ones, on both ", tags$b("FWER"), " and ",
                 tags$b("Statistical power", .noWS = "outside"), "."),
         tags$li("Then untick ", tags$i("Right-skewed data"), " and run again.")),
       answer = tags$div(
         tags$p(tags$b("Skewness does not break the t-based tests here; it is where the
                       rank-based tests pay off.")),
         tags$p("With equal group sizes and equal spreads the t statistic is robust to
                skewness: the FWER of Dunnett and Tukey is the same as with normal data,
                and so is their ", tags$b("statistical power", .noWS = "outside"),
                " (it is even slightly higher)."),
         tags$p("The rank-based procedures — Steel (arms against a control) or
                Steel–Dwass (every pair), and Dunn — make no normality assumption, and
                on these skewed data they are by far the most powerful in the table:
                about 0.64 for Steel against 0.42 for Dunnett, and 0.53 for Steel–Dwass
                against 0.31 for Tukey. With normal data the order reverses (Steel about
                0.35 against Dunnett 0.40; Steel–Dwass 0.24 against Tukey 0.29)."),
         tags$p("The reason is not that the mean becomes unreliable: its variance is the
                spread squared divided by the group size whatever the shape, which is why
                the t-based numbers hardly move. It is that ranks become far more
                efficient here. Most of the values are bunched at the low end of the
                distribution, so a shift pushes most of the treated values past the
                control's, and that is exactly what a rank test counts. For
                exponential-shaped data the Wilcoxon test is asymptotically about three
                times as efficient as the t test; for normal data it is about 0.95 times."),
         tags$p("Steel and Steel–Dwass use the same max-t logic, with the joint
                distribution of the rank statistics; Dunn is the exception, correcting
                its p-values one at a time instead."))),

  list(id = "s5", title = "5. Many groups",
       settings = list(k = 10, n_eff = 2, sd_ratio = 1, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("Ten groups: 45 comparisons for every pair, 9 against a control."),
         tags$li("With ", tags$i("every pair of groups"), ", compare the ",
                 tags$b("Statistical power"), " of Tukey with Bonferroni."),
         tags$li("Then switch the family of hypotheses to ", tags$i("each arm against a control"),
                 " and compare Dunnett's ", tags$b("Statistical power"),
                 " with Tukey's."),
         tags$li("In both families of hypotheses, look at the ", tags$b("FDR", .noWS = "outside"), ".")),
       answer = tags$div(
         tags$p(tags$b("With many groups, choosing the family of hypotheses matters far
                       more than the correction.")),
         tags$p("Tukey still beats Bonferroni, but only a little: with ten groups most
                pairs of comparisons share no group, so there is little correlation to
                use."),
         tags$p(tags$b("Too many tests hide the real effects."), " The real effects are
                the same in both families of hypotheses, but asking for every pair means 45 comparisons
                instead of 9. The more tests you run, the higher the bar each one must
                clear, and the harder it becomes to detect a real difference: ",
                tags$b("statistical power"), " falls. That is why Dunnett finds clearly
                more of the real effects than Tukey. Test only the comparisons you need:
                if you only care about the comparisons against the control, do not ask
                for every pair."),
         tags$p("The FDR is very low, especially with every pair: when a method reports
                many real differences, an occasional false one is only a small share of
                its list.")))
)
step_titles <- vapply(steps, function(s) s$title, character(1))

## ---------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML(sprintf("
    .well { background:#F7FAFF; border-color:%s33; }
    h4 { color:%s; }
    .answer-box { border:2px solid %s; border-radius:5px; padding:10px 14px;
                  background:#FFF; margin-top:10px; }
    .step-box { border-left:5px solid %s; padding:6px 12px; background:#F2F9FF; }
    .nav-pills { display:flex; flex-wrap:wrap; gap:5px; }
    .nav-pills::before, .nav-pills::after { display:none; }
    .nav-pills > li, .nav-pills > li + li { float:none; margin:0;
                                            flex:0 0 calc(33.333%% - 4px); }
    .nav-pills > li > a { text-align:left; }
  ", crukdark, crukdark, crukpink, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "FWER and statistical power: Dunnett and Tukey")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Several groups are compared in one experiment. The <b>family of
            hypotheses</b>, the comparisons tested together, is either each arm against a
            <b>control</b> (Dunnett's family) or <b>every pair</b> (Tukey's family).
            One outcome is measured on <b>independent</b> individuals, the same number in
            every group. By default it is <b>normally distributed with equal variance</b> in
            every group; the controls on the left can give a different spread to the
            groups with an effect, or make the data right-skewed, to see which methods
            still hold. Because the comparisons share groups, their test statistics are correlated,
            and a procedure that uses that correlation needs a smaller correction than one
            that does not (hence, <b>greater statistical power</b>). The table reports how often each method makes
            <b style='color:%s'>any false claim</b> and what share of the real differences
            it finds (<b>statistical power</b>).", crukdark))),

  do.call(tabsetPanel, c(
    ## Free exploration first: with three tabs per row the pills fall into
    ## two lines, "Free, 1, 2" then "3, 4, 5"
    list(id = "guide", type = "pills", selected = "Free exploration",
         tabPanel("Free exploration",
                  div(class = "step-box",
                      HTML("Set the design on the left and press <b>Run simulation</b>.
                            Or choose one of the numbered steps.")))),
    lapply(steps, function(s)
      tabPanel(s$title,
               div(class = "step-box", h4(s$title),
                   actionButton(paste0("apply_", s$id), "Apply these settings",
                                class = "btn-sm", style = "margin-bottom:6px;"),
                   s$look))))),
  hr(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      actionButton("run", "Run simulation", class = "btn-primary",
                   style = sprintf("background:%s;border-color:%s;width:100%%;font-size:15px;",
                                   crukdark, crukdark)),
      br(), br(),
      h4("Family of hypotheses"),
      radioButtons("family", NULL,
                   choices = c("Each arm against a control (Dunnett)" = "dunnett",
                               "Every pair of groups (Tukey)" = "tukey"),
                   selected = "dunnett"),
      uiOutput("fam_text"),

      h4("Statistical design"),
      sliderInput("k", "Number of groups", min = 3, max = 10, value = 5, step = 1),
      sliderInput("n_arm", "Observations per group", min = 4, max = 30, value = 10, step = 1),
      uiOutput("n_arm_text"),
      sliderInput("n_eff", "Groups carrying a real effect", min = 0, max = 9, value = 2, step = 1),
      uiOutput("eff_text"),
      sliderInput("d", "Effect size: difference in means", min = 0.2, max = 3, value = 1.5, step = 0.1),
      uiOutput("d_text"),

      h4("Assumptions"),
      sliderInput("sd_ratio", "SD of the groups with an effect",
                  min = 0.25, max = 6, value = 1, step = 0.25),
      uiOutput("sd_text"),
      checkboxInput("skew", "Right-skewed data instead of normal", value = FALSE),
      uiOutput("skew_text"),

      h4("Simulation"),
      sliderInput("alpha", "α — the family-wise error rate",
                  min = 0.01, max = 0.10, value = 0.05, step = 0.01),
      sliderInput("runs", "Number of simulated experiments",
                  min = 100, max = 10000, value = 5000, step = 100),
      helpText("5000 experiments give FWERs to about ±0.006 and take about 6–10 seconds
                with 5 groups, but about 40 seconds for every pair of 10 groups. Lower
                the number to explore quickly, or raise it to 10 000 for sharper numbers."),
      ## a fixed seed makes a run reproducible: the same settings give the same
      ## numbers, so the values quoted in the steps can be checked exactly
      numericInput("seed", "Random seed (empty = a new draw every run)", value = 1),
      downloadButton("save", "Download CSV", style = "width:100%;")
    ),

    mainPanel(
      width = 9,
      h4(textOutput("scenario")),
      DTOutput("tab"),
      helpText(HTML("Each comparison gives a test statistic from simulated data on
                k groups. Write <b>m</b> for the number of comparisons in the family of hypotheses,
                of which <b>m<sub>0</sub></b> are TRUE nulls (the two groups
                have the same mean) and <b>m&nbsp;&minus;&nbsp;m<sub>0</sub></b> are
                FALSE. In one experiment <b>V</b> of the TRUE nulls are rejected
                (false positives), <b>S</b> of the FALSE nulls are rejected
                (true positives), and <b>R</b> = V + S comparisons are reported.<br>
                <b>FWER</b> := P(V &gt; 0), estimated by the share of experiments with
                at least one false positive.
                <b>Statistical power</b>: the average of S/(m&nbsp;&minus;&nbsp;m<sub>0</sub>), the
                share of the real differences found. <b>Any statistical power</b>: the chance that
                S &gt; 0. <b>False positives (mean)</b> and <b>True positives (mean)</b>:
                the average <i>counts</i> V and S per experiment, not proportions.
                <b>FDR</b>: of the comparisons a method reports, the average share that
                are false (counted as 0 when it reports nothing).<br>
                The first block of methods corrects the raw p-values and uses nothing
                else. The next three assume progressively less: equal variances and
                normality, then unequal variances (Welch statistics), then no normality
                (ranks within each pair). All of these are max-t procedures:
                the joint methods read their critical value from the joint distribution
                of the comparisons, correlations included, while the Šidák versions use
                the max-t value computed as if the comparisons were independent. The last
                block, Dunn, is not a max-t procedure: it ranks all the data together
                and simply corrects the resulting p-values with Bonferroni or Holm.")),
      uiOutput("method_box"),
      br(),
      ## the plot gets the full width: with 12-14 methods it needs room for the legend
      fluidRow(
        column(12, h4("Error rate against statistical power"), plotOutput("tradeoff", height = "500px"),
               uiOutput("tradeoff_note"))),
      br(),
      fluidRow(
        column(12, h4("What is going on"), uiOutput("reveal_ui"), uiOutput("answer"),
               ## breathing room so the answer box is never flush with the window edge
               div(style = "height:60px;"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$k, {
    updateSliderInput(session, "n_eff", max = input$k - 1,
                      value = min(input$n_eff, input$k - 1))
  })

  ## spell out what the SD ratio does, with the group SDs it actually produces
  ## the effect size is measured in the CONTROL's SDs, which the SD slider does
  ## not touch, so the two sliders stay independent
  output$d_text <- renderUI({
    ne <- min(input$n_eff, input$k - 1)
    if (ne == 0)
      return(helpText("The difference between the mean of a group with an effect and
                       the control's mean. No group carries an effect at the moment, so
                       this slider does nothing."))
    base <- sprintf("The difference between the mean of each group with an effect and
                     the control&#39;s mean: at %s, their mean is <b>%s</b> against the
                     control&#39;s 0.",
                    as.character(input$d), as.character(input$d))
    helpText(HTML(base))
  })

  output$sd_text <- renderUI({
    ne <- min(input$n_eff, input$k - 1)
    helpText(HTML(
      if (ne == 0)
        "Every group without an effect has SD 1. No group carries an effect at the
         moment, so all the SDs are 1 and this slider does nothing: raise
         <i>Groups carrying a real effect</i> to make the variances differ."
      else if (input$sd_ratio == 1)
        sprintf("At 1 every group has SD 1, so the equal-variance assumption holds.
                 Move it and the %d group%s carrying an effect take that SD while the
                 other %d, the control included, stay at 1.",
                ne, if (ne == 1) "" else "s", input$k - ne)
      else
        sprintf("The %d group%s carrying an effect %s SD <b>%s</b>; the other %d,
                 the control included, stay at <b>1</b>. So the group%s with an effect
                 %s <b>%s variable</b> than the rest.",
                ne, if (ne == 1) "" else "s", if (ne == 1) "has" else "share",
                as.character(input$sd_ratio), input$k - ne,
                if (ne == 1) "" else "s", if (ne == 1) "is" else "are",
                if (input$sd_ratio > 1) "more" else "less")))
  })

  ## spell out what the skew box does, and when it bites
  ## the rank statistics are bounded: with very few observations the largest
  ## possible |z| can fall below the max-t threshold, so no rejection is possible
  output$n_arm_text <- renderUI({
    n <- input$n_arm
    zmax <- (n / 2) * sqrt(12 / (2 * n + 1))
    cv <- tryCatch(critical_values(input$k, n, input$family, input$alpha)$rank_single,
                   error = function(e) NA_real_)
    if (!is.na(cv) && zmax <= cv)
      div(style = sprintf("color:%s; font-weight:bold; margin-top:4px;", crukpink),
          sprintf("With only %d observations per group the largest possible rank statistic is %.2f, below the threshold of %.2f: Steel%s can never reject here. Raise the group size.",
                  n, zmax, cv, if (input$family == "dunnett") "" else "–Dwass"))
    else helpText("Every group has the same number of observations.")
  })

  output$skew_text <- renderUI({
    helpText(HTML(if (!isTRUE(input$skew))
      "Unticked: the data are <b>normal</b>, symmetric around each group mean.
       Ticked: they become <b>right-skewed</b>: most values bunched low, with a long
       tail of occasional large ones, as for biomarker concentrations or cell counts
       (exponential distribution). Only the shape changes; the means and SDs stay as
       set above."
    else
      "Now <b>right-skewed</b>: most values bunched low, with a long tail of occasional
       large ones, as for biomarker concentrations or cell counts (exponential
       distribution). The means and SDs are unchanged. The t-based methods assume
       normal data; the rank-based ones (Steel, Steel&ndash;Dwass, Dunn) do not.",
      if (isTRUE(input$skew) && input$sd_ratio != 1 && min(input$n_eff, input$k - 1) > 0)
        " <b>With skewed data and unequal spreads a group with an effect no longer differs
          from the control by a simple shift</b>: the means differ by the effect size, but the
          rest of the distribution does not move with them. The rank-based tests then answer
          a different question — whether one group tends to give the larger value — so
          compare their statistical power with care."
      else ""))
  })

  ## every effect group is shifted by the same amount, in the same direction
  output$eff_text <- renderUI({
    k <- input$k; ne <- min(input$n_eff, k - 1)
    if (ne == 0)
      return(helpText("No group carries an effect: every comparison is a TRUE null."))
    helpText(HTML(sprintf(
      "The last %d group%s all have their mean shifted by the <b>same amount</b> (the
       effect size set below) and in the <b>same direction</b> (upwards); every other
       group, the control included, has mean 0. These are also the groups whose SD the
       <i>SD of the groups with an effect</i> slider changes.",
      ne, if (ne == 1) "" else "s")))
  })

  output$fam_text <- renderUI({
    m <- if (input$family == "dunnett") input$k - 1 else choose(input$k, 2)
    helpText(sprintf("%d comparisons among %d groups.", m, input$k))
  })

  lapply(steps, function(s) {
    observeEvent(input[[paste0("apply_", s$id)]], {
      st <- s$settings
      updateSliderInput(session, "k", value = st$k)
      updateSliderInput(session, "n_eff", value = st$n_eff)
      updateSliderInput(session, "sd_ratio", value = st$sd_ratio)
      updateCheckboxInput(session, "skew", value = st$skew)
      if (!is.null(st$d)) updateSliderInput(session, "d", value = st$d)
      updateSliderInput(session, "n_arm", value = if (is.null(st$n_arm)) 10 else st$n_arm)
      ## the numbers quoted in the steps are at alpha = 0.05
      updateSliderInput(session, "alpha", value = 0.05)
      showNotification("Settings applied — now press 'Run simulation'.",
                       type = "message", duration = 3)
    }, ignoreInit = TRUE)
  })

  ## switching the family changes which methods exist, so it re-runs at once:
  ## otherwise the table kept showing the other family's methods until Run
  sim <- eventReactive(list(input$run, input$family), ignoreNULL = FALSE, {
    withProgress(message = "Simulating experiments...", value = 0.4, {
      simulate_mc(k = input$k, n_arm = input$n_arm, d = input$d,
                  n_eff = min(input$n_eff, input$k - 1),
                  sd_ratio = input$sd_ratio, skew = input$skew,
                  family = input$family, alpha = input$alpha, n_runs = input$runs,
                  seed = if (is.null(input$seed) || is.na(input$seed)) NULL else input$seed)
    })
  })

  output$scenario <- renderText({
    s <- sim()$settings
    sprintf("%d groups of %d | %s: %d comparisons, %d of them true nulls | %s | %s | %d experiments | α = %.2f",
            s$k, s$n_arm,
            if (s$family == "dunnett") "against a control" else "all pairs",
            s$comparisons, s$true_nulls,
            if (s$sd_ratio == 1) "equal variances"
            else sprintf("SD %s in the groups with an effect", format(s$sd_ratio)),
            if (s$skew) "skewed data" else "normal data",
            s$runs, s$alpha)
  })

  picked <- reactiveVal(NULL)
  observeEvent(input$pick_method, picked(input$pick_method))
  ## a method that the other family does not run should not stay on screen
  observeEvent(input$family, {
    if (!is.null(picked()) && !(picked() %in% method_sets[[input$family]]))
      picked(NULL)
  })

  output$method_box <- renderUI({
    mm <- picked()
    if (is.null(mm))
      return(div(class = "answer-box",
                 HTML("<b>Click the name of any test above.</b> Each one is described with the same
                       points used in the lecture: the <b>statistical model</b> it assumes, the
                       <b>hypotheses</b> it tests, whether it is a <b>max-<i>t</i> procedure</b>, its
                       <b>test statistics and their joint law</b> <i>D</i>, and <b>how that joint law is
                       used</b> — <b>jointly</b>, reading one threshold from <i>D</i> and so using
                       the correlations, or <b>marginally</b>, taking each comparison on its own.
                       <i>D</i> is the distribution the whole set of test statistics follows together,
                       correlations included.")))
    it <- method_info[[mm]]
    if (is.null(it))
      return(div(class = "answer-box", HTML(sprintf("<b>%s</b>", mm))))
    div(class = "answer-box",
        h4(mm, style = sprintf("margin-top:0; color:%s;", crukdark)),
        HTML(sprintf(
          "<b>1 Statistical model.</b> %s<br>
           <b>2 Hypotheses.</b> %s<br>
           <b>3 Max-<i>t</i> procedure.</b> <span style='color:%s'><b>%s</b></span><br>
           <b>4 Test statistics and their joint law.</b> %s<br>
           <b>5 How the joint law is used.</b> <span style='color:%s'><b>%s</b></span>: %s.<br>
           <b>6 Two-sided needed?</b> %s<br>
           <b>7 In practice.</b> %s%s",
          it$model, it$hyp, crukpink, it$maxt, it$stat, crukpink, it$used, it$crit,
          it$side, it$note,
          if (is.null(it$corr)) "" else
            sprintf("<br><b>8 Correction.</b> <span style='color:%s'><b>%s</b></span> on the raw <i>p</i>-values",
                    crukpink, it$corr))))
  })

  output$tab <- renderDT({
    d <- sim()$summary
    d$Group[duplicated(d$Group)] <- ""
    ## the method name is a link: clicking it fills the explanation box below
    d$Method <- sprintf(
      paste0('<a href="#" style="color:%s; text-decoration:underline dotted;"',
             ' onclick="Shiny.setInputValue(&quot;pick_method&quot;, &quot;%s&quot;,',
             ' {priority: &quot;event&quot;}); return false;">%s</a>'),
      crukdark, as.character(d$Method), as.character(d$Method))
    ## paging = FALSE matters: DT shows only 10 rows by default and dom = "t"
    ## hides the pager, so the last methods would silently disappear
    datatable(d, rownames = FALSE, selection = "none", escape = FALSE,
              options = list(dom = "t", ordering = FALSE, paging = FALSE),
              colnames = c("", "Method", "FWER", "Statistical power", "Any statistical power",
                           "False positives (mean)", "True positives (mean)",
                           "FDR")) |>
      formatStyle("Group", fontWeight = "bold", color = crukdark) |>
      formatRound(c("FWER", "Power", "AnyPower", "FDR"), 3) |>
      formatRound(c("FalsePositives", "TruePositives"), 2) |>
      formatStyle("FWER", color = crukdark, fontWeight = "bold") |>
      formatStyle("Power", color = crukpink, fontWeight = "bold")
  })

  ## how to read the plot, in ordinary-sized text rather than a small caption
  output$tradeoff_note <- renderUI({
    d <- sim()$summary; a <- sim()$settings$alpha
    txt <- if (all(is.na(d$Power)))
      sprintf("No group carries an effect here, so there is nothing to find and statistical power is
               not defined. Each bar is one method&#39;s FWER: the chance of at least one
               false claim. The dashed line marks &alpha; = %s; a bar above it means that
               method is <b>not</b> controlling the FWER.", a)
    else
      sprintf("Each point is one method, named in the legend; the shape shows its block in the table. <b>Left to right</b>: how often it makes at least
               one false claim (FWER). <b>Bottom to top</b>: the share of the real differences
               it finds (statistical power). The best methods sit <b>high up and to the left</b>.
               The dashed line marks &alpha; = %s: a point to its right means that method is
               <b>not</b> controlling the FWER.", a)
    div(style = "font-size:14px; color:#333; margin-top:6px;", HTML(txt))
  })

  output$tradeoff <- renderPlot({
    d <- sim()$summary
    a <- sim()$settings$alpha
    if (all(is.na(d$Power)))
      return(ggplot(d, aes(Method, FWER, fill = Method)) +
               geom_col(width = 0.7) +
               geom_hline(yintercept = a, linetype = 2, colour = "grey30") +
               scale_fill_manual(values = method_cols) +
               labs(x = NULL, y = "FWER") +
               theme_minimal(base_size = 12) +
               theme(axis.text.x = element_text(angle = 25, hjust = 1),
                     legend.position = "none"))
    ## with 12-14 methods the text labels collide, so name them in a legend;
    ## colour and shape share the title, so ggplot merges them into one key
    ggplot(d, aes(FWER, Power, colour = Method, shape = Method)) +
      geom_vline(xintercept = a, linetype = 2, colour = "grey40") +
      geom_point(size = 3.6, stroke = 1.1) +
      scale_colour_manual(values = method_cols) +
      scale_shape_manual(values = method_shapes) +
      coord_cartesian(xlim = c(0, max(0.12, max(d$FWER) * 1.15)), ylim = c(0, 1.05)) +
      labs(x = "FWER", y = "Statistical power", colour = NULL, shape = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "right", legend.text = element_text(size = 11),
            legend.key.height = unit(0.9, "lines"))
  })

  current_step <- reactive({
    ## length 0 before the tab menu has reported its value
    i <- match(input$guide, step_titles)
    if (length(i) == 0 || is.na(i)) NULL else steps[[i]]
  })

  output$reveal_ui <- renderUI({
    if (is.null(current_step()))
      helpText("Choose one of the numbered steps above. Each keeps its explanation
                hidden until you have run the simulation and read the numbers.")
    else
      actionButton("reveal", "Show the explanation", class = "btn-sm",
                   style = sprintf("background:%s;border-color:%s;color:#FFF;", crukpink, crukpink))
  })

  revealed <- reactiveVal(FALSE)
  observeEvent(input$reveal, { req(input$reveal > 0); revealed(TRUE) })
  observeEvent(input$guide, revealed(FALSE))

  output$answer <- renderUI({
    st <- current_step()
    if (is.null(st) || !revealed()) return(NULL)
    div(class = "answer-box", st$answer)
  })

  ## ---- saving ------------------------------------------------------------
  output$save <- downloadHandler(
    filename = function() sprintf("dunnett-tukey-%s.csv", format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      res <- sim()
      write.csv(cbind(res$settings[rep(1, nrow(res$summary)), ], res$summary),
                file, row.names = FALSE)
    })
}

shinyApp(ui, server)
