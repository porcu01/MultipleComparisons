## Simulation study: power of FWER corrections as the number of genes grows.
##
## For each family size N and each proportion of truly-affected genes, we
## simulate 100,000 independent trials.  The N genes are statistically
## independent.  Each gene is a two-sided two-sample
## t-test, 5 units per arm, equal variances; affected genes have |ES| = 2.
## The p-value of that test is exactly the two-sided tail of a noncentral t
## with df = 8 and ncp = ES / sqrt(2/5), so we draw p-values directly rather
## than the 10 raw normals per gene.  Null genes have uniform p-values.
##
## Reported quantity: the proportion of the ES = 2 genes that reach
## significance at alpha = 0.05 after each correction (average power).
## A significant result counts as a positive whatever its direction: a
## significant test with the opposite sign is also counted.

set.seed(20260904)

alpha  <- 0.05
ES     <- 2
n_arm  <- 5
df     <- 2 * n_arm - 2
ncp    <- ES / sqrt(2 / n_arm)
Ns     <- c(10, 100, 1000, 10000, 100000)
fracs  <- c(0.10, 0.50)
trials <- 100000

## Holm and Hochberg both compare the ordered p-values with alpha/(N - i + 1);
## Holm steps down (stop at the first failure), Hochberg steps up (stop at the
## last success).  Both reject the k smallest p-values, so we only need to know
## how many of those k belong to the affected genes.
run_cell <- function(N, frac, trials) {
  m1   <- round(N * frac)
  m0   <- N - m1
  crit <- alpha / (N:1)                      # thresholds for the ordered p-values
  bonf_thr  <- alpha / N
  sidak_thr <- 1 - (1 - alpha)^(1 / N)

  hits <- matrix(0, nrow = trials, ncol = 4,
                 dimnames = list(NULL, c("Bonferroni", "Holm", "Sidak", "Hochberg")))

  for (b in seq_len(trials)) {
    ## Any significant affected gene counts as a hit, whatever the sign of
    ## its test statistic.
    t_alt <- rt(m1, df = df, ncp = ncp)
    p_alt <- 2 * pt(-abs(t_alt), df = df)
    p     <- c(p_alt, runif(m0))
    hit   <- c(rep(TRUE, m1), rep(FALSE, m0))

    ord    <- order(p)
    ps     <- p[ord]
    alt_s  <- hit[ord]
    pass   <- ps <= crit

    ## Holm: reject up to the last index before the first failure.
    k_holm <- if (pass[1]) {
      fail <- which(!pass)
      if (length(fail)) fail[1] - 1L else N
    } else 0L

    ## Hochberg: reject up to the largest index that still passes.
    k_hoch <- {
      ok <- which(pass)
      if (length(ok)) ok[length(ok)] else 0L
    }

    cum_alt <- cumsum(alt_s)
    hits[b, "Bonferroni"] <- sum(p_alt < bonf_thr)
    hits[b, "Sidak"]      <- sum(p_alt < sidak_thr)
    hits[b, "Holm"]       <- if (k_holm > 0L) cum_alt[k_holm] else 0L
    hits[b, "Hochberg"]   <- if (k_hoch > 0L) cum_alt[k_hoch] else 0L
  }

  colMeans(hits) / m1
}

## The results are cached: the (seeded, hence reproducible) simulation only
## runs when the CSV is missing.  Delete the CSV, or set refresh = TRUE, to
## run it again.
out     <- "fwer-power-simulation.csv"
refresh <- FALSE

if (!refresh && file.exists(out)) {
  res <- read.csv(out, check.names = FALSE)
  message("Using cached results in ", out, " (delete it to re-simulate).")
} else {
  res <- do.call(rbind, lapply(fracs, function(f) {
    do.call(rbind, lapply(Ns, function(N) {
      t0  <- proc.time()[["elapsed"]]
      pow <- run_cell(N, f, trials)
      message(sprintf("frac = %.0f%%  N = %6d  done in %5.1f s", 100 * f, N,
                      proc.time()[["elapsed"]] - t0))
      data.frame(frac = f, N = N, t(pow), check.names = FALSE)
    }))
  }))
  write.csv(res, out, row.names = FALSE)
}

print(res, digits = 3)
