library(shiny)
library(shinyWidgets)
library(DT)
library(ggplot2)

crukdark  <- "#00007E"
cruklight <- "#009CEE"
crukpink  <- "#FF0087"

method_labels <- c(
  none        = "No correction",
  bonferroni  = "Bonferroni",
  holm        = "Holm",
  sidak       = "Šidák",
  BH          = "Benjamini-Hochberg",
  BY          = "Benjamini-Yekutieli",
  storey      = "Storey"
)
method_cols <- c(
  "No correction"        = "grey45",
  "Bonferroni"           = crukdark,
  "Holm"                 = cruklight,
  "Šidák"                = "#5B8C00",
  "Benjamini-Hochberg"   = crukpink,
  "Benjamini-Yekutieli"  = "#7B2D8E",
  "Storey"               = "#E08A00"
)

## a pair of genes is negatively correlated when exactly one of the two loads
## the opposite way, so flipping k of m0 genes makes k(m0 - k) of the
## m0(m0 - 1)/2 pairs negative. Only the TRUE nulls matter here: the conditions
## BH and Šidák need constrain the true-null p-values alone, so negative
## dependence among the genes that really are differentially expressed cannot
## cause a false discovery. k(m0 - k) grows in steps, so not every count is
## reachable.
flips_for_pairs <- function(m0, want) {
  k <- 0:floor(m0 / 2)
  k[which.min(abs(k * (m0 - k) - want))]
}
neg_pair_n   <- function(m0, k) k * (m0 - k)
neg_pair_pct <- function(m0, k) if (m0 < 2) 0 else 100 * k * (m0 - k) / choose(m0, 2)
max_neg_pairs <- function(m0) floor(m0 / 2) * ceiling(m0 / 2)

## Šidák: 1 - (1 - p)^N, written so that it stays accurate for tiny p
p_sidak <- function(p) pmin(1, -expm1(length(p) * log1p(-p)))

## Storey's adaptive BH: the BH adjusted p-values scaled by the estimated
## proportion of TRUE nulls. With lambda = 0.5, every p-value above 0.5 is
## counted as coming from a TRUE null; the estimate is capped at 1, so Storey
## never rejects less than BH. The +1 is Storey, Taylor and Siegmund's
## conservative form: without it the estimate can be exactly 0, and Storey
## would then reject every gene whatever its p-value.
p_storey <- function(p, lambda = 0.5) {
  n   <- length(p)
  pi0 <- min(1, (1 + sum(p > lambda)) / (n * (1 - lambda)))
  pmin(1, pi0 * p.adjust(p, "BH"))
}

## ---------------------------------------------------------------------------
## Simulation
##
## One experiment = N genes measured in two groups of n_arm samples each.
## m0 genes are TRUE nulls (no differential expression), the remaining
## m1 = N - m0 are truly differentially expressed with effect size d.
##
## Each gene is tested with a two-sample t statistic, simulated directly from
## its sampling distribution: a numerator z_i ~ N(mu_i, 1) with mu_i = 0 for a
## true null and mu_i = d * sqrt(n_arm / 2) otherwise, divided by that gene's
## own variance estimate, t_i = z_i / sqrt(w_i), w_i ~ chi^2_df / df,
## df = 2(n_arm - 1). The numerators are made dependent through one latent factor:
##
##     z_i = s_i * sqrt(|rho|) * F + sqrt(1 - |rho|) * e_i
##
## with F and e_i independent standard normals. All loadings s_i = +1 gives
## equicorrelated genes with correlation +|rho| (co-regulated genes). Taking
## s_i = +1 on one half of the genes and -1 on the other half gives +|rho|
## within a half and -|rho| between halves: the genes push each other in
## opposite directions.
## ---------------------------------------------------------------------------
simulate_experiments <- function(N, m0, d, n_arm, rho, neg_pairs, alpha, q,
                                 n_runs, sided = "two", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m1  <- N - m0
  ncp <- d * sqrt(n_arm / 2)
  df  <- 2 * (n_arm - 1)
  mu  <- c(rep(0, m0), rep(ncp, m1))

  ## r is the strength: corr(z_i, z_j) = s_i * s_j * r, so pairs inside a block
  ## correlate +r and pairs across the two blocks correlate -r. As in app 2 the
  ## flips are confined to the TRUE nulls (see the note on the helpers above).
  r     <- abs(rho)
  n_neg <- flips_for_pairs(m0, neg_pairs)
  s     <- rep(1, N)
  if (n_neg > 0) s[sample.int(m0, n_neg)] <- -1

  methods <- names(method_labels)
  per_run <- lapply(seq_len(n_runs), function(i) {
    f <- if (r > 0) rnorm(1) else 0
    z <- mu + s * sqrt(r) * f + sqrt(1 - r) * rnorm(N)
    w <- rchisq(N, df) / df                  # each gene estimates its own variance
    tt <- z / sqrt(w)
    ## A two-sided test reads |t|, which is blind to the sign of a gene's
    ## loading: for a mean-zero gene, flipping the sign leaves the distribution
    ## of |t| unchanged. A one-sided test does see the sign, so it is the only
    ## way negative dependence reaches the p-values at all.
    p <- if (sided == "one") pt(tt, df = df, lower.tail = FALSE)
         else 2 * pt(abs(tt), df = df, lower.tail = FALSE)

    vapply(methods, function(mm) {
      padj <- switch(mm,
                     none   = p,
                     sidak  = p_sidak(p),
                     storey = p_storey(p),
                     p.adjust(p, method = mm))
      ## FWER procedures are judged at alpha, FDR procedures at the FDR level q
      lvl  <- if (mm %in% c("BH", "BY", "storey")) q else alpha
      rej  <- padj <= lvl
      V <- sum(rej[seq_len(m0)])                       # false positives
      S <- if (m1 > 0) sum(rej[(m0 + 1):N]) else 0     # true positives
      R <- V + S
      c(V = V, S = S, R = R,
        FDP = if (R == 0) 0 else V / R,                # Q = V/R, 0 when R = 0
        any_false = as.numeric(V > 0),
        power = if (m1 > 0) S / m1 else NA_real_)
    }, numeric(6))
  })

  arr <- simplify2array(per_run)            # 6 x methods x runs
  summ <- data.frame(
    Method = unname(method_labels[methods]),
    Level  = ifelse(methods %in% c("BH", "BY", "storey"), q, alpha),
    FWER   = apply(arr["any_false", , , drop = FALSE], 2, mean),
    FDR    = apply(arr["FDP", , , drop = FALSE], 2, mean),
    Power  = apply(arr["power", , , drop = FALSE], 2, mean),
    Rejections = apply(arr["R", , , drop = FALSE], 2, mean),
    FalsePositives = apply(arr["V", , , drop = FALSE], 2, mean),
    row.names = NULL, stringsAsFactors = FALSE
  )
  summ$Method <- factor(summ$Method, levels = unname(method_labels))

  fdp <- data.frame(
    Method = factor(rep(unname(method_labels[methods]), each = n_runs),
                    levels = unname(method_labels)),
    FDP    = as.vector(t(arr["FDP", , ]))
  )

  ## one experiment kept for the p-value picture
  f <- if (r > 0) rnorm(1) else 0
  z1 <- mu + s * sqrt(r) * f + sqrt(1 - r) * rnorm(N)
  tt1 <- z1 / sqrt(rchisq(N, df) / df)
  p1 <- if (sided == "one") pt(tt1, df = df, lower.tail = FALSE)
        else 2 * pt(abs(tt1), df = df, lower.tail = FALSE)

  list(summary = summ, fdp = fdp,
       example = data.frame(p = p1,
                            truth = factor(c(rep("TRUE null", m0),
                                             rep("FALSE null", m1)),
                                           levels = c("TRUE null", "FALSE null"))),
       settings = data.frame(N = N, m0 = m0, m1 = m1, pi0 = m0 / N,
                             effect = d, n_arm = n_arm, df = df, sided = sided,
                             dependence = if (r == 0) "independent"
                                          else if (n_neg == 0) "positive" else "negative",
                             neg_pairs = if (r > 0) neg_pair_n(m0, n_neg) else 0,
                             n_neg = n_neg, rho = r,
                             alpha = alpha, q = q, runs = n_runs))
}

## ---------------------------------------------------------------------------
## Guided steps: settings to apply, what to look at, and the hidden answer
## ---------------------------------------------------------------------------
steps <- list(
  list(
    id = "step1", title = "1. All genes null",
    settings = list(N = 1000, pi0 = 100, dependence = "independent", d = 2),
    look = tags$ul(
      tags$li("Every gene is a TRUE null: there is nothing to find."),
      tags$li("Compare the FWER and the FDR columns: they are equal in every row."),
      tags$li("Look at the histogram of the p-values from a single experiment."),
      tags$li("The ", tags$b("Statistical power"), " column is empty: with no real
               effects there is nothing to find.")
    ),
    answer = tags$div(
      tags$p(tags$b("Under the global null the FDR and the FWER are the same thing.")),
      tags$p("With no true effects every rejection is false, so V/R = 1 whenever
              R > 0 and 0 otherwise. The expected proportion of false discoveries
              is then exactly the probability of making at least one: FDR = FWER."),
      tags$p("The histogram of the p-values is flat: when a null hypothesis is TRUE
              its p-value is ", tags$b("uniformly distributed between 0 and 1", .noWS = "outside"),
              ", so any interval of p-values of width 0.05 catches about 5% of the
              genes. That is why, without correction, about 5% of these 1000 genes —
              some 50 — come out significant at 0.05 although none has a real effect."),
      tags$p("This is the one configuration where controlling the FDR is as strict
              as controlling the FWER. Otherwise FDR ≤ FWER, and strictly below it
              as soon as some of the rejections are true discoveries — which is why
              FDR methods can afford to reject more.")
    )
  ),
  list(
    id = "step2", title = "2. A few real effects",
    settings = list(N = 1000, pi0 = 90, dependence = "independent", d = 2),
    look = tags$ul(
      tags$li("10% of the genes are truly differentially expressed."),
      tags$li("Compare the ", tags$b("Statistical power"), " of Bonferroni and Benjamini-Hochberg."),
      tags$li("Check that the FDR of BH sits near π₀ × q, not at q."),
      tags$li("Now look at the ", tags$b("Storey"), " row: its FDR sits closer to q.")
    ),
    answer = tags$div(
      tags$p(tags$b("The two kinds of methods control errors in different ways.")),
      tags$p("Bonferroni tests every gene at α/N: the threshold shrinks with the
              number of genes, so real effects are missed. BH compares the i-th
              smallest p-value with i·q/N: the threshold grows along the ranking,
              so a gene with strong evidence is judged against a far milder bar."),
      tags$p("BH controls the FDR at no more than π₀·q, rather than q. With 90% true
              nulls and high power the realised FDR sits close to 0.9q. That slack is the price of not
              knowing π₀."),
      tags$p("Storey estimates π₀ from the p-values themselves and divides the
              threshold by it, so it spends the whole budget: its FDR sits closer to
              q and it rejects a little more than BH. With 90% true nulls there is
              only 10% to recover — step 6 makes the gain obvious.")
    )
  ),
  list(
    id = "step3", title = "3. More genes",
    settings = list(N = 10000, pi0 = 90, dependence = "independent", d = 2),
    look = tags$ul(
      tags$li("Same proportions as step 2, ten times as many genes."),
      tags$li("Watch the ", tags$b("False positives"), " column for the uncorrected row,
               and BH's statistical power and FDR."),
      tags$li("Then set N back to 1000 and compare.")
    ),
    answer = tags$div(
      tags$p(tags$b("The FWER saturates; the FDR does not move.")),
      tags$p("Uncorrected, the FWER is 1 - (1-α)^m₀. Already with the 900 true nulls
              of step 2 it is 1 to the last decimal, so with 9000 it cannot rise any
              further — what
              grows is the number of false positives, from about 45 to about 450.
              Bonferroni's statistical power collapses as N grows."),
      tags$p("The FDR of BH is a proportion, not a probability. It stays flat at
              about π₀q whatever the number of genes, and the statistical power barely moves.
              That scale-free behaviour is the whole reason the FDR was invented
              for genomics.")
    )
  ),
  list(
    id = "step4", title = "4. Co-regulated genes",
    settings = list(N = 1000, pi0 = 90, dependence = "positive", rho = 0.5, d = 2),
    look = tags$ul(
      tags$li("Genes now share a latent factor: they move up and down together."),
      tags$li("Look at the mean FDR of BH — then look at how far its dots reach to the right in the false discovery proportion plot."),
      tags$li("Look at the FDR of ", tags$b("Storey"), " against q.")
    ),
    answer = tags$div(
      tags$p(tags$b("The average is still controlled; the single experiment is not.")),
      tags$p("The mean FDP stays below q, and the theory says it should: a single
              factor with same-sign loadings is the canonical example of positive
              regression dependence (PRDS), under which Benjamini and Yekutieli
              proved BH valid. For two-sided tests the guarantee needs an extra
              condition on the correlation structure, which this one-factor model
              satisfies; the textbook result also assumes one shared variance
              estimate, whereas here each gene has its own. BH is conservative on
              this screen either way. The point of this step is the spread, not the mean."),
      tags$p("The experiments no longer resemble each other. Four in ten now have no
              false discovery at all, while one in ten has more than twice the share
              you allowed, and a few are almost all wrong — those are the dots far to
              the right. One experiment can be perfectly clean and the next a mess,
              because the shared factor pushes all the genes in the same direction."),
      tags$p("The FDR is an expectation over repeated experiments. It is not a
              guarantee about the experiment on your desk — and you will usually run
              your study only once."),
      tags$p(tags$b("Storey, unlike BH, is not protected here."), " Its FDR climbs
              above q (about 0.08 against q = 0.05). It estimates the proportion of
              TRUE nulls by counting the p-values above 0.5 in the one experiment it
              is given; when every gene moves with the shared factor, that count
              swings from experiment to experiment, and whenever it comes out too low
              the threshold is too lenient. BH estimates nothing, so the same
              dependence leaves it conservative instead.")
    )
  ),
  list(
    id = "step5", title = "5. Genes pushing in opposite directions",
    settings = list(N = 1000, pi0 = 90, dependence = "negative", rho = 0.5,
                    neg_pairs = 202500, d = 2, sided = "one"),
    look = tags$ul(
      tags$li("A share of the genes now loads negatively on the shared factor: when
               one block goes up, the other goes down."),
      tags$li("This step switches the test to ", tags$b("one-sided"), ": a two-sided
               test reads |t| and cannot see the sign of a loading, so the error rates
               would look exactly like step 4."),
      tags$li("Compare the FDR of BH with that of Benjamini-Yekutieli, and compare
               their statistical power.")
    ),
    answer = tags$div(
      tags$p(tags$b("Negative dependence is the one case that BH's proof does not cover.")),
      tags$p("BH's guarantee has been proved for independent tests and for positive regression
              dependence (PRDS). With genes pushing in opposite directions neither
              applies, so the guarantee is no longer automatic."),
      tags$p("BY divides the threshold by c(N) = Σ(1/i) and is valid under
              ", tags$b("any"), " dependence — but c(N) is about 7.5 at N = 1000 and
              about 10 at N = 10,000, so the loss of power is substantial. In
              practice BH is usually kept; the point is to know which assumption
              you are relying on."),
      tags$p(tags$b("Why only one-sided tests?"), " Take two genes, A and B, that
              the treatment does not affect: both are TRUE nulls. They share a hidden
              factor, for example a pathway activity that varies from sample to sample.
              When the factor is high in a sample, A goes up and B goes down. That is
              the negative correlation. In one experiment, three things can happen by
              chance:"),
      tags$table(class = "table table-condensed",
                 style = "width:auto; margin-bottom:8px;",
        tags$thead(tags$tr(tags$th("chance event"), tags$th("t of gene A"),
                           tags$th("t of gene B"), tags$th("two-sided"),
                           tags$th("one-sided (up only)"))),
        tags$tbody(
          tags$tr(tags$td("treated samples happen to have ", tags$b("high"), " factor"),
                  tags$td(tags$b("large +")), tags$td(tags$b("large −")),
                  tags$td(tags$b("both"), " rejected"), tags$td("only A rejected")),
          tags$tr(tags$td("treated samples happen to have ", tags$b("low"), " factor"),
                  tags$td(tags$b("large −")), tags$td(tags$b("large +")),
                  tags$td(tags$b("both"), " rejected"), tags$td("only B rejected")),
          tags$tr(tags$td("the factor is balanced between the groups"),
                  tags$td("≈ 0"), tags$td("≈ 0"), tags$td("neither"), tags$td("neither")))),
      tags$p(tags$b("Two-sided."), " The same chance imbalance in the samples creates a
              false positive in A and in B at the same time; B simply goes in the
              opposite direction. The two-sided test does not care about direction, so
              both false positives fall in the same experiment. That is exactly what
              happens with positive correlation (where B would go up with A), so the
              sign of the correlation makes no difference."),
      tags$p(tags$b("One-sided (looking only for up-regulation)."), " A high-factor
              experiment gives a false positive in A only, and a low-factor experiment
              gives one in B only. The false positives now fall in different
              experiments, so they are spread over more of them — which is the
              structure no theorem covers for BH."),
      tags$p("To summarise: what makes false positives coincide is the shared fluke in
              the samples, not the sign of the correlation. A two-sided test catches the
              fluke in both genes whichever direction each one goes. A one-sided test
              catches it only in the gene that goes up."),
      tags$p(tags$em("One last thing: even with genuine negative dependence BH still
              comes out conservative on this screen. A one-factor model is not a
              counterexample; the point is that no theorem covers this case, not that
              BH fails here."))
    )
  ),
  list(
    id = "step6", title = "6. Many real effects",
    settings = list(N = 1000, pi0 = 50, dependence = "independent", d = 2),
    look = tags$ul(
      tags$li("Half the genes are now truly differentially expressed."),
      tags$li("Compare the ", tags$b("FDR"), " of Benjamini-Hochberg and Storey: BH
               spends half the budget, Storey all of it."),
      tags$li("Then their ", tags$b("Rejections"), " and ", tags$b("Statistical power",
               .noWS = "outside"), ": Storey finds more, though power is already near
               its ceiling here.")
    ),
    answer = tags$div(
      tags$p(tags$b("The fewer the true nulls, the more Storey gains.")),
      tags$p("BH is calibrated for the worst case, as if every gene were a TRUE null.
              Its FDR therefore lands near π₀·q: about 0.5q here, so it spends only
              half the budget it is allowed."),
      tags$p("Storey counts how many p-values fall above 0.5 and uses that to estimate
              the proportion of TRUE nulls. It then divides the threshold by that
              estimate, so its FDR sits close to q and it rejects more genes. With
              90% true nulls (step 2) the gain was small; with 50% it is clear."),
      tags$p("The price is the estimate itself: it needs many tests, and genes that
              all move together make it unreliable. BH, which estimates nothing, is
              the safer choice when the tests are few or strongly correlated.")
    )
  )
)

step_titles <- vapply(steps, function(s) s$title, character(1))

## ---------------------------------------------------------------------------
## UI
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
    .nav-pills > li, .nav-pills > li + li { float:none; margin:0; }
    .nav-pills > li:first-child { flex:0 0 100%%; }
    .nav-pills > li:not(:first-child) { flex:0 0 calc(33.333%% - 4px); }
    .nav-pills > li > a { text-align:left; }
  ", crukdark, crukdark, crukpink, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "FWER and FDR: Bonferroni, Holm, Šidák, BH, BY, Storey")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Each simulated experiment compares <b>two groups</b> of samples across
            <b>N genes</b>, testing every gene for differential expression. Expression is
            <b>normally distributed</b> with the <b>same variance</b> in the two groups, so every
            gene gives a <b>two-sample Student t-test</b> on 2(n&nbsp;&minus;&nbsp;1) degrees of
            freedom. Of those genes,
            <b style='color:%s'>m<sub>0</sub> are true nulls</b> (the two group means are equal) and the rest are
            <b style='color:%s'>truly differentially expressed</b> (their means differ by the effect
            size set on the left). Genes may vary independently or be <b>correlated</b>, in the same
            or in opposite directions. The experiment is repeated many
            times. We can then measure how often a correction lets <b>at least one</b> false positive
            through (<b style='color:%s'>FWER</b>, a probability) and, <b>averaged over experiments</b>,
            what proportion of the genes it calls significant are false
            (<b style='color:%s'>FDR</b>, an expected proportion).", "grey45", crukpink, crukdark, crukpink))),

  ## top menu: the guided steps
  do.call(tabsetPanel, c(
    list(id = "guide", type = "pills", selected = "Free exploration",
         tabPanel("Free exploration",
                  div(class = "step-box",
                      HTML("Set the design on the left and press <b>Run simulation</b>.
                            Or choose one of the numbered steps.")))),
    lapply(steps, function(s) {
      tabPanel(s$title,
               div(class = "step-box",
                   h4(s$title),
                   actionButton(paste0("apply_", s$id), "Apply these settings",
                                class = "btn-sm", style = "margin-bottom:6px;"),
                   s$look))
    })
  )),

  hr(),

  sidebarLayout(
    ## ---- left: the parameters that drive FWER and FDR ----------------------
    sidebarPanel(
      width = 3,
      actionButton("run", "Run simulation", class = "btn-primary",
                   style = sprintf("background:%s;border-color:%s;width:100%%;font-size:15px;", crukdark, crukdark)),
      br(), br(),
      h4("Statistical design"),
      sliderTextInput("N", "Number of genes",
                      choices = c(100, 500, 1000, 5000, 10000, 20000),
                      selected = 1000, grid = TRUE),
      sliderInput("pi0", "TRUE null genes (% of N)",
                  min = 0, max = 100, value = 90, step = 1, post = " %"),
      uiOutput("m0_text"),
      sliderInput("n_arm", "Observations per group", min = 3, max = 30, value = 10, step = 1),
      sliderInput("d", "Effect size: difference in means",
                  min = 0.2, max = 3, value = 2, step = 0.1),
      uiOutput("d_text"),
      radioButtons("sided", "Test",
                   choices = c("Two-sided" = "two",
                               "One-sided (upper tail)" = "one"),
                   selected = "two", inline = TRUE),

      h4("Dependence between genes"),
      radioButtons("dependence", NULL,
                   choices = c("Independent"                    = "independent",
                               "Positive (co-regulated)"        = "positive",
                               "Negative (opposite directions)" = "negative"),
                   selected = "independent"),
      conditionalPanel("input.dependence != 'independent'",
        sliderInput("rho", "Strength of the dependence (|ρ|)",
                    min = 0.05, max = 0.9, value = 0.5, step = 0.05)),
      conditionalPanel("input.dependence == 'negative' && input.sided == 'one'",
        sliderInput("neg_pairs", "TRUE-null gene pairs with negative correlation",
                    min = 0, max = 400, value = 0, step = 1)),
      conditionalPanel("input.dependence == 'negative' && input.sided == 'two'",
        helpText(HTML("Negative correlation can only be set for one-sided tests: with
                  two-sided tests a gene's loading sign changes nothing below, because
                  the tests read |t|. <b>Why?</b> See step 5
                  (<i>Genes pushing in opposite directions</i>), then <i>Show the
                  explanation</i>."))),
      uiOutput("dep_text"),

      h4("Simulation"),
      radioButtons("goal", "Which methods to show?",
                   choices = c("FWER methods (α)" = "fwer",
                               "FDR methods (q)"  = "fdr",
                               "Both kinds, side by side" = "both"),
                   selected = "both"),
      conditionalPanel("input.goal != 'fdr'",
        sliderInput("alpha", "α — the family-wise error rate (Bonferroni, Holm, Šidák)",
                    min = 0.01, max = 0.10, value = 0.05, step = 0.01)),
      conditionalPanel("input.goal != 'fwer'",
        sliderInput("q", "q — the false discovery rate (BH, BY, Storey)",
                    min = 0.01, max = 0.25, value = 0.05, step = 0.01)),
      sliderInput("runs", "Number of simulated experiments",
                  min = 100, max = 10000, value = 5000, step = 100),
      helpText("5000 experiments pin down a rate near 5% to about ±0.006, less
                precisely for rates near a half, and take a few seconds with 1000
                genes, longer with many more genes. Lower the number to explore quickly, or raise it to
                10 000 for sharper numbers."),
      ## a fixed seed makes a run reproducible: the same settings give the same
      ## numbers, so the values quoted in the steps can be checked exactly
      numericInput("seed", "Random seed (empty = a new draw every run)", value = 1),
      downloadButton("save", "Download CSV", style = "width:100%;")
    ),

    ## ---- right: the simulation results -------------------------------------
    mainPanel(
      width = 9,
      tagList(
        fluidRow(
          column(12, h4(textOutput("scenario")), DTOutput("tab"),
                 helpText(HTML("Each gene gives a two-sample t statistic on
                2(n&nbsp;&minus;&nbsp;1) degrees of freedom, drawn from its sampling
                distribution rather than from simulated expression values. Write
                <b>N</b> for the number of genes, of which <b>m<sub>0</sub></b> are TRUE
                nulls (the two group means are equal) and
                <b>N&nbsp;&minus;&nbsp;m<sub>0</sub></b> are FALSE; the proportion of TRUE
                nulls is <b>&pi;<sub>0</sub> = m<sub>0</sub>/N</b>. In one experiment
                <b>V</b> of the TRUE nulls are rejected (false positives), <b>S</b> of
                the FALSE nulls are rejected (true positives), and <b>R</b> = V + S
                genes are reported.<br>
                <b>FWER</b> := P(V &gt; 0), estimated by the share of experiments with
                at least one false positive.
                <b>FDR</b>: of the genes a method reports, the average share that are
                false (counted as 0 when it reports nothing).
                <b>Statistical power</b>: the average of
                S/(N&nbsp;&minus;&nbsp;m<sub>0</sub>), the share of the real effects
                found. <b>Rejections (mean)</b> and <b>False positives (mean)</b>: the
                average <i>counts</i> R and V per experiment, not proportions.<br>
                <b>Bonferroni</b>, <b>Holm</b> and <b>Šidák</b> control the FWER at α. <b>BH</b>,
                <b>BY</b> and <b>Storey</b> control the FDR at q. <b>BH</b> compares the
                i-th smallest p-value with the threshold <b>i&middot;q/N</b> and rejects
                every gene up to the largest i whose p-value is below it. <b>BY</b> uses
                the same rule with the threshold divided by c(N) = &Sigma; 1/i, and
                <b>Storey</b> with the threshold divided by &pi;&#770;<sub>0</sub>, the
                estimated proportion of TRUE nulls.")),
                 br())
        ),
        fluidRow(
          column(6, h4("Observed error rates"), plotOutput("rates", height = "260px")),
          column(6, h4("False discovery proportion, experiment by experiment"),
                 plotOutput("fdp", height = "260px"))
        ),
        br(),
        fluidRow(
          column(6, h4("p-values from a single experiment"), plotOutput("pv", height = "230px")),
          column(6,
                 h4("What is going on"),
                 uiOutput("reveal_ui"),
                 uiOutput("answer"))
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  m0_now <- reactive(round(input$N * input$pi0 / 100))

  ## a value set by a step's button is kept here until the range has caught up:
  ## sent on its own into the old range it would be snapped to that range's step
  pending_neg <- reactiveVal(NULL)
  observeEvent(m0_now(), {
    mx  <- max_neg_pairs(m0_now())
    cur <- if (!is.null(pending_neg())) pending_neg()
           else if (is.null(input$neg_pairs)) 0 else input$neg_pairs
    pending_neg(NULL)
    updateSliderInput(session, "neg_pairs", max = max(mx, 1),
                      value = neg_pair_n(m0_now(), flips_for_pairs(m0_now(), min(cur, mx))),
                      step = 1)
  })

  ## only counts k(m0 - k) can be built, so the slider jumps to the nearest one
  observeEvent(input$neg_pairs, {
    req(is.null(pending_neg()), m0_now() >= 2)
    np <- neg_pair_n(m0_now(), flips_for_pairs(m0_now(), input$neg_pairs))
    if (np != input$neg_pairs) updateSliderInput(session, "neg_pairs", value = np)
  }, ignoreInit = TRUE)

  output$dep_text <- renderUI({
    if (input$dependence == "independent")
      return(helpText("All genes vary independently."))
    if (input$dependence == "positive")
      return(helpText(sprintf("The shared factor correlates every pair of test statistics at about +%.2f.", input$rho)))
    ## with a two-sided test the slider is hidden and the panel above explains why
    if (input$sided == "two") return(NULL)
    req(m0_now() >= 2)
    m0 <- m0_now()
    k  <- flips_for_pairs(m0, input$neg_pairs)
    np <- neg_pair_n(m0, k)
    ## pairs are unordered, and a pair is negative when it has one gene from each
    ## block, so the count is k(m0 - k); say how it is reached
    pairs_txt <- sprintf("Pairs are unordered (gene A with gene B is the same pair as gene B
                          with gene A), so the %d TRUE-null genes form %d × %d / 2 = %d %s.",
                         m0, m0, m0 - 1, choose(m0, 2), if (m0 == 2) "pair" else "pairs")
    ## from the number of negative pairs to the number of genes k that move
    ## against the shared factor: solve k(m0 - k) = pairs
    reach <- unique(neg_pair_n(m0, 0:floor(m0 / 2)))
    reach_txt <- paste0(paste(head(reach, 5), collapse = ", "),
                        if (length(reach) > 5) ", …" else "")
    body <- if (np == 0)
      "No pair is negative: every gene moves with the shared factor."
    else sprintf(
      "If k genes move against the shared factor and the other %d − k with it, a pair
       is negative exactly when it has one gene from each group, so there are
       k × (%d − k) negative pairs. For %d negative %s, solving k × (%d − k) = %d gives
       k = (%d − √(%d² − 4 × %d)) / 2 = %d: <b>%d %s against the shared factor</b>
       (%.1f%% of all pairs are negative)%s.",
      m0, m0, np, if (np == 1) "pair" else "pairs", m0, np, m0, m0, np, k,
      k, if (k == 1) "gene moves" else "genes move", max(neg_pair_pct(m0, k), 0.1),
      if (np == 1) ", and this pair involves both TRUE-null genes"
      else sprintf(", and these pairs together involve all %d TRUE-null genes", m0))
    ## the two blocks: +|rho| within each, -|rho| between them
    if (np > 0 && choose(k, 2) + choose(m0 - k, 2) > 0) body <- paste(body, sprintf(
      "The TRUE-null genes thus form two blocks, of %d and %d genes: <b>positively
       correlated within each block</b> (%d + %d = %d pairs) and <b>negatively correlated
       between them</b> (%d %s).",
      k, m0 - k, choose(k, 2), choose(m0 - k, 2), choose(k, 2) + choose(m0 - k, 2),
      np, if (np == 1) "pair" else "pairs"))
    body <- paste(body, sprintf(
      "Only counts of this form can be built (%s), so the slider jumps to the nearest one.",
      reach_txt))
    helpText(HTML(paste(pairs_txt, body,
      "With one-sided tests, only negatively correlated pairs can break BH's and Šidák's guarantee; Bonferroni, Holm and BY hold regardless. Storey is a separate matter: it can fail under positive dependence too (step 4).")))
  })

  output$d_text <- renderUI({
    helpText(sprintf("The difference between the mean of a differentially expressed
                      gene and the mean of a TRUE null gene: at %s, their means differ
                      by %s.", format(input$d), format(input$d)))
  })

  output$m0_text <- renderUI({
    helpText(sprintf("m₀ = %d TRUE nulls (%.0f%% of N = %d); m₁ = %d truly DE genes",
                     m0_now(), input$pi0, input$N, input$N - m0_now()))
  })

  ## guided steps write their settings into the sidebar
  lapply(steps, function(s) {
    observeEvent(input[[paste0("apply_", s$id)]], {
      st <- s$settings
      updateSliderTextInput(session, "N", selected = st$N)
      updateSliderInput(session, "pi0", value = st$pi0)
      updateRadioButtons(session, "dependence", selected = st$dependence)
      if (!is.null(st$rho))       updateSliderInput(session, "rho", value = st$rho)
      ## the range, the step and the value go together, so the value cannot be
      ## snapped to the previous range's step
      if (!is.null(st$neg_pairs)) {
        m0_st  <- round(st$N * st$pi0 / 100)
        neg_st <- neg_pair_n(m0_st, flips_for_pairs(m0_st, st$neg_pairs))
        updateSliderInput(session, "neg_pairs", max = max(max_neg_pairs(m0_st), 1),
                          step = 1, value = neg_st)
        if (m0_st != m0_now()) pending_neg(neg_st)
      }
      if (!is.null(st$d))   updateSliderInput(session, "d", value = st$d)
      ## the sidedness of the test is part of a step's setup: step 5 needs a
      ## one-sided test for the negative structure to reach the p-values at all,
      ## and the other steps must be put back to two-sided
      updateRadioButtons(session, "sided",
                         selected = if (is.null(st$sided)) "two" else st$sided)
      ## the numbers quoted in the answers assume the rest of the scenario too
      updateSliderInput(session, "n_arm", value = 10)
      updateSliderInput(session, "alpha", value = 0.05)
      updateSliderInput(session, "q", value = 0.05)
      updateRadioButtons(session, "goal", selected = "both")
      showNotification("Settings applied - press Run simulation", type = "message",
                       duration = 3)
    }, ignoreInit = TRUE)
  })

  sim <- eventReactive(input$run, ignoreNULL = FALSE, {
    withProgress(message = "Simulating experiments...", value = 0.4, {
      simulate_experiments(
        N = input$N, m0 = m0_now(), d = input$d, n_arm = input$n_arm,
        rho = if (input$dependence == "independent") 0 else input$rho,
        neg_pairs = if (input$dependence == "negative" && input$sided == "one" &&
                        !is.null(input$neg_pairs)) input$neg_pairs else 0,
        alpha = input$alpha, q = input$q, n_runs = input$runs,
        sided = input$sided,
        seed = if (is.null(input$seed) || is.na(input$seed)) NULL else input$seed
      )
    })
  })

  ## which procedures belong to the chosen goal
  shown <- reactive({
    switch(input$goal,
      fwer = c("No correction", "Bonferroni", "Holm", "Šidák"),
      fdr  = c("No correction", "Benjamini-Hochberg", "Benjamini-Yekutieli", "Storey"),
      unname(method_labels))
  })

  output$scenario <- renderText({
    s <- sim()$settings
    sprintf("%d genes, %d TRUE nulls (%.0f%%), %d DE genes | %s test | %s dependence%s | %d experiments | α = %.2f, q = %.2f",
            s$N, s$m0, 100 * s$pi0, s$m1, if (s$sided == "one") "one-sided" else "two-sided",
            if (s$rho == 0) "no" else if (s$neg_pairs == 0) "positive" else "negative",
            if (s$rho == 0) "" else if (s$neg_pairs == 0)
              sprintf(" (|ρ| = %.2f%s)", s$rho,
                      if (input$dependence == "negative")
                        "; a two-sided test cannot see the sign, so no pair is negative" else "")
              else sprintf(" (|ρ| = %.2f, %s negative TRUE-null pairs)", s$rho,
                           format(s$neg_pairs, big.mark = ",", scientific = FALSE)),
            s$runs, s$alpha, s$q)
  })

  output$tab <- renderDT({
    d <- sim()$summary
    d <- d[d$Method %in% shown(), , drop = FALSE]
    datatable(d, rownames = FALSE, options = list(dom = "t", ordering = FALSE),
              colnames = c("Method", "Level used", "FWER", "FDR", "Statistical power",
                           "Rejections (mean)", "False positives (mean)")) |>
      formatRound(c("FWER", "Power"), 3) |>
      formatRound("Level", 2) |>
      formatRound("FDR", 3) |>
      formatRound("Rejections", 1) |>
      formatRound("FalsePositives", 2) |>
      formatStyle("FDR", fontWeight = "bold", color = crukpink) |>
      formatStyle("FWER", color = crukdark)
  })

  output$rates <- renderPlot({
    d <- sim()$summary
    d <- d[d$Method %in% shown(), , drop = FALSE]
    dd <- rbind(data.frame(Method = d$Method, Rate = "FWER", Value = d$FWER),
                data.frame(Method = d$Method, Rate = "FDR",  Value = d$FDR))
    dd$Rate <- factor(dd$Rate, levels = c("FWER", "FDR"))
    ggplot(dd, aes(Method, Value, fill = Rate)) +
      geom_col(position = position_dodge(0.75), width = 0.7) +
      {if (input$goal != "fdr")
         geom_hline(yintercept = sim()$settings$alpha, linetype = 2, colour = crukdark)} +
      {if (input$goal != "fwer")
         geom_hline(yintercept = sim()$settings$q, linetype = 3, colour = crukpink)} +
      scale_fill_manual(values = c(FWER = crukdark, FDR = crukpink)) +
      labs(x = NULL, y = NULL,
           caption = switch(input$goal,
             fwer = "dashed: α",
             fdr  = "dotted: q",
             "dashed: α (FWER methods)    dotted: q (FDR methods)")) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            legend.position = "top", legend.title = element_blank())
  })

  ## one row per method, so nothing overlaps: the box is the spread of the FDP
  ## over the experiments, the diamond its mean (which is the FDR)
  output$fdp <- renderPlot({
    d <- sim()$fdp
    d <- d[d$Method %in% shown(), ]
    d$Method <- factor(d$Method, levels = rev(intersect(unname(method_labels), shown())))
    m <- aggregate(FDP ~ Method, d, mean)
    ggplot(d, aes(FDP, Method, colour = Method, fill = Method)) +
      {if (input$goal != "fwer")
         geom_vline(xintercept = sim()$settings$q, linetype = 2, colour = crukpink)} +
      geom_boxplot(alpha = 0.25, outlier.size = 0.6, outlier.alpha = 0.4, width = 0.6) +
      geom_point(data = m, shape = 23, size = 3.2, colour = "black", stroke = 0.6) +
      scale_colour_manual(values = method_cols) +
      scale_fill_manual(values = method_cols) +
      ## square root, not log: the FDP is exactly 0 in most experiments for the
      ## FWER methods, and a log axis would drop those points altogether.
      ## coord_trans, not scale_x_sqrt: the scale would transform the data before
      ## the box statistics, so the whiskers would be drawn on the square roots
      ## coord_transform in ggplot2 >= 4.0, coord_trans before that
      (if (exists("coord_transform")) coord_transform(x = "sqrt") else coord_trans(x = "sqrt")) +
      scale_x_continuous(breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) +
      labs(x = "False discovery proportion (V / R) in one experiment", y = NULL,
           caption = paste("square-root scale   box: where most experiments fall   diamond: mean = FDR",
                           if (input$goal != "fwer") "  dashed: q" else "")) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
  })

  output$pv <- renderPlot({
    d <- sim()$example
    ggplot(d, aes(p, fill = truth)) +
      geom_histogram(bins = 40, position = "stack") +
      scale_fill_manual(values = c("TRUE null" = "grey70", "FALSE null" = crukpink)) +
      labs(x = "p-value", y = "genes") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top", legend.title = element_blank())
  })

  ## the phenomenon stays hidden until the button is pressed
  current_step <- reactive({
    ## length 0 before the tab menu has reported its value
    i <- match(input$guide, step_titles)
    if (length(i) == 0 || is.na(i)) NULL else steps[[i]]
  })

  output$reveal_ui <- renderUI({
    st <- current_step()
    if (is.null(st)) {
      helpText("Pick one of the numbered steps in the menu at the top: each one
                hides the phenomenon until you have run the simulation and
                looked at the numbers yourself.")
    } else {
      actionButton("reveal", "Show the explanation", class = "btn-sm",
                   style = sprintf("background:%s;border-color:%s;color:#FFF;", crukpink, crukpink))
    }
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
    filename = function() sprintf("fdr-simulation-%s.csv", format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      res <- sim()
      out <- res$summary[res$summary$Method %in% shown(), , drop = FALSE]
      write.csv(cbind(res$settings[rep(1, nrow(out)), ], out),
                file, row.names = FALSE)
    }

  )
}

shinyApp(ui, server)
