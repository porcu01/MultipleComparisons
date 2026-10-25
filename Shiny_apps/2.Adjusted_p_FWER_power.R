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
  holm_sidak  = "Holm–Šidák"
)
method_cols <- c("No correction" = "grey45", "Bonferroni" = crukdark,
                 "Holm" = cruklight, "Šidák" = "#5B8C00",
                 "Holm–Šidák" = "#7B2D8E")
## the four corrections often land within a few thousandths of each other, so a
## colour alone leaves points hidden underneath one another
method_shapes <- c("No correction" = 4, "Bonferroni" = 16, "Holm" = 17,
                   "Šidák" = 15, "Holm–Šidák" = 18)

## the run counts are shown as 100,000 rather than 1e+05
runs_label <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
runs_value <- function(x) as.numeric(gsub(",", "", x))

adjust_p <- function(p, method) {
  n <- length(p)
  if (method == "none") p
  ## -expm1(k*log1p(-p)) instead of 1-(1-p)^k: the direct form rounds 1-p to 1
  ## and returns an adjusted p of exactly 0 for p below about 1e-16
  else if (method == "sidak") pmin(-expm1(n * log1p(-p)), 1)
  else if (method == "holm_sidak") {
    o <- order(p); ro <- order(o)
    pmin(cummax(-expm1((n - seq_len(n) + 1) * log1p(-p[o])))[ro], 1)
  } else p.adjust(p, method = method)
}

## ---------------------------------------------------------------------------
## One experiment = expression measured for n genes on the SAME 2 * n_arm
## samples, n_arm per group, then a two-sample t-test per gene.
## The population is jointly normal: for gene i in sample j,
##     x_ij = d * 1{j treated, gene i FALSE null}
##            + s_i * sqrt(|rho|) * F_j + sqrt(1 - |rho|) * e_ij
## with one latent factor F_j per SAMPLE (a pathway, a batch), so every gene has
## unit variance and corr(gene i, gene j) = s_i s_j |rho|, and d is in population
## SD units. Because the genes share the samples, their pooled variance
## estimates are correlated exactly as they would be in a real experiment, so
## the t statistics inherit |rho| almost exactly from about n_arm = 10 upwards
## (some attenuation remains at the smallest sample sizes, where the pooled
## variances are too noisy to track each other).
## ---------------------------------------------------------------------------
## a pair of genes is negatively correlated when exactly one of the two loads
## the opposite way, so flipping k of m0 genes makes k(m0 - k) of the
## m0(m0 - 1)/2 pairs negative. Only the TRUE nulls matter here: the condition
## Sidak and Holm-Sidak need (PUOD) constrains the true-null p-values alone, so
## negative dependence among the genes that really are differentially expressed
## cannot cause a false positive. The flips are therefore confined to the true
## nulls. k(m0 - k) grows in steps, so not every count is reachable.
flips_for_pairs <- function(m0, want) {
  k <- 0:floor(m0 / 2)
  k[which.min(abs(k * (m0 - k) - want))]
}
neg_pair_n   <- function(m0, k) k * (m0 - k)
neg_pair_pct <- function(m0, k) if (m0 < 2) 0 else 100 * k * (m0 - k) / choose(m0, 2)
max_neg_pairs <- function(m0) floor(m0 / 2) * ceiling(m0 / 2)

simulate_fwer <- function(n, m0, d, n_arm, rho, neg_pairs, alpha, n_runs,
                          sided = "two", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m1  <- n - m0
  df  <- 2 * (n_arm - 1)
  N   <- 2 * n_arm
  grp <- rep(c(0, 1), each = n_arm)
  ## the treatment shifts the FALSE nulls in the second group by d SDs
  shift <- outer(c(rep(0, m0), rep(d, m1)), grp)

  ## cap below 1: at |rho| = 1 every gene is an exact +/- copy of the factor
  r     <- min(abs(rho), 0.999)
  n_neg <- flips_for_pairs(m0, neg_pairs)
  s     <- rep(1, n)
  if (n_neg > 0) s[sample.int(m0, n_neg)] <- -1
  load  <- s * sqrt(r)

  methods <- names(method_labels)
  per_run <- lapply(seq_len(n_runs), function(i) {
    f  <- rnorm(N)                      # one factor value per sample
    X  <- outer(load, f) + sqrt(1 - r) * matrix(rnorm(n * N), n, N) + shift
    mA <- rowMeans(X[, grp == 0, drop = FALSE])
    mB <- rowMeans(X[, grp == 1, drop = FALSE])
    sp <- sqrt((rowSums((X[, grp == 0, drop = FALSE] - mA)^2) +
                rowSums((X[, grp == 1, drop = FALSE] - mB)^2)) / df)
    tt <- (mB - mA) / (sp * sqrt(2 / n_arm))
    p  <- if (sided == "one") pt(tt, df = df, lower.tail = FALSE)
          else 2 * pt(abs(tt), df = df, lower.tail = FALSE)

    vapply(methods, function(mm) {
      rej <- adjust_p(p, mm) <= alpha
      V <- sum(rej[seq_len(m0)])
      S <- if (m1 > 0) sum(rej[(m0 + 1):n]) else 0
      c(V = V, S = S,
        any_false = as.numeric(V > 0),
        ## false discovery proportion of this experiment: V / (V + S), counted
        ## as 0 when nothing was rejected. Its average over runs is the FDR.
        fdp       = if (V + S > 0) V / (V + S) else 0,
        power     = if (m1 > 0) S / m1 else NA_real_,
        any_power = if (m1 > 0) as.numeric(S > 0) else NA_real_)
    }, numeric(6))
  })

  arr <- simplify2array(per_run)
  summ <- data.frame(
    Method = unname(method_labels[methods]),
    FWER   = apply(arr["any_false", , , drop = FALSE], 2, mean),
    Power  = apply(arr["power", , , drop = FALSE], 2, mean),
    AnyPower = apply(arr["any_power", , , drop = FALSE], 2, mean),
    FalsePositives = apply(arr["V", , , drop = FALSE], 2, mean),
    TruePositives  = apply(arr["S", , , drop = FALSE], 2, mean),
    FDR    = apply(arr["fdp", , , drop = FALSE], 2, mean),
    row.names = NULL, stringsAsFactors = FALSE
  )
  summ$Method <- factor(summ$Method, levels = unname(method_labels))
  list(summary = summ,
       settings = data.frame(n = n, m0 = m0, m1 = m1, effect = d, n_arm = n_arm,
                             df = df, rho = r,
                             neg_pairs = if (r > 0) neg_pair_n(m0, n_neg) else 0,
                             sided = sided, alpha = alpha, runs = n_runs))
}

## ---------------------------------------------------------------------------
steps <- list(
  list(id = "s1", title = "1. The error rate grows with the family of hypotheses",
       settings = list(n = 20, pi0 = 100, rho = 0, d = 2, runs = 25000),
       look = tags$ul(
         tags$li("Every hypothesis is a TRUE null, so any rejection is a mistake."),
         tags$li("Read the ", tags$b("FWER"), " of ", tags$b("No correction"),
                 " and compare it with α."),
         tags$li("Then raise n from 20 to 100 and run again.")),
       answer = tags$div(
         tags$p(tags$b("Testing at α per hypothesis does not hold the family of hypotheses at α.")),
         tags$p("With m₀ independent true nulls the chance of at least one false
                positive is 1 − (1 − α)^m₀: about 0.64 at n = 20 and 0.99 at
                n = 100, for α = 0.05. All four corrections stay at or below α instead,
                which is what controlling the FWER means. Šidák and
                Holm–Šidák sit exactly at α when the tests are independent, so the
                simulation can show them a hair above 0.050: with 25 000 experiments
                anything within about ±0.003 is simulation noise."))),

  list(id = "s2", title = "2. Holm beats Bonferroni, for free",
       ## 25 000 runs: at 1000 the Monte Carlo error on the FWER is about 0.007,
       ## wide enough that Holm can read above alpha by chance
       settings = list(n = 50, pi0 = 80, rho = 0, d = 2, runs = 25000),
       look = tags$ul(
         tags$li("Compare the ", tags$b("FWER"), " of Bonferroni and Holm — both at or
                  below α, up to simulation noise."),
         tags$li("Now compare their ", tags$b("Statistical power", .noWS = "outside"), ".")),
       answer = tags$div(
         tags$p(tags$b("Same guarantee, more rejections.")),
         tags$p("Bonferroni tests every hypothesis at α/n. Holm uses the same
                threshold for the smallest p-value, then relaxes it to α/(n−1),
                α/(n−2) and so on as hypotheses are rejected."),
         tags$p("So Holm rejects everything Bonferroni rejects and sometimes more:
                it is uniformly more powerful, at no cost in error control and
                under no extra assumption. There is never a reason to prefer Bonferroni
                on grounds of ", tags$b("statistical power", .noWS = "outside"), "."))),

  list(id = "s3", title = "3. What Šidák adds, and what it assumes",
       ## 50 000 runs so that a third-decimal difference in power is real
       settings = list(n = 50, pi0 = 80, rho = 0, d = 2, runs = 50000),
       look = tags$ul(
         tags$li("Compare ", tags$b("Šidák"), " with ", tags$b("Bonferroni", .noWS = "outside"),
                 ", and ", tags$b("Holm–Šidák"), " with ", tags$b("Holm", .noWS = "outside"), "."),
         tags$li("The differences are small — look at the third decimal of ", tags$b("Statistical power", .noWS = "outside"), ".")),
       answer = tags$div(
         tags$p(tags$b("A small gain, bought with an assumption.")),
         tags$p("Šidák tests at 1 − (1−α)^(1/n) instead of α/n. That is
                slightly larger — 0.001025 against 0.001000 at n = 50, α = 0.05 —
                so Šidák is uniformly a little more powerful than Bonferroni, and
                Holm–Šidák a little more powerful than Holm. With 50 000 experiments
                the gain is visible in the third decimal: Bonferroni 0.695 against
                Šidák 0.698, Holm 0.713 against Holm–Šidák 0.717."),
         tags$p("The price is an assumption: Šidák needs the TRUE-null p-values to be at
                least as likely to be all non-significant together as they would be if
                they were independent (a sufficient condition), whereas Bonferroni and
                Holm need nothing at all.
                Independence is the textbook case, but it is not the condition itself —
                step 4 breaks independence and Šidák survives, while step 5 shows the
                kind of dependence that really does break it."))),

  list(id = "s4", title = "4. When the tests are not independent",
       ## same run count as step 3, so the two tables can be compared directly
       settings = list(n = 50, pi0 = 80, rho = 0.6, d = 2, runs = 50000),
       look = tags$ul(
         tags$li("The genes now share one source of variation, so their tests
                  move together."),
         tags$li("Watch the ", tags$b("FWER"), " of all four corrections: it drops
                  well below α."),
         tags$li("Now compare the ", tags$b("Statistical power", .noWS = "outside"), " with step 3, where the genes
                  were independent — and notice how little it has moved.")),
       answer = tags$div(
         tags$p(tags$b("Positive dependence makes every correction conservative.")),
         tags$p("Bonferroni and Holm rest on Boole's inequality, which holds whatever
                the dependence, so they remain valid — just further below α than
                before. Šidák's independence assumption is violated, but for
                positive dependence like this it stays valid too, and also becomes
                conservative."),
         tags$p("Notice what barely changes: in expectation Bonferroni's and Šidák's
                ", tags$b("statistical power", .noWS = "outside"), " is exactly what it was in step 3. Their threshold is fixed at α/n
                (or its Šidák twin), and the shared factor leaves each gene's own test
                exactly as it was — only their joint behaviour changes. Holm and
                Holm–Šidák lose a thousandth or two in expectation (Holm about 0.714 →
                0.711 over a very large number of experiments), because their step-down
                relaxation depends on how many hypotheses get rejected, and dependence
                makes rejections arrive in clumps. That loss is below what 50 000
                experiments can resolve, so here the last decimal may not move at all."),
         tags$p("What is lost is the ", tags$b("statistical power", .noWS = "outside"), " you could have had. These methods now
                spend only about half of their 0.05 error budget — partly because only
                40 of the 50 genes are TRUE nulls (in step 3 they already sat near 0.04),
                and mostly because of the dependence — and they cannot reclaim the rest, because they are built to be safe whatever the
                dependence: a method that knew how strongly the genes move together
                could afford a looser threshold and turn that unused budget into real
                ", tags$b("statistical power", .noWS = "outside"), ". The clumping does show
                up in one column: ", tags$b("Any ", tags$b("statistical power", .noWS = "outside"), ""), " falls a little, because
                results now tend to arrive together or not at all."))),

  list(id = "s5", title = "5. When FWER control actually fails",
       settings = list(n = 2, pi0 = 100, rho = 0.9, d = 2, neg_pairs = 1,
                       sided = "one", alpha = 0.10, runs = 200000),
       look = tags$ul(
         tags$li("Two genes, both TRUE nulls, pushed hard in ",
                 tags$b("opposite directions"), " (|ρ| = 0.9, the one pair negative)."),
         tags$li("One-sided tests at α = 0.10, and 200 000 experiments — the breach
                  is tiny, so it needs both to be seen at all (about 15–20 seconds)."),
         tags$li("Look at the ", tags$b("Šidák"), " FWER against α = 0.10.")),
       answer = tags$div(
         tags$p(tags$b("Šidák and Holm–Šidák rise above α: about 0.102 against a
                       promised 0.100, while Bonferroni and Holm sit right on α = 0.100.")),
         tags$p("Šidák needs the TRUE-null p-values to be at least as likely to be all
                non-significant together as under independence. Pushing the two genes
                in opposite directions reverses that inequality, and the guarantee goes
                with it. Bonferroni and Holm are unmoved: they assume nothing about
                dependence."),
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
                experiments, so more experiments contain at least one. The FWER goes up,
                and that is what breaks Šidák."),
         tags$p("To summarise: what makes false positives coincide is the shared fluke in
                the samples, not the sign of the correlation. A two-sided test catches the
                fluke in both genes whichever direction each one goes. A one-sided test
                catches it only in the gene that goes up."),
         tags$p("In numbers, for two genes, both TRUE nulls, Šidák at α = 0.10, normal
                statistics:"),
         tags$table(class = "table table-condensed",
                    style = "width:auto; margin-bottom:8px;",
           tags$thead(tags$tr(tags$th("correlation"), tags$th("one-sided FWER"),
                              tags$th("two-sided FWER"))),
           tags$tbody(
             tags$tr(tags$td("0"),     tags$td("0.1000"), tags$td("0.1000")),
             tags$tr(tags$td("+0.9"),  tags$td("0.0698"), tags$td("0.0721")),
             tags$tr(tags$td(tags$b("−0.9")), tags$td(tags$b("0.1026 (above α)")),
                     tags$td(tags$b("0.0721 (same as +0.9)"))))),
         tags$p("This is Šidák's inequality: for normal statistics with any
                correlations, the |t| values are at least as likely to be all small
                together as if they were independent (it also holds for t statistics
                sharing one variance estimate). So two-sided Šidák controls the FWER
                whatever the correlation, without needing independence. It is a
                property of normal and t statistics: for other statistics a two-sided
                test does not guarantee it."))),

  list(id = "s6", title = "6. Many tests, few real effects",
       ## 25 000 runs so that the quoted numbers are stable (a few seconds)
       settings = list(n = 200, pi0 = 95, rho = 0, d = 2, n_arm = 5, runs = 25000),
       look = tags$ul(
         tags$li("200 genes, of which 10 really are differentially expressed, and only
                  five samples per group — the setting used in the lecture."),
         tags$li("Look at ", tags$b("Statistical power", .noWS = "outside"), " — the share of the real effects found."),
         tags$li("Then look at ", tags$b("Any ", tags$b("statistical power", .noWS = "outside"), ""), " — the chance of finding even
                  one of them.")),
       answer = tags$div(
         tags$p(tags$b("FWER control is expensive when the family of hypotheses is large.")),
         tags$p("Guarding against even one false positive across 200 tests forces a
                threshold near α/200. Even an effect as large as two standard
                deviations is then found only about 5% of the time: the ", tags$b("statistical power", .noWS = "outside"), " column
                falls from roughly 0.79 without correction to about 0.05."),
         tags$p("And the ", tags$b("Any statistical power"), " of the four corrections is
                about 0.42 — with ten real effects in the data, more than half of the
                experiments find none of them."),
         tags$p("Meanwhile their ", tags$b("FDR"), " sits just under 0.04. That is an average
                over all experiments, with those in which a method reports nothing counted
                as 0 — and here more than half report nothing. In the experiments where
                they do report something, about one reported gene in twelve is wrong. So
                false discoveries are rare, and that is exactly what you pay for — the
                chance of detecting a gene that really is differentially expressed is
                low too."),
         tags$p("This is the Achilles' heel of the FWER, and it is why the false
                discovery rate was invented: if you can accept that a controlled
                proportion of your findings are wrong, the thresholds can relax — a
                little here, a great deal when the real effects are many. That is the
                subject of app 4.")))
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
    /* Bootstrap's clearfix ::before/::after would become flex items and eat a gap */
    .nav-pills::before, .nav-pills::after { display:none; }
    .nav-pills > li, .nav-pills > li + li { float:none; margin:0; }
    .nav-pills > li:first-child { flex:0 0 100%%; }
    .nav-pills > li:not(:first-child) { flex:0 0 calc(33.333%% - 4px); }
    .nav-pills > li > a { text-align:left; }
  ", crukdark, crukdark, crukpink, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "FWER and statistical power: Bonferroni, Holm, Šidák, Holm–Šidák")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Think of a <b>gene expression</b> study: <b>n genes</b> measured in
            <b>treated</b> and <b>control</b> samples, one hypothesis per gene.
            <b>In the population</b> the genes are <b>jointly normally
            distributed</b>, with the <b>same variance</b> in the two groups and a mean
            that differs between them only if the gene really responds to the treatment.
            Each simulated experiment draws <b>one set of samples</b> — the <b>same</b>
            individuals measured for <b>every</b> gene, in two <b>independent</b> groups
            of equal size (set in <i>Samples per group</i>) — and then runs a
            <b>two-sample Student t-test</b> (one- or two-sided) on each gene, comparing
            the difference between the two sample means with the variability within the
            groups. Because the genes share the samples, both the measurements <b>and</b>
            their variance estimates are correlated, as in a real experiment. Of the n genes,
            <b style='color:%s'>m<sub>0</sub> are not differentially expressed</b>
            (true nulls) and the rest, <b style='color:%s'>n − m<sub>0</sub></b>, really
            are (they all have the same true difference between groups, set in
            <i>Effect size of the false nulls</i>).
            Genes are <b>co-regulated</b> through one shared source of variation — a
            pathway, a batch — so for every pair of genes the correlation is <b>0</b> when
            <i>Strength |ρ|</i> is 0 (genes vary independently), <b>+|ρ|</b> when the two
            move in the same direction with it, and <b>−|ρ|</b> when the shared factor
            pushes one gene up while it pushes the other down. The slider <i>TRUE-null gene
            pairs with negative correlation</i> sets how many such pairs there are among the
            TRUE nulls.
            Repeating the
            experiment many times (set in <i>Number of simulated experiments</i>) shows
            how often each correction calls <b>at least one</b> unchanged gene
            differentially expressed (<b style='color:%s'>FWER</b>) and what share of the
            real effects it finds (<b>statistical power</b>).", crukdark, crukpink, crukdark))),

  fluidRow(
    column(6, h4("One gene: the test statistic"),
           plotOutput("marginal", height = "250px")),
    column(6, h4("Two TRUE-null genes: the joint distribution"),
           plotOutput("bivariate", height = "250px"))),
  br(),

  do.call(tabsetPanel, c(
    list(id = "guide", type = "pills",
         tabPanel("Free exploration",
                  div(class = "step-box",
                      HTML("Set the parameters on the left and press <b>Run simulation</b>.
                            Or choose one of the numbered steps above.")))),
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
      h4("The family of hypotheses"),
      sliderTextInput("n", "n: number of hypotheses",
                      choices = c(2, 5, 10, 20, 50, 100, 200, 500), selected = 50, grid = TRUE),
      sliderInput("pi0", "TRUE nulls (% of n)", min = 0, max = 100, value = 80, step = 5, post = " %"),
      uiOutput("m0_text"),

      h4("Experiment"),
      sliderInput("d", "Effect size of the false nulls (SD units)",
                  min = 0.2, max = 3, value = 2, step = 0.1),
      sliderInput("n_arm", "Samples per group", min = 3, max = 30, value = 10, step = 1),
      radioButtons("sided", "Test", choices = c("Two-sided" = "two",
                                                "One-sided (upper tail)" = "one"),
                   selected = "two"),
      sliderInput("alpha", "α — the FWER we want to hold",
                  min = 0.01, max = 0.10, value = 0.05, step = 0.01),
      ## a FWER that sits just above alpha needs many runs to be seen at all:
      ## at 2000 runs the Monte Carlo error is about 0.005, as wide as the effect
      sliderTextInput("runs", "Number of simulated experiments",
                      choices = runs_label(c(50, 100, 250, 500, 1000, 2000, 5000,
                                             10000, 25000, 50000, 100000, 200000)),
                      selected = runs_label(10000), grid = TRUE),
      ## a fixed seed makes a run reproducible: the same settings give the same
      ## numbers, so the values quoted in the steps can be checked exactly
      numericInput("seed", "Random seed (empty = a new draw every run)", value = 1),

      h4("Dependence between tests"),
      sliderInput("rho", "Strength |ρ|  (0 = independent)",
                  min = 0, max = 0.9, value = 0, step = 0.05),
      ## only offered for one-sided tests: with two-sided tests the sign of a
      ## gene's loading cancels out of every quantity below, because |t| has an
      ## even density, so the slider would be a no-op
      conditionalPanel("input.rho > 0 && input.sided == 'one'",
        sliderInput("neg_pairs", "TRUE-null gene pairs with negative correlation",
                    min = 0, max = 400, value = 0, step = 1),
        uiOutput("neg_text")),
      conditionalPanel("input.rho > 0 && input.sided == 'two'",
        helpText(HTML("Negative correlation can only be set for one-sided tests: with
                  two-sided tests Šidák holds at any correlation in this model, and
                  flipping a gene changes nothing below. <b>Why?</b> See step 5
                  (<i>When FWER control actually fails</i>), then <i>Show the
                  explanation</i>."))),

      downloadButton("save", "Download CSV", style = "width:100%;")
    ),

    mainPanel(
      width = 9,
      h4(textOutput("scenario")),
      DTOutput("tab"),
      helpText(HTML("Each gene gives a two-sample t statistic on 2 × (samples per group − 1)
                degrees of freedom, computed from simulated expression data.
                In one experiment <b>V</b> of the m<sub>0</sub> TRUE nulls are rejected
                (false positives), <b>S</b> of the n−m<sub>0</sub> FALSE nulls are
                rejected (true positives).<br>
                <b>FWER</b> := P(V &gt; 0), estimated by the share of experiments with
                at least one false positive.
                <b>Statistical power</b>: the average of S/(n−m<sub>0</sub>), the share of the real
                effects found. <b>Any statistical power</b>: the chance that S &gt; 0.
                <b>False positives (mean)</b> and <b>True positives (mean)</b>: the average
                <i>counts</i> V and S per experiment, not proportions.
                <b>FDR</b>: of the genes a method reports, the average share that are
                false (counted as 0 when it reports nothing).")),
      br(),
      ## the plot gets the full width, with a readable note on how to read it
      fluidRow(
        column(12, h4("Error rate against statistical power"), plotOutput("tradeoff", height = "500px"),
               uiOutput("tradeoff_note"))),
      br(),
      fluidRow(
        column(12, h4("What is going on"),
               uiOutput("reveal_ui"), uiOutput("answer"),
               ## breathing room so the answer box is never flush with the window edge
               div(style = "height:60px;"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## floor(x + 0.5): round() is half-to-even, so e.g. n = 10 at 25% gave m0 = 2
  m0_now <- reactive(floor(input$n * input$pi0 / 100 + 0.5))

  ## the number of pairs depends on n, so the slider is rebuilt when n changes
  ## the slider is declared in the UI so that it always exists, even while the
  ## panel is hidden; only its range follows m0. Rebuilding it with renderUI
  ## suspended it when hidden and reset the value on every m0 change.
  ## a value set by a step's button is kept here until the range has caught up:
  ## sent on its own into the old range, it was snapped to that range's step
  ## (1 became 0 when the step was 8) and the step lost its negative pair
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
  ## and always shows the number of pairs the simulation will use
  observeEvent(input$neg_pairs, {
    req(is.null(pending_neg()), m0_now() >= 2)
    np <- neg_pair_n(m0_now(), flips_for_pairs(m0_now(), input$neg_pairs))
    if (np != input$neg_pairs) updateSliderInput(session, "neg_pairs", value = np)
  }, ignoreInit = TRUE)

  ## k(n - k) moves in steps, so say which count the simulation actually used
  output$neg_text <- renderUI({
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
       (%.0f%% of all pairs are negative)%s.",
      m0, m0, np, if (np == 1) "pair" else "pairs", m0, np, m0, m0, np, k,
      k, if (k == 1) "gene moves" else "genes move", neg_pair_pct(m0, k),
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
                   "Only negatively correlated pairs can break Šidák's and Holm–Šidák's guarantee; Bonferroni and Holm hold regardless.")))
  })

  output$m0_text <- renderUI(
    helpText(sprintf("m₀ = %d TRUE nulls, n − m₀ = %d FALSE nulls",
                     m0_now(), input$n - m0_now())))

  lapply(steps, function(s) {
    observeEvent(input[[paste0("apply_", s$id)]], {
      st <- s$settings
      updateSliderTextInput(session, "n", selected = st$n)
      updateSliderInput(session, "pi0", value = st$pi0)
      updateSliderInput(session, "rho", value = st$rho)
      if (!is.null(st$d)) updateSliderInput(session, "d", value = st$d)
      ## the range, the step and the value go together, so the value cannot be
      ## snapped to the previous range's step
      m0_st  <- floor(st$n * st$pi0 / 100 + 0.5)
      mx_st  <- max_neg_pairs(m0_st)
      neg_st <- if (is.null(st$neg_pairs)) 0 else st$neg_pairs
      updateSliderInput(session, "neg_pairs", max = max(mx_st, 1),
                        step = 1, value = neg_st)
      if (m0_st != m0_now()) pending_neg(neg_st)
      updateRadioButtons(session, "sided",
                         selected = if (is.null(st$sided)) "two" else st$sided)
      updateSliderTextInput(session, "runs",
                            selected = runs_label(if (is.null(st$runs)) 10000 else st$runs))
      updateSliderInput(session, "n_arm",
                        value = if (is.null(st$n_arm)) 10 else st$n_arm)
      updateSliderInput(session, "alpha",
                        value = if (is.null(st$alpha)) 0.05 else st$alpha)
      showNotification("Settings applied — now press Run simulation.",
                       type = "message", duration = 3)
    }, ignoreInit = TRUE)
  })

  sim <- eventReactive(input$run, ignoreNULL = FALSE, {
    runs <- runs_value(input$runs)
    seed <- if (is.null(input$seed) || is.na(input$seed)) NULL else input$seed
    withProgress(message = sprintf("Simulating %s experiments on %s genes...",
                                   runs_label(runs),
                                   format(input$n, big.mark = ",", scientific = FALSE)),
                 value = 0.4, {
      simulate_fwer(n = input$n, m0 = m0_now(), d = input$d, n_arm = input$n_arm,
                    rho = input$rho,
                    neg_pairs = if (input$rho > 0 && input$sided == "one" &&
                                    !is.null(input$neg_pairs))
                                  input$neg_pairs else 0,
                    alpha = input$alpha, n_runs = runs, sided = input$sided,
                    seed = seed)
    })
  })

  output$scenario <- renderText({
    s <- sim()$settings
    sprintf("%d genes, %d TRUE nulls (%.0f%%), %d FALSE nulls | %s test | %s | %s experiments | α = %.2f",
            s$n, s$m0, 100 * s$m0 / s$n, s$m1,
            if (s$sided == "one") "one-sided" else "two-sided",
            if (s$rho == 0) "independent"
            else sprintf("|ρ| = %.2f%s", s$rho,
                         if (s$neg_pairs > 0)
                           sprintf(", %d negative TRUE-null %s", s$neg_pairs,
                                   if (s$neg_pairs == 1) "pair" else "pairs") else ""),
            runs_label(s$runs), s$alpha)
  })

  output$tab <- renderDT({
    d <- sim()$summary
    datatable(d, rownames = FALSE, options = list(dom = "t", ordering = FALSE),
              colnames = c("Method", "FWER", "Statistical power", "Any statistical power",
                           "False positives (mean)", "True positives (mean)",
                           "FDR")) |>
      formatRound(c("FWER", "FDR", "Power", "AnyPower"), 3) |>
      formatRound(c("FalsePositives", "TruePositives"), 2) |>
      formatStyle("FWER", color = crukdark, fontWeight = "bold") |>
      formatStyle("Power", color = crukpink, fontWeight = "bold")
  })

  ## the marginal picture: what one gene's t statistic looks like under each
  ## hypothesis, with the uncorrected rejection region shaded
  ## both top plots read the LAST RUN, not the live sliders: otherwise moving a
  ## slider redraws them while the table below still describes the previous run
  output$marginal <- renderPlot({
    st   <- sim()$settings
    df   <- st$df
    ncp  <- st$effect * sqrt(st$n_arm / 2)
    crit <- if (st$sided == "one") qt(1 - st$alpha, df)
            else qt(1 - st$alpha / 2, df)
    ## take the limits from quantiles, not a fixed margin: with few samples per
    ## group the t tails are heavy and a fixed margin cut the curves off visibly
    lo   <- min(-4.5, qt(5e-4, df))
    hi   <- max(4.5, qt(1 - 5e-4, df), qt(0.995, df, ncp = ncp))
    xs   <- sort(c(seq(lo, hi, length.out = 900), crit, -crit))
    dd   <- rbind(
      data.frame(x = xs, dens = dt(xs, df), Gene = "unchanged (TRUE null)"),
      data.frame(x = xs, dens = dt(xs, df, ncp = ncp),
                 Gene = "differentially expressed"))
    ## with no FALSE nulls there is no effect curve to draw, and with no TRUE
    ## nulls no null curve and no rejection region
    if (st$m1 == 0) dd <- dd[dd$Gene == "unchanged (TRUE null)", ]
    if (st$m0 == 0) dd <- dd[dd$Gene == "differentially expressed", ]
    rej  <- dd[dd$Gene == "unchanged (TRUE null)" &
               (dd$x >= crit | (st$sided == "two" & dd$x <= -crit)), ]
    ## the two tails must be separate polygons: one geom_area over both would
    ## join them across the middle and shade the acceptance region instead
    rej$tail <- ifelse(rej$x > 0, "upper", "lower")
    ggplot(dd, aes(x, dens, colour = Gene, linetype = Gene)) +
      geom_area(data = rej, aes(x, dens, group = tail), inherit.aes = FALSE,
                fill = "grey45", alpha = 0.5) +
      geom_line(linewidth = 1.1) +
      geom_vline(xintercept = if (st$sided == "one") crit else c(-crit, crit),
                 linetype = 3, colour = "grey40") +
      scale_colour_manual(values = c("unchanged (TRUE null)" = crukdark,
                                     "differentially expressed" = crukpink),
                          breaks = c("unchanged (TRUE null)",
                                     "differentially expressed")) +
      scale_linetype_manual(values = c("unchanged (TRUE null)" = "solid",
                                       "differentially expressed" = "longdash"),
                            breaks = c("unchanged (TRUE null)",
                                       "differentially expressed")) +
      labs(x = "t statistic", y = NULL,
           caption = if (st$m0 == 0)
             sprintf("No TRUE nulls in this run. Differentially expressed: noncentral t on %d df", df)
           else sprintf(
             "TRUE null: Student t on %d df%s. Shaded area = \u03b1 = %.2f, the chance one TRUE-null gene is rejected without correction",
             df, if (st$m1 == 0) "" else "; differentially expressed: noncentral t",
             st$alpha)) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top", legend.title = element_blank(),
            axis.text.y = element_blank(), panel.grid.minor = element_blank())
  })

  ## the joint picture: a pair of TRUE-null genes, drawn from the same one-factor
  ## model the simulation uses, so the cloud tilts exactly as |rho| says
  output$bivariate <- renderPlot({
    ## a fixed illustration of the three situations, the same whatever the
    ## sliders say: rho = 0, +0.65 and -0.65, with 30 samples per group so that
    ## the t statistics carry very nearly the genes' correlation
    st  <- list(n_arm = 30)
    r   <- 0.65
    df  <- 2 * (st$n_arm - 1)
    np  <- 500
    ## a still picture, not a new draw per redraw - but restore the RNG stream
    ## afterwards, or the next simulation becomes deterministic
    old_seed <- if (exists(".Random.seed", .GlobalEnv))
                  get(".Random.seed", .GlobalEnv) else NULL
    on.exit(if (is.null(old_seed))
              suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
            else assign(".Random.seed", old_seed, .GlobalEnv))
    set.seed(7)
    ## drawn the same way the simulation draws: raw expression for two genes on
    ## the same samples, then a two-sample t per gene. The label is the GENE-level
    ## correlation.
    N <- 2 * st$n_arm; grp <- rep(c(0, 1), each = st$n_arm)
    pair <- function(sgn, lab) {
      rr  <- if (sgn == 0) 0 else r       # sgn = 0 draws the independent panel
      ld  <- c(1, sign(sgn + (sgn == 0))) * sqrt(rr)
      tt  <- vapply(seq_len(np), function(i) {
        f  <- rnorm(N)
        X  <- outer(ld, f) + sqrt(1 - rr) * matrix(rnorm(2 * N), 2, N)
        mA <- rowMeans(X[, grp == 0, drop = FALSE])
        mB <- rowMeans(X[, grp == 1, drop = FALSE])
        sp <- sqrt((rowSums((X[, grp == 0, drop = FALSE] - mA)^2) +
                    rowSums((X[, grp == 1, drop = FALSE] - mB)^2)) / df)
        (mB - mA) / (sp * sqrt(2 / st$n_arm))
      }, numeric(2))
      data.frame(t1 = tt[1, ], t2 = tt[2, ], Pair = lab)
    }
    ## the three cases together, so the contrast is visible whatever the sliders say
    dd <- rbind(pair(0,  "independent (\u03c1 = 0)"),
                pair(1,  sprintf("same direction (\u03c1 = +%.2f)", r)),
                pair(-1, sprintf("opposite directions (\u03c1 = \u2212%.2f)", r)))
    dd$Pair <- factor(dd$Pair, levels = unique(dd$Pair))
    ## with few samples per group the t tails are heavy enough that a handful of
    ## outliers would squash all three clouds into the same blob, so clip
    lim_t <- max(4, stats::quantile(abs(c(dd$t1, dd$t2)), 0.995))
    ggplot(dd, aes(t1, t2)) +
      geom_hline(yintercept = 0, colour = "grey85") +
      geom_vline(xintercept = 0, colour = "grey85") +
      geom_point(colour = crukdark, alpha = 0.35, size = 1.2) +
      facet_wrap(~ Pair) +
      labs(x = "t statistic, gene i", y = "t statistic, gene j",
           caption = "\u03c1: correlation between the two genes' expression values") +
      coord_equal(xlim = c(-lim_t, lim_t), ylim = c(-lim_t, lim_t)) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  })

  ## the FWER axis of the bar chart: wide enough for alpha and the corrections,
  ## not for No correction, which is often near 1
  bar_cap <- function(d, a)
    min(1.05, max(0.1, 1.3 * a, 1.3 * max(d$FWER[d$Method != "No correction"])))

  ## how to read the plot, in ordinary-sized text rather than a small caption
  output$tradeoff_note <- renderUI({
    d <- sim()$summary; st <- sim()$settings; a <- st$alpha
    ## these are simulated proportions, so a small crossing of the line is noise
    se2 <- 2 * sqrt(a * (1 - a) / st$runs)
    noise <- sprintf("These are simulated rates, so they wobble: with %s experiments a
                      point within about %.3f of the line is <b>not</b> shown to be
                      failing. Raise the number of experiments to sharpen them.",
                     runs_label(st$runs), se2)
    txt <- if (all(is.na(d$Power))) {
      base <- sprintf("No gene is truly differentially expressed here, so there is nothing
               to find and statistical power is not defined. Each bar is one method&#39;s FWER: the
               chance of at least one false claim. The dashed line marks &alpha; = %.2f.
               %s", a, noise)
      unc <- d$FWER[d$Method == "No correction"]
      if (unc > bar_cap(d, a))
        sprintf("%s The <i>No correction</i> bar is cut off by the axis: its value is %.3f.",
                base, unc)
      else base
    } else {
      unc <- d[d$Method == "No correction", ]
      xmax <- max(0.12, max(d$FWER[d$Method != "No correction"]) * 1.25)
      tail <- if (unc$FWER > xmax)
        sprintf("<b><i>No correction</i> is off the right of the plot: FWER %.3f,
                 statistical power %.3f.</b>", unc$FWER, unc$Power)
      else
        sprintf("<i>No correction</i> fits on the plot here (FWER %.3f), so it is shown
                 too.", unc$FWER)
      sprintf("Each point is one method, named in the legend; the shape repeats the colour
               so that overlapping points can still be told apart. <b>Left to right</b>: how
               often it makes at least one false claim (FWER). <b>Bottom to top</b>: the
               share of the real effects it finds (statistical power). The best method is the
               <b>highest point left of the dashed line</b> &mdash; sitting far to the left
               is not a virtue in itself, it is unspent error budget. The dashed line marks
               &alpha; = %.2f. %s %s", a, noise, tail)
    }
    div(style = "font-size:14px; color:#333; margin-top:6px;", HTML(txt))
  })

  output$tradeoff <- renderPlot({
    d <- sim()$summary; st <- sim()$settings

    ## with no true effects there is no power to plot: show the error rates alone
    if (all(is.na(d$Power))) {
      ## No correction can be near 1 and would flatten the four corrections against
      ## the alpha line, so cap the axis -- and keep every label inside the panel,
      ## or a clipped bar loses the number that explains it
      cap <- bar_cap(d, st$alpha)
      d$lab_y <- pmin(d$FWER, cap * 0.94)
      ## a label pulled down inside a clipped bar is written in white
      d$lab_col <- ifelse(d$FWER > cap * 0.94, "white", "grey20")
      return(ggplot(d, aes(Method, FWER, fill = Method)) +
               geom_col(width = 0.7) +
               geom_text(aes(label = sprintf("%.3f", FWER), y = lab_y, colour = lab_col),
                         vjust = -0.6, size = 4) +
               scale_colour_identity() +
               geom_hline(yintercept = st$alpha, linetype = 2, colour = "grey30") +
               scale_fill_manual(values = method_cols) +
               coord_cartesian(ylim = c(0, cap)) +
               labs(x = NULL, y = "FWER") +
               theme_minimal(base_size = 13) +
               theme(axis.text.x = element_text(angle = 20, hjust = 1),
                     legend.position = "none"))
    }

    ## "No correction" usually sits near FWER = 1 and would stretch the axis until
    ## the four corrections overlap in the left margin -- but not always, so drop
    ## it only when it really is off the panel
    xmax <- max(0.12, max(d$FWER[d$Method != "No correction"]) * 1.25)
    dc   <- if (d$FWER[d$Method == "No correction"] > xmax)
              d[d$Method != "No correction", ] else d
    ggplot(dc, aes(FWER, Power, colour = Method, shape = Method)) +
      geom_vline(xintercept = st$alpha, linetype = 2, colour = "grey40") +
      geom_point(size = 3.6, stroke = 1.1) +
      scale_colour_manual(values = method_cols) +
      scale_shape_manual(values = method_shapes) +
      labs(x = "FWER", y = "Statistical power", colour = NULL, shape = NULL) +
      coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
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
      helpText("Choose one of the numbered steps above. Each hides its explanation
                behind a button, so you can run the simulation and read the numbers
                first.")
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
    filename = function() sprintf("fwer-power-%s.csv", format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      res <- sim()
      out <- res$summary
      names(out) <- c("Method", "FWER", "Statistical power", "Any statistical power",
                      "False positives (mean)", "True positives (mean)", "FDR")
      write.csv(cbind(res$settings[rep(1, nrow(out)), ], out, check.names = FALSE),
                file, row.names = FALSE)
    })
}

shinyApp(ui, server)
