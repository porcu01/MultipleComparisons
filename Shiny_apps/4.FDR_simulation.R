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
  BH          = "Benjamini-Hochberg",
  BY          = "Benjamini-Yekutieli"
)
method_cols <- c(
  "No correction"        = "grey45",
  "Bonferroni"           = crukdark,
  "Holm"                 = cruklight,
  "Benjamini-Hochberg"   = crukpink,
  "Benjamini-Yekutieli"  = "#7B2D8E"
)

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
simulate_experiments <- function(N, m0, d, n_arm, rho, neg_share, alpha, q,
                                 n_runs, sided = "two", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m1  <- N - m0
  ncp <- d * sqrt(n_arm / 2)
  df  <- 2 * (n_arm - 1)
  mu  <- c(rep(0, m0), rep(ncp, m1))

  ## r is the strength, neg_share the proportion of genes loading negatively:
  ## corr(z_i, z_j) = s_i * s_j * r, so pairs inside a block correlate +r and
  ## pairs across the two blocks correlate -r.
  r     <- abs(rho)
  n_neg <- round(N * neg_share / 100)
  s     <- rep(1, N)
  ## the sign is assigned at random, independently of which genes are null,
  ## so the two blocks are not confounded with the truth
  if (n_neg > 0) s[sample.int(N, n_neg)] <- -1

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
      padj <- if (mm == "none") p else p.adjust(p, method = mm)
      ## FWER procedures are judged at alpha, FDR procedures at the FDR level q
      lvl  <- if (mm %in% c("BH", "BY")) q else alpha
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
    Level  = ifelse(methods %in% c("BH", "BY"), q, alpha),
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
                                          else if (neg_share == 0) "positive" else "negative",
                             neg_share = neg_share, n_neg = n_neg, rho = r,
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
      tags$li("Compare the FWER and the FDR columns for ", tags$b("No correction"), ".")
    ),
    answer = tags$div(
      tags$p(tags$b("Under the global null the FDR and the FWER are the same thing.")),
      tags$p("With no true effects every rejection is false, so V/R = 1 whenever
              R > 0 and 0 otherwise. The expected proportion of false discoveries
              is then exactly the probability of making at least one: FDR = FWER."),
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
      tags$li("Compare the ", tags$b("Power"), " of Bonferroni and Benjamini-Hochberg."),
      tags$li("Check that the FDR of BH sits near π₀ × q, not at q.")
    ),
    answer = tags$div(
      tags$p(tags$b("The two families of methods control errors in different ways.")),
      tags$p("Bonferroni tests every gene at α/N: the threshold shrinks with the
              number of genes, so real effects are missed. BH compares the i-th
              smallest p-value with i·q/N: the threshold grows along the ranking,
              so a gene with strong evidence is judged against a far milder bar."),
      tags$p("BH controls the FDR at no more than π₀·q, rather than q. With 90% true
              nulls and high power the realised FDR sits close to 0.9q. That slack is the price of not
              knowing π₀ — and the reason Storey's estimator exists.")
    )
  ),
  list(
    id = "step3", title = "3. More genes",
    settings = list(N = 10000, pi0 = 90, dependence = "independent", d = 2),
    look = tags$ul(
      tags$li("Same proportions as step 2, ten times as many genes."),
      tags$li("Watch the ", tags$b("False positives"), " column for the uncorrected row,
               and the BH power and FDR."),
      tags$li("Then set N back to 1000 and compare.")
    ),
    answer = tags$div(
      tags$p(tags$b("The FWER saturates; the FDR does not move.")),
      tags$p("Uncorrected, the FWER is 1 - (1-α)^m₀. With 900 true nulls it is
              already 1 to the last decimal, so it cannot rise any further — what
              grows is the number of false positives, from about 45 to about 450.
              Bonferroni's power collapses as N grows."),
      tags$p("The FDR of BH is a proportion, not a probability. It stays flat at
              about π₀q whatever the number of genes, and the power barely moves.
              That scale-free behaviour is the whole reason the FDR was invented
              for genomics.")
    )
  ),
  list(
    id = "step4", title = "4. Co-regulated genes",
    settings = list(N = 1000, pi0 = 90, dependence = "positive", rho = 0.5, d = 2),
    look = tags$ul(
      tags$li("Genes now share a latent factor: they move up and down together."),
      tags$li("Look at the mean FDR of BH — then look at the spread of the FDP histogram.")
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
      tags$p("The histogram widens sharply. One experiment can have an FDP of zero,
              while the next has an FDP far above q, because the shared factor pushes
              all the genes in the same direction."),
      tags$p("The FDR is an expectation over repeated experiments. It is not a
              guarantee about the experiment on your desk — and you will usually run
              your study only once.")
    )
  ),
  list(
    id = "step5", title = "5. Genes pushing in opposite directions",
    settings = list(N = 1000, pi0 = 90, dependence = "negative", rho = 0.5,
                    neg_share = 50, d = 2, sided = "one"),
    look = tags$ul(
      tags$li("A share of the genes now loads negatively on the shared factor: when
               one block goes up, the other goes down."),
      tags$li("This step switches the test to ", tags$b("one-sided"), ": a two-sided
               test reads |t| and cannot see the sign of a loading, so the error rates
               would look exactly like step 4."),
      tags$li("Compare the FDR of BH with that of Benjamini-Yekutieli, and compare
               their power.")
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
      tags$p(tags$em("Two things worth knowing. First, this is why the test here is
              one-sided: a two-sided test reads |t|, so flipping the sign of a null
              gene's loading leaves the whole distribution of its p-values unchanged,
              and the negative structure would be invisible. Second, even with genuine
              negative dependence BH still comes out conservative on this screen. A
              one-factor model is not a counterexample; the point is that no theorem
              covers this case, not that BH fails here."))
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
  ", crukdark, crukdark, crukpink, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "FWER and FDR in a gene expression experiment")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Each simulated experiment compares <b>two groups</b> of samples across
            <b>N genes</b>, testing every gene for differential expression. Of those genes,
            <b style='color:%s'>m<sub>0</sub> are true nulls</b> (no real difference) and the rest are
            <b style='color:%s'>truly differentially expressed</b>. The experiment is repeated many
            times. We can then measure how often a correction lets <b>at least one</b> false positive
            through (<b style='color:%s'>FWER</b>, a probability) and, <b>averaged over experiments</b>,
            what proportion of the genes it calls significant are false
            (<b style='color:%s'>FDR</b>, an expected proportion).", "grey45", crukpink, crukdark, crukpink))),

  ## top menu: the guided steps
  do.call(tabsetPanel, c(
    list(id = "guide", type = "pills",
         tabPanel("Free exploration",
                  div(class = "step-box",
                      p("Set the parameters on the left and press ", tags$b("Run simulation"),
                        ". Or pick a numbered step above to be guided to a particular
                         phenomenon.")))),
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
      h4("Genes"),
      sliderTextInput("N", "N: number of genes",
                      choices = c(100, 500, 1000, 5000, 10000, 20000),
                      selected = 1000, grid = TRUE),
      sliderInput("pi0", "m₀: TRUE null genes (% of N)",
                  min = 0, max = 100, value = 90, step = 1, post = " %"),
      uiOutput("m0_text"),

      h4("Dependence between genes"),
      radioButtons("dependence", NULL,
                   choices = c("Independent"                    = "independent",
                               "Positive (co-regulated)"        = "positive",
                               "Negative (opposite directions)" = "negative"),
                   selected = "independent"),
      conditionalPanel("input.dependence != 'independent'",
        sliderInput("rho", "Strength of the dependence (|ρ|)",
                    min = 0.05, max = 0.9, value = 0.5, step = 0.05)),
      conditionalPanel("input.dependence == 'negative'",
        sliderInput("neg_share", "Proportion of genes moving in the opposite direction",
                    min = 5, max = 50, value = 50, step = 5, post = " %")),
      uiOutput("dep_text"),

      h4("Experiment"),
      sliderInput("d", "Effect size of the differentially expressed genes (SD units)",
                  min = 0.2, max = 3, value = 2, step = 0.1),
      sliderInput("n_arm", "Samples per group", min = 3, max = 30, value = 10, step = 1),
      radioButtons("sided", "Test",
                   choices = c("Two-sided" = "two",
                               "One-sided (upper tail)" = "one"),
                   selected = "two", inline = TRUE),
      radioButtons("goal", "What do you want to control?",
                   choices = c("FWER only (α)" = "fwer",
                               "FDR only (q)"  = "fdr",
                               "Both, side by side" = "both"),
                   selected = "both"),
      conditionalPanel("input.goal != 'fdr'",
        sliderInput("alpha", "α — level of the FWER methods (Bonferroni, Holm)",
                    min = 0.01, max = 0.10, value = 0.05, step = 0.01)),
      conditionalPanel("input.goal != 'fwer'",
        sliderInput("q", "q — the FDR level you are willing to accept (BH, BY)",
                    min = 0.01, max = 0.25, value = 0.05, step = 0.01)),
      sliderInput("runs", "Number of simulated experiments", min = 50, max = 1000,
                  value = 1000, step = 50),

      downloadButton("save", "Download CSV", style = "width:100%;")
    ),

    ## ---- right: the simulation results -------------------------------------
    mainPanel(
      width = 9,
      tagList(
        fluidRow(
          column(12, h4(textOutput("scenario")), DTOutput("tab"),
                 helpText("Each gene is tested with a two-sample t statistic on
                           2(n-1) degrees of freedom, drawn from its sampling
                           distribution rather than simulated from expression values."),
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

  output$dep_text <- renderUI({
    switch(input$dependence,
      independent = helpText("All genes vary independently."),
      positive    = helpText(sprintf("The shared factor correlates every pair at +%.2f (slightly less for the t statistics, which carry each gene's own variance estimate).", input$rho)),
      negative    = tagList(
        helpText(sprintf("%d%% of the genes load the opposite way: +%.2f within each block, -%.2f between the two blocks.",
                         input$neg_share, input$rho, input$rho)),
        if (input$sided == "two")
          div(style = sprintf("color:%s; font-weight:bold; font-size:12px;", crukpink),
              "With a two-sided test the sign of a loading is invisible to the null genes, so the error rates match positive dependence. Choose a one-sided test to see the negative structure.")))
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
      if (!is.null(st$neg_share)) updateSliderInput(session, "neg_share", value = st$neg_share)
      if (!is.null(st$d))   updateSliderInput(session, "d", value = st$d)
      ## the sidedness of the test is part of a step's setup: step 5 needs a
      ## one-sided test for the negative structure to reach the p-values at all,
      ## and the other steps must be put back to two-sided
      updateRadioButtons(session, "sided",
                         selected = if (is.null(st$sided)) "two" else st$sided)
      showNotification("Settings applied - press Run simulation", type = "message",
                       duration = 3)
    }, ignoreInit = TRUE)
  })

  sim <- eventReactive(input$run, ignoreNULL = FALSE, {
    withProgress(message = "Simulating experiments...", value = 0.4, {
      simulate_experiments(
        N = input$N, m0 = m0_now(), d = input$d, n_arm = input$n_arm,
        rho = if (input$dependence == "independent") 0 else input$rho,
        neg_share = if (input$dependence == "negative") input$neg_share else 0,
        alpha = input$alpha, q = input$q, n_runs = input$runs,
        sided = input$sided
      )
    })
  })

  ## which procedures belong to the chosen goal
  shown <- reactive({
    switch(input$goal,
      fwer = c("No correction", "Bonferroni", "Holm"),
      fdr  = c("No correction", "Benjamini-Hochberg", "Benjamini-Yekutieli"),
      unname(method_labels))
  })

  output$scenario <- renderText({
    s <- sim()$settings
    sprintf("%d genes, %d TRUE nulls (%.0f%%), %d DE genes | %s test | %s dependence%s | %d experiments | α = %.2f, q = %.2f",
            s$N, s$m0, 100 * s$pi0, s$m1, if (s$sided == "one") "one-sided" else "two-sided",
            if (s$rho == 0) "no" else if (s$neg_share == 0) "positive" else "negative",
            if (s$rho == 0) "" else if (s$neg_share == 0) sprintf(" (|ρ| = %.2f)", s$rho)
              else sprintf(" (|ρ| = %.2f, %d%% opposite)", s$rho, s$neg_share),
            s$runs, s$alpha, s$q)
  })

  output$tab <- renderDT({
    d <- sim()$summary
    d <- d[d$Method %in% shown(), , drop = FALSE]
    datatable(d, rownames = FALSE, options = list(dom = "t", ordering = FALSE),
              colnames = c("Method", "Level used", "FWER", "FDR", "Power",
                           "Rejections (mean)", "False positives (mean)")) |>
      formatRound(c("FWER", "Power"), 2) |>
      formatRound("Level", 2) |>
      formatRound("FDR", 3) |>
      formatRound(c("Rejections", "FalsePositives"), 1) |>
      formatStyle("FDR", fontWeight = "bold", color = crukpink) |>
      formatStyle("FWER", color = crukdark)
  })

  output$rates <- renderPlot({
    d <- sim()$summary
    d <- d[d$Method %in% shown(), , drop = FALSE]
    dd <- rbind(data.frame(Method = d$Method, Rate = "FWER", Value = d$FWER),
                data.frame(Method = d$Method, Rate = "FDR",  Value = d$FDR))
    ggplot(dd, aes(Method, Value, fill = Rate)) +
      geom_col(position = position_dodge(0.75), width = 0.7) +
      {if (input$goal != "fdr")
         geom_hline(yintercept = sim()$settings$alpha, linetype = 2, colour = crukdark)} +
      {if (input$goal != "fwer")
         geom_hline(yintercept = sim()$settings$q, linetype = 3, colour = crukpink)} +
      scale_fill_manual(values = c(FWER = crukdark, FDR = crukpink)) +
      labs(x = NULL, y = NULL,
           caption = "dashed: α (FWER methods)    dotted: q (FDR methods)") +
      coord_cartesian(ylim = c(0, 1)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            legend.position = "top", legend.title = element_blank())
  })

  output$fdp <- renderPlot({
    d <- sim()$fdp
    d <- d[d$Method %in% intersect(shown(),
             c("No correction", "Benjamini-Hochberg", "Benjamini-Yekutieli",
               "Bonferroni", "Holm")), ]
    ggplot(d, aes(FDP, fill = Method)) +
      geom_histogram(bins = 30, alpha = 0.65, position = "identity") +
      {if (input$goal != "fwer")
         geom_vline(xintercept = sim()$settings$q, linetype = 2, colour = crukpink)} +
      scale_fill_manual(values = method_cols) +
      labs(x = "False discovery proportion (V / R) in one experiment", y = "experiments",
           caption = "the FDR is the mean of this distribution") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top", legend.title = element_blank())
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
    i <- match(input$guide, step_titles)
    if (is.na(i)) NULL else steps[[i]]
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
      write.csv(cbind(res$settings[rep(1, nrow(res$summary)), ], res$summary),
                file, row.names = FALSE)
    }

  )
}

shinyApp(ui, server)
