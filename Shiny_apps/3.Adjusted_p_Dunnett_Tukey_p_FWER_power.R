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
p_sidak      <- function(p) pmin(1 - (1 - p)^length(p), 1)
p_holm_sidak <- function(p) {
  n <- length(p); o <- order(p); ro <- order(o)
  pmin(cummax(1 - (1 - p[o])^(n - seq_len(n) + 1))[ro], 1)
}

## ---------------------------------------------------------------------------
## One simulated experiment: k groups of n_arm observations.
## Group 1 is the control. `n_eff` of the remaining groups carry a real effect
## of size d (in units of the control SD). Variances may differ between groups,
## and the data may be normal or right-skewed.
## ---------------------------------------------------------------------------
make_data <- function(k, n_arm, d, n_eff, sd_ratio, skew) {
  mu  <- c(0, rep(0, k - 1))
  if (n_eff > 0) mu[seq(k - n_eff + 1, k)] <- d
  sds <- seq(1, sd_ratio, length.out = k)        # 1 when sd_ratio = 1
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
## (Dunnett) or the number of groups (Tukey) still in play.
critical_values <- function(k, n_arm, family, alpha) {
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
      else qmvt(1 - alpha, tail = "both.tails", df = nu,
                corr = matrix(0.5, j, j) + diag(0.5, j))$quantile
    }, numeric(1))
  } else {
    single <- qtukey(1 - alpha, nmeans = k, df = nu) / sqrt(2)
    ## step-down: with g groups still involved, use the range over g means
    stepdown <- vapply(seq_len(m), function(j) {
      g <- max(2, ceiling((1 + sqrt(1 + 8 * j)) / 2))
      qtukey(1 - alpha, nmeans = g, df = nu) / sqrt(2)
    }, numeric(1))
  }
  rank_single <- qmvnorm(1 - alpha, tail = "both.tails", corr = R)$quantile
  list(nu = nu, pairs = pairs, single = single, stepdown = stepdown,
       rank_single = rank_single, m = m)
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
              "Welch + Holm (T2/T3)", "Steel (ranks)"),
  tukey   = c("Bonferroni", "Holm", "Šidák", "Holm–Šidák",
              "Tukey", "Holm–Tukey",
              "Games–Howell", "Steel–Dwass (ranks)"))

method_cols <- c(
  "Bonferroni" = crukdark, "Holm" = cruklight,
  "Šidák" = "#5B8C00", "Holm–Šidák" = "#7B2D8E",
  "Dunnett" = crukpink, "Holm–Dunnett" = "#B3006B",
  "Tukey" = crukpink, "Holm–Tukey" = "#B3006B",
  "Welch + Holm (T2/T3)" = "#E08A00", "Games–Howell" = "#E08A00",
  "Steel (ranks)" = "#00897B", "Steel–Dwass (ranks)" = "#00897B")

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

  ## Welch-based, each comparison with its own variance estimate and df
  w   <- t(apply(pairs, 1, function(ij) welch_pair(y, g, ij[1], ij[2])))
  p_w <- 2 * pt(abs(w[, "t"]), df = w[, "df"], lower.tail = FALSE)
  if (family == "dunnett")
    out[["Welch + Holm (T2/T3)"]] <- p.adjust(p_w, "holm") <= alpha
  else
    out[["Games–Howell"]] <-
      ptukey(sqrt(2) * abs(w[, "t"]), nmeans = k, df = w[, "df"],
             lower.tail = FALSE) <= alpha

  ## rank-based max-t
  z <- apply(pairs, 1, function(ij) rank_pair(y, g, ij[1], ij[2]))
  out[[if (family == "dunnett") "Steel (ranks)" else "Steel–Dwass (ranks)"]] <-
    abs(z) > cv$rank_single

  list(rej = out, is_null = fp$is_null)
}

simulate_mc <- function(k, n_arm, d, n_eff, sd_ratio, skew, family, alpha,
                        n_runs, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cv <- critical_values(k, n_arm, family, alpha)

  runs <- lapply(seq_len(n_runs), function(i)
    one_experiment(k, n_arm, d, n_eff, sd_ratio, skew, family, alpha, cv))

  methods <- names(runs[[1]]$rej)
  is_null <- runs[[1]]$is_null
  m1 <- sum(!is_null)

  summ <- do.call(rbind, lapply(methods, function(mm) {
    V <- vapply(runs, function(r) sum(r$rej[[mm]][r$is_null]), numeric(1))
    S <- vapply(runs, function(r) sum(r$rej[[mm]][!r$is_null]), numeric(1))
    data.frame(Method = mm, FWER = mean(V > 0),
               Power = if (m1 > 0) mean(S / m1) else NA_real_,
               FalsePositives = mean(V), TruePositives = mean(S))
  }))
  summ$Method <- factor(summ$Method, levels = methods)

  list(summary = summ,
       settings = data.frame(k = k, n_arm = n_arm, comparisons = cv$m,
                             effect = d, n_eff = n_eff, sd_ratio = sd_ratio,
                             skew = skew, family = family, alpha = alpha,
                             runs = n_runs, true_nulls = sum(is_null)))
}

## ---------------------------------------------------------------------------
steps <- list(
  list(id = "s1", title = "1. Using the design beats correcting blindly",
       settings = list(k = 5, n_eff = 0, sd_ratio = 1, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("No group differs: every comparison is a true null."),
         tags$li("Compare the ", tags$b("FWER"), " of Bonferroni with that of Dunnett or Tukey.")),
       answer = tags$div(
         tags$p(tags$b("All of them control the FWER, but not equally tightly.")),
         tags$p("Bonferroni and Holm treat the comparisons as if nothing were known
                about them. Dunnett and Tukey use the fact that the comparisons share
                groups, so their statistics are correlated — +1/2 between any two
                Dunnett contrasts in a balanced design — and the critical value is
                read off that joint distribution."),
         tags$p("Because they use the real joint distribution, Dunnett and Tukey land
                almost exactly on α, while Bonferroni sits below it. That unused
                margin is wasted power, which the next step makes visible."))),

  list(id = "s2", title = "2. The power that correlation buys",
       settings = list(k = 5, n_eff = 2, sd_ratio = 1, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("Two of the groups now carry a real effect."),
         tags$li("Compare ", tags$b("Power"), " down the column."),
         tags$li("Then compare Dunnett with Holm–Dunnett (or Tukey with Holm–Tukey).")),
       answer = tags$div(
         tags$p(tags$b("Two separate gains, and they add up.")),
         tags$p("The first gain is the joint distribution: Dunnett and Tukey are more
                powerful than Bonferroni and Šidák because they stop correcting
                for multiplicity that the correlation has already removed."),
         tags$p("The second is the step-down. Holm–Dunnett and Holm–Tukey apply the
                same max-t idea sequentially, relaxing the critical value as
                hypotheses are rejected. They reject everything the single-step
                version rejects and sometimes more — again for free."))),

  list(id = "s3", title = "3. When the variances differ",
       settings = list(k = 5, n_eff = 2, sd_ratio = 4, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("The groups now have different variances (the last has 4× the SD of the first)."),
         tags$li("Watch the ", tags$b("FWER"), " of Dunnett and Tukey."),
         tags$li("Compare them with the Welch-based method (T2/T3 or Games–Howell).")),
       answer = tags$div(
         tags$p(tags$b("Equal variances is a real assumption, not a formality.")),
         tags$p("Dunnett and Tukey pool all the groups into one variance estimate.
                When the variances actually differ, that pooled estimate is wrong for
                each individual comparison — too large for the tight groups and too
                small for the spread-out ones — and the error rate drifts away from
                α, often above it."),
         tags$p("The Welch-based procedures give every comparison its own variance
                estimate and its own Satterthwaite degrees of freedom. They cost a
                little power when the variances really are equal, and they are the
                safe choice when you do not know."))),

  list(id = "s4", title = "4. When the data are not normal",
       settings = list(k = 5, n_eff = 2, sd_ratio = 1, skew = TRUE, d = 1.5),
       look = tags$ul(
         tags$li("The data are now strongly right-skewed rather than normal."),
         tags$li("Compare the rank-based method with the t-based ones, on both
                  ", tags$b("FWER"), " and ", tags$b("Power"), ".")),
       answer = tags$div(
         tags$p(tags$b("The t-based tests survive skewness better than people expect.")),
         tags$p("With equal group sizes the t statistic is fairly robust to skewness,
                so Dunnett and Tukey usually keep their error rate near α here.
                What suffers is power: the mean is a poor summary of a skewed
                distribution, so a shift of a given size is harder to detect."),
         tags$p("The rank-based procedures — Steel, Steel–Dwass, Dunn — make no
                normality assumption and often gain power on skewed data. They use
                the same max-t logic, with the joint distribution of the rank
                statistics; Dunn is the exception, correcting its p-values one at a
                time instead."))),

  list(id = "s5", title = "5. Many groups",
       settings = list(k = 10, n_eff = 2, sd_ratio = 1, skew = FALSE, d = 1.5),
       look = tags$ul(
         tags$li("Ten groups. For Tukey that is 45 comparisons; for Dunnett, 9."),
         tags$li("Switch the family on the left and compare the two."),
         tags$li("Look at how far Bonferroni now falls behind Tukey.")),
       answer = tags$div(
         tags$p(tags$b("The more comparisons share groups, the more the design is worth.")),
         tags$p("With 45 pairwise comparisons among 10 groups, Bonferroni divides α
                by 45. But those 45 statistics are built from only 10 group means, so
                they are heavily correlated and there are nowhere near 45 independent
                questions. Tukey accounts for that; Bonferroni cannot."),
         tags$p("Choosing the right family matters just as much. If only the
                comparisons against a control are of interest, asking for all 45
                pairs throws away power for questions you never meant to ask — which
                is the whole reason Dunnett's test exists.")))
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
  ", crukdark, crukdark, crukpink, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "Dunnett and Tukey compared with Bonferroni and Holm: FWER and power")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Several groups are compared in one experiment — either each arm
            against a <b>control</b> (Dunnett's family) or <b>every pair</b> (Tukey's family).
            Because the comparisons share groups, their test statistics are correlated,
            and a procedure that uses that correlation needs a smaller correction than one
            that does not. The table reports how often each method makes
            <b style='color:%s'>any false claim</b> and what share of the real differences
            it finds.", crukdark))),

  do.call(tabsetPanel, c(
    list(id = "guide", type = "pills",
         tabPanel("Free exploration",
                  div(class = "step-box",
                      HTML("Set the design on the left and press <b>Run simulation</b>.
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
      h4("Which comparisons"),
      radioButtons("family", NULL,
                   choices = c("Each arm against a control (Dunnett)" = "dunnett",
                               "Every pair of groups (Tukey)" = "tukey"),
                   selected = "dunnett"),
      uiOutput("fam_text"),

      h4("The design"),
      sliderInput("k", "Number of groups", min = 3, max = 10, value = 5, step = 1),
      sliderInput("n_arm", "Observations per group", min = 4, max = 30, value = 10, step = 1),
      sliderInput("n_eff", "Groups carrying a real effect", min = 0, max = 9, value = 2, step = 1),
      sliderInput("d", "Effect size (SD units)", min = 0.2, max = 3, value = 1.5, step = 0.1),

      h4("Assumptions"),
      sliderInput("sd_ratio", "Largest / smallest SD  (1 = equal variances)",
                  min = 1, max = 6, value = 1, step = 0.5),
      checkboxInput("skew", "Right-skewed data instead of normal", value = FALSE),

      h4("Simulation"),
      sliderInput("alpha", "α — the family-wise error rate",
                  min = 0.01, max = 0.10, value = 0.05, step = 0.01),
      sliderInput("runs", "Number of simulated experiments",
                  min = 50, max = 1000, value = 300, step = 50),
      helpText("Tukey with many groups is the slow case: 300 runs take a few seconds."),
      downloadButton("save", "Download CSV", style = "width:100%;")
    ),

    mainPanel(
      width = 9,
      h4(textOutput("scenario")),
      DTOutput("tab"),
      helpText("Dunnett, Tukey and the rank-based methods read their critical value
                from the joint distribution of the comparisons. Bonferroni, Šidák
                and Holm treat the comparisons one at a time."),
      br(),
      fluidRow(
        column(6, h4("Error rate against power"), plotOutput("tradeoff", height = "320px")),
        column(6, h4("What is going on"), uiOutput("reveal_ui"), uiOutput("answer"))
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
      showNotification("Settings applied — now press 'Run simulation'.",
                       type = "message", duration = 3)
    }, ignoreInit = TRUE)
  })

  sim <- eventReactive(input$run, ignoreNULL = FALSE, {
    withProgress(message = "Simulating experiments...", value = 0.4, {
      simulate_mc(k = input$k, n_arm = input$n_arm, d = input$d,
                  n_eff = min(input$n_eff, input$k - 1),
                  sd_ratio = input$sd_ratio, skew = input$skew,
                  family = input$family, alpha = input$alpha, n_runs = input$runs)
    })
  })

  output$scenario <- renderText({
    s <- sim()$settings
    sprintf("%d groups of %d | %s: %d comparisons, %d of them true nulls | %s | %s | %d experiments | α = %.2f",
            s$k, s$n_arm,
            if (s$family == "dunnett") "against a control" else "all pairs",
            s$comparisons, s$true_nulls,
            if (s$sd_ratio == 1) "equal variances" else sprintf("SD ratio %.1f", s$sd_ratio),
            if (s$skew) "skewed data" else "normal data",
            s$runs, s$alpha)
  })

  output$tab <- renderDT({
    d <- sim()$summary
    datatable(d, rownames = FALSE, options = list(dom = "t", ordering = FALSE),
              colnames = c("Method", "FWER", "Power",
                           "False positives (mean)", "True positives (mean)")) |>
      formatRound(c("FWER", "Power"), 3) |>
      formatRound(c("FalsePositives", "TruePositives"), 2) |>
      formatStyle("FWER", color = crukdark, fontWeight = "bold") |>
      formatStyle("Power", color = crukpink, fontWeight = "bold")
  })

  output$tradeoff <- renderPlot({
    d <- sim()$summary
    a <- sim()$settings$alpha
    if (all(is.na(d$Power)))
      return(ggplot(d, aes(Method, FWER, fill = Method)) +
               geom_col(width = 0.7) +
               geom_hline(yintercept = a, linetype = 2, colour = "grey30") +
               scale_fill_manual(values = method_cols) +
               labs(x = NULL, y = "FWER",
                    caption = "no real effects here, so power is undefined") +
               theme_minimal(base_size = 12) +
               theme(axis.text.x = element_text(angle = 25, hjust = 1),
                     legend.position = "none"))
    ggplot(d, aes(FWER, Power, colour = Method)) +
      geom_vline(xintercept = a, linetype = 2, colour = "grey40") +
      geom_point(size = 4) +
      geom_text(aes(label = Method), vjust = -1.0, size = 3.2, show.legend = FALSE) +
      scale_colour_manual(values = method_cols) +
      coord_cartesian(xlim = c(0, max(0.12, max(d$FWER) * 1.15)), ylim = c(0, 1.05)) +
      labs(x = "FWER", y = "power",
           caption = "dashed: the nominal α. Methods to the right of it are not controlling the FWER.") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  })

  current_step <- reactive({
    i <- match(input$guide, step_titles)
    if (is.na(i)) NULL else steps[[i]]
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
