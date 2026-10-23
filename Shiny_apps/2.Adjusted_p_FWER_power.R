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

adjust_p <- function(p, method) {
  n <- length(p)
  if (method == "none") p
  else if (method == "sidak") pmin(1 - (1 - p)^n, 1)
  else if (method == "holm_sidak") {
    o <- order(p); ro <- order(o)
    pmin(cummax(1 - (1 - p[o])^(n - seq_len(n) + 1))[ro], 1)
  } else p.adjust(p, method = method)
}

## ---------------------------------------------------------------------------
## One experiment = n hypotheses tested in two groups of n_arm samples.
## m0 are TRUE nulls. Statistics are two-sample t: a numerator z_i ~ N(mu_i, 1)
## divided by that test's own variance estimate, t_i = z_i / sqrt(w_i),
## w_i ~ chi^2_df / df, df = 2(n_arm - 1).
## Dependence enters the numerators through one latent factor:
##     z_i = mu_i + s_i * sqrt(|rho|) * F + sqrt(1 - |rho|) * e_i
## so corr(z_i, z_j) = s_i s_j |rho|.
## ---------------------------------------------------------------------------
simulate_fwer <- function(n, m0, d, n_arm, rho, neg_share, alpha, n_runs,
                          sided = "two", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m1  <- n - m0
  ncp <- d * sqrt(n_arm / 2)
  df  <- 2 * (n_arm - 1)
  mu  <- c(rep(0, m0), rep(ncp, m1))

  r     <- abs(rho)
  n_neg <- round(n * neg_share / 100)
  s     <- rep(1, n)
  if (n_neg > 0) s[sample.int(n, n_neg)] <- -1

  methods <- names(method_labels)
  per_run <- lapply(seq_len(n_runs), function(i) {
    f  <- if (r > 0) rnorm(1) else 0
    z  <- mu + s * sqrt(r) * f + sqrt(1 - r) * rnorm(n)
    tt <- z / sqrt(rchisq(n, df) / df)
    p  <- if (sided == "one") pt(tt, df = df, lower.tail = FALSE)
          else 2 * pt(abs(tt), df = df, lower.tail = FALSE)

    vapply(methods, function(mm) {
      rej <- adjust_p(p, mm) <= alpha
      V <- sum(rej[seq_len(m0)])
      S <- if (m1 > 0) sum(rej[(m0 + 1):n]) else 0
      c(V = V, S = S,
        any_false = as.numeric(V > 0),
        power     = if (m1 > 0) S / m1 else NA_real_,
        any_power = if (m1 > 0) as.numeric(S > 0) else NA_real_)
    }, numeric(5))
  })

  arr <- simplify2array(per_run)
  summ <- data.frame(
    Method = unname(method_labels[methods]),
    FWER   = apply(arr["any_false", , , drop = FALSE], 2, mean),
    Power  = apply(arr["power", , , drop = FALSE], 2, mean),
    AnyPower = apply(arr["any_power", , , drop = FALSE], 2, mean),
    FalsePositives = apply(arr["V", , , drop = FALSE], 2, mean),
    TruePositives  = apply(arr["S", , , drop = FALSE], 2, mean),
    row.names = NULL, stringsAsFactors = FALSE
  )
  summ$Method <- factor(summ$Method, levels = unname(method_labels))
  list(summary = summ,
       settings = data.frame(n = n, m0 = m0, m1 = m1, effect = d, n_arm = n_arm,
                             df = df, rho = r, neg_share = neg_share,
                             sided = sided, alpha = alpha, runs = n_runs))
}

## ---------------------------------------------------------------------------
steps <- list(
  list(id = "s1", title = "1. The error rate grows with the family",
       settings = list(n = 20, pi0 = 100, rho = 0, d = 2),
       look = tags$ul(
         tags$li("Every hypothesis is a TRUE null, so any rejection is a mistake."),
         tags$li("Read the ", tags$b("FWER"), " of ", tags$b("No correction"),
                 " and compare it with α."),
         tags$li("Then raise n from 20 to 100 and run again.")),
       answer = tags$div(
         tags$p(tags$b("Testing at α per hypothesis does not hold the family at α.")),
         tags$p("With m₀ independent true nulls the chance of at least one false
                positive is 1 − (1 − α)^m₀: about 0.64 at n = 20 and 0.99 at
                n = 100, for α = 0.05. All four corrections sit at or below α
                instead, which is what controlling the FWER means."))),

  list(id = "s2", title = "2. Holm beats Bonferroni, for free",
       settings = list(n = 50, pi0 = 80, rho = 0, d = 2),
       look = tags$ul(
         tags$li("Compare the ", tags$b("FWER"), " of Bonferroni and Holm — both at or below α."),
         tags$li("Now compare their ", tags$b("Power"), ".")),
       answer = tags$div(
         tags$p(tags$b("Same guarantee, more rejections.")),
         tags$p("Bonferroni tests every hypothesis at α/n. Holm uses the same
                threshold for the smallest p-value, then relaxes it to α/(n−1),
                α/(n−2) and so on as hypotheses are rejected."),
         tags$p("So Holm rejects everything Bonferroni rejects and sometimes more:
                it is uniformly more powerful, at no cost in error control and
                under no extra assumption. There is never a reason to prefer
                Bonferroni on statistical grounds."))),

  list(id = "s3", title = "3. What Šidák adds, and what it assumes",
       settings = list(n = 50, pi0 = 80, rho = 0, d = 2),
       look = tags$ul(
         tags$li("Compare ", tags$b("Šidák"), " with ", tags$b("Bonferroni"),
                 ", and ", tags$b("Holm–Šidák"), " with ", tags$b("Holm"), "."),
         tags$li("The differences are small — look at the third decimal of Power.")),
       answer = tags$div(
         tags$p(tags$b("A small gain, bought with an assumption.")),
         tags$p("Šidák tests at 1 − (1−α)^(1/n) instead of α/n. That is
                slightly larger — 0.00512 against 0.00500 at n = 10, α = 0.05 —
                so Šidák is uniformly a little more powerful than Bonferroni, and
                Holm–Šidák a little more powerful than Holm."),
         tags$p("The price is that the gain is justified by assuming the tests are
                independent, whereas Bonferroni and Holm need no such assumption.
                In exchange for a third-decimal improvement you take on a condition
                you usually cannot check. Step 4 shows what happens when it fails."))),

  list(id = "s4", title = "4. When the tests are not independent",
       settings = list(n = 50, pi0 = 80, rho = 0.6, d = 2),
       look = tags$ul(
         tags$li("The tests now share a latent factor, so they move together."),
         tags$li("Watch the ", tags$b("FWER"), " of all four corrections."),
         tags$li("Compare the power with step 3, where the tests were independent.")),
       answer = tags$div(
         tags$p(tags$b("Positive dependence makes every correction conservative.")),
         tags$p("Bonferroni and Holm rest on Boole's inequality, which holds whatever
                the dependence, so they remain valid — just further below α than
                before. Šidák's independence assumption is violated, but for
                two-sided tests on this kind of positively dependent structure it
                stays valid too, and also becomes conservative."),
         tags$p("The practical consequence is loss of power: when the tests overlap,
                there are effectively fewer than n independent questions, and all
                four methods are correcting for more multiplicity than is really
                there."))),

  list(id = "s5", title = "5. Many tests, few real effects",
       settings = list(n = 200, pi0 = 95, rho = 0, d = 2),
       look = tags$ul(
         tags$li("200 hypotheses, of which 10 are truly false."),
         tags$li("Look at ", tags$b("Power"), " — the share of the real effects found."),
         tags$li("Then look at ", tags$b("Any power"), " — the chance of finding at least one.")),
       answer = tags$div(
         tags$p(tags$b("FWER control is expensive when the family is large.")),
         tags$p("Guarding against even one false positive across 200 tests forces a
                threshold near α/200. Real effects of moderate size are then missed
                most of the time, and the power column falls accordingly."),
         tags$p("This is the Achilles' heel of the FWER, and it is why the false
                discovery rate was invented: if you can accept that a controlled
                proportion of your findings are wrong, the thresholds relax
                enormously. That is the subject of app 4.")))
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
                 "FWER and power: Bonferroni, Holm, Šidák, Holm–Šidák")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Each simulated experiment tests <b>n hypotheses</b> in two groups of
            samples. Of them, <b style='color:grey'>m<sub>0</sub> are true nulls</b> and the
            rest are <b style='color:%s'>truly false</b>. Repeating the experiment many times
            shows how often each correction lets <b>at least one</b> false positive through
            (<b style='color:%s'>FWER</b>) and what share of the real effects it finds
            (<b>power</b>).", crukpink, crukdark))),

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
      h4("The family"),
      sliderTextInput("n", "n: number of hypotheses",
                      choices = c(5, 10, 20, 50, 100, 200, 500), selected = 50, grid = TRUE),
      sliderInput("pi0", "TRUE nulls (% of n)", min = 0, max = 100, value = 80, step = 5, post = " %"),
      uiOutput("m0_text"),

      h4("Dependence between tests"),
      sliderInput("rho", "Strength |ρ|  (0 = independent)",
                  min = 0, max = 0.9, value = 0, step = 0.05),
      conditionalPanel("input.rho > 0",
        sliderInput("neg_share", "Proportion loading the opposite way",
                    min = 0, max = 50, value = 0, step = 5, post = " %")),

      h4("Experiment"),
      sliderInput("d", "Effect size of the false nulls (SD units)",
                  min = 0.2, max = 3, value = 2, step = 0.1),
      sliderInput("n_arm", "Samples per group", min = 3, max = 30, value = 10, step = 1),
      radioButtons("sided", "Test", choices = c("Two-sided" = "two",
                                                "One-sided (upper tail)" = "one"),
                   selected = "two", inline = TRUE),
      sliderInput("alpha", "α — the family-wise error rate",
                  min = 0.01, max = 0.10, value = 0.05, step = 0.01),
      sliderInput("runs", "Number of simulated experiments",
                  min = 50, max = 2000, value = 1000, step = 50),
      downloadButton("save", "Download CSV", style = "width:100%;")
    ),

    mainPanel(
      width = 9,
      h4(textOutput("scenario")),
      DTOutput("tab"),
      helpText("Each test is a two-sample t statistic on 2(n−1) degrees of freedom,
                drawn from its sampling distribution rather than simulated from data.
                \"Any power\" is the chance of finding at least one of the real effects."),
      br(),
      fluidRow(
        column(6, h4("Error rate against power"), plotOutput("tradeoff", height = "300px")),
        column(6, h4("What is going on"),
               uiOutput("reveal_ui"), uiOutput("answer"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  m0_now <- reactive(round(input$n * input$pi0 / 100))

  output$m0_text <- renderUI(
    helpText(sprintf("m₀ = %d true nulls, m₁ = %d truly false",
                     m0_now(), input$n - m0_now())))

  lapply(steps, function(s) {
    observeEvent(input[[paste0("apply_", s$id)]], {
      st <- s$settings
      updateSliderTextInput(session, "n", selected = st$n)
      updateSliderInput(session, "pi0", value = st$pi0)
      updateSliderInput(session, "rho", value = st$rho)
      if (!is.null(st$d)) updateSliderInput(session, "d", value = st$d)
      updateSliderInput(session, "neg_share",
                        value = if (is.null(st$neg_share)) 0 else st$neg_share)
      updateRadioButtons(session, "sided",
                         selected = if (is.null(st$sided)) "two" else st$sided)
      showNotification("Settings applied — now press 'Run simulation'.",
                       type = "message", duration = 3)
    }, ignoreInit = TRUE)
  })

  sim <- eventReactive(input$run, ignoreNULL = FALSE, {
    withProgress(message = "Simulating experiments...", value = 0.4, {
      simulate_fwer(n = input$n, m0 = m0_now(), d = input$d, n_arm = input$n_arm,
                    rho = input$rho,
                    neg_share = if (input$rho > 0) input$neg_share else 0,
                    alpha = input$alpha, n_runs = input$runs, sided = input$sided)
    })
  })

  output$scenario <- renderText({
    s <- sim()$settings
    sprintf("%d hypotheses, %d true nulls (%.0f%%), %d false | %s test | %s | %d experiments | α = %.2f",
            s$n, s$m0, 100 * s$m0 / s$n, s$m1,
            if (s$sided == "one") "one-sided" else "two-sided",
            if (s$rho == 0) "independent"
            else sprintf("|ρ| = %.2f%s", s$rho,
                         if (s$neg_share > 0) sprintf(", %d%% opposite", s$neg_share) else ""),
            s$runs, s$alpha)
  })

  output$tab <- renderDT({
    d <- sim()$summary
    datatable(d, rownames = FALSE, options = list(dom = "t", ordering = FALSE),
              colnames = c("Method", "FWER", "Power", "Any power",
                           "False positives (mean)", "True positives (mean)")) |>
      formatRound(c("FWER", "Power", "AnyPower"), 3) |>
      formatRound(c("FalsePositives", "TruePositives"), 2) |>
      formatStyle("FWER", color = crukdark, fontWeight = "bold") |>
      formatStyle("Power", color = crukpink, fontWeight = "bold")
  })

  output$tradeoff <- renderPlot({
    d <- sim()$summary
    ## with no true effects there is no power to plot: show the error rates alone
    if (all(is.na(d$Power)))
      return(ggplot(d, aes(Method, FWER, fill = Method)) +
               geom_col(width = 0.7) +
               geom_hline(yintercept = sim()$settings$alpha, linetype = 2, colour = "grey30") +
               scale_fill_manual(values = method_cols) +
               labs(x = NULL, y = "FWER",
                    caption = "no true effects in this scenario, so power is undefined") +
               theme_minimal(base_size = 12) +
               theme(axis.text.x = element_text(angle = 20, hjust = 1),
                     legend.position = "none"))
    ggplot(d, aes(FWER, Power, colour = Method)) +
      geom_vline(xintercept = sim()$settings$alpha, linetype = 2, colour = "grey40") +
      geom_point(size = 4) +
      geom_text(aes(label = Method), vjust = -1.1, size = 3.4, show.legend = FALSE) +
      scale_colour_manual(values = method_cols) +
      labs(x = "FWER", y = "power",
           caption = "dashed: the nominal α. Good methods sit left of it and as high as possible.") +
      coord_cartesian(xlim = c(0, max(0.12, max(d$FWER) * 1.05)), ylim = c(0, 1)) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right", legend.title = element_blank())
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
    filename = function() sprintf("fwer-power-%s.csv", format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      res <- sim()
      write.csv(cbind(res$settings[rep(1, nrow(res$summary)), ], res$summary),
                file, row.names = FALSE)
    })
}

shinyApp(ui, server)
