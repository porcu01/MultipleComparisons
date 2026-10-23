library(shiny)
library(DT)
library(ggplot2)

crukdark  <- "#00007E"
cruklight <- "#009CEE"
crukpink  <- "#FF0087"

## ---------------------------------------------------------------------------
## Adjusted p-values
##
## p.adjust covers Bonferroni and Holm. Sidak and Holm-Sidak are not in
## p.adjust, so they are written out here:
##   Sidak       (single-step):  1 - (1 - p_i)^n
##   Holm-Sidak  (step-down):    running max of 1 - (1 - p_(i))^(n - i + 1)
## The running maximum is what enforces monotonicity of the adjusted p-values,
## i.e. the step-down stopping rule.
## ---------------------------------------------------------------------------

## 1 - (1 - p)^k written so that it stays accurate for tiny p (1 - p would round
## to 1 below about 1e-16), and returns p itself when k = 1
sidak_adj <- function(p, k) ifelse(k == 1, p, -expm1(k * log1p(-p)))

## round away floating-point noise, so that 3 x 0.05 is exactly 0.15
clean <- function(x) signif(x, 12)

adjust_p <- function(p, method) {
  n <- length(p)
  clean(switch(method,
    bonferroni = pmin(1, n * p),
    holm       = p.adjust(p, method = "holm"),
    sidak      = sidak_adj(p, n),
    holm_sidak = {
      o  <- order(p); ro <- order(o)
      pmin(cummax(sidak_adj(p[o], n - seq_len(n) + 1))[ro], 1)
    }
  ))
}

## a hypothesis is rejected when its adjusted p-value is at most alpha; the tiny
## tolerance keeps a value exactly on the level from being retained by rounding
rejects <- function(adj, alpha) adj <= alpha * (1 + 1e-9)

methods_all <- c(
  "Bonferroni"          = "bonferroni",
  "Holm"                = "holm",
  "Šidák"               = "sidak",
  "Holm–Šidák"          = "holm_sidak"
)
single_step <- c("bonferroni", "sidak")

## how each method builds its adjusted p-value, shown step by step
step_formula <- function(method, n) {
  switch(method,
    bonferroni = sprintf("p̃ᵢ = min(1, %d × pᵢ)", n),
    holm       = sprintf("p̃(ᵢ) = min(1, max over j ≤ i of  (%d − j + 1) × p(ⱼ))", n),
    sidak      = sprintf("p̃ᵢ = 1 − (1 − pᵢ)^%d", n),
    holm_sidak = sprintf("p̃(ᵢ) = max over j ≤ i of  1 − (1 − p(ⱼ))^(%d − j + 1)", n)
  )
}

## the threshold each ordered p-value is compared with (the rule on the raw p)
step_threshold <- function(method, n, i, alpha) {
  clean(switch(method,
    bonferroni = alpha / n,
    holm       = alpha / (n - i + 1),
    sidak      = 1 - (1 - alpha)^(1 / n),
    holm_sidak = 1 - (1 - alpha)^(1 / (n - i + 1))
  ))
}

fmt <- function(x) sprintf("%.4g", x)

## the arithmetic of a single step, written out with the actual numbers:
## the true product first, then any cap at 1, then any running maximum
step_calc <- function(method, n, i, p, adj) {
  k <- if (method %in% single_step) n else n - i + 1
  if (method %in% c("bonferroni", "holm")) {
    raw  <- k * p
    head <- if (method == "bonferroni") sprintf("%d × %s", n, fmt(p))
            else sprintf("(%d − %d + 1) × %s", n, i, fmt(p))
    s <- sprintf("%s = %s", head, fmt(raw))
    if (raw > 1) s <- paste0(s, ", capped at 1")
    own <- clean(min(1, raw))
  } else {
    own  <- clean(sidak_adj(p, k))
    head <- if (method == "sidak") sprintf("1 − (1 − %s)^%d", fmt(p), n)
            else sprintf("1 − (1 − %s)^(%d − %d + 1)", fmt(p), n, i)
    s <- sprintf("%s = %s", head, fmt(own))
  }
  if (adj > own * (1 + 1e-9))
    s <- paste0(s, sprintf(", raised to the previous adjusted value, %s (adjusted p-values cannot decrease)",
                           fmt(adj)))
  s
}

## read the box into tokens and sort them into valid p-values and problems
number_re <- "^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$"
parse_input <- function(txt) {
  tok <- strsplit(trimws(txt), "[,;[:space:]]+")[[1]]
  tok <- tok[nzchar(tok)]
  is_num <- grepl(number_re, tok)
  val <- rep(NA_real_, length(tok))
  val[is_num] <- as.numeric(tok[is_num])
  list(tok = tok, val = val,
       not_number   = tok[!is_num],
       out_of_range = tok[is_num & (val < 0 | val > 1)])
}

## a p-value written with a decimal comma (0,05) would be split into 0 and 5,
## so flag "0," followed by digits that are not themselves a decimal number
has_decimal_comma <- function(txt) {
  grepl("(^|[^0-9.eE+-])0,[0-9]+($|[^0-9.])", txt)
}

plural <- function(k, one, many) if (k == 1) one else many

## tables show 50 rows per page; beyond that add the row count and page buttons
rows_per_page <- 50
table_dom <- function(d) if (nrow(d) > rows_per_page) "tip" else "t"

example_sets <- list(
  "Two hypotheses (Thuesen data regression)" = "0.00000001, 0.048",
  "Ten genes, a few small p-values"          = "0.0001, 0.0021, 0.013, 0.021, 0.036, 0.11, 0.24, 0.38, 0.61, 0.92",
  "Twenty p-values"                          = paste(sprintf("%.3f", sort(c(0.012, 0.031, 0.049, 0.067, 0.091,
                                                                           0.12, 0.18, 0.21, 0.26, 0.33,
                                                                           0.37, 0.44, 0.51, 0.58, 0.63,
                                                                           0.70, 0.77, 0.84, 0.90, 0.97))),
                                                     collapse = ", ")
)

## ---------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML(sprintf("
    .well { background:#F7FAFF; border-color:%s33; }
    h4 { color:%s; }
    .note { border-left:5px solid %s; padding:6px 12px; background:#F2F9FF; }
  ", crukdark, crukdark, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "Adjusted p-values: how each correction is calculated")),

  div(style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("Enter the raw p-values of a family of tests. For each test, the
            <b>adjusted</b> p-value under each correction is reported, and the arithmetic
            that produces it is shown step by step for the correction you choose. A
            hypothesis is rejected when its adjusted p-value is at most the level α you
            choose, so each correction's decisions can be read from its own column. All
            four corrections control the <b style='color:%s'>family-wise error rate
            (FWER)</b>, the probability of <i>at least one</i> false positive anywhere in
            the family, at level α: Bonferroni and Holm whatever the dependence between
            the tests, Šidák and Holm–Šidák when the tests are independent or positively
            dependent (positive orthant dependence, POD).", crukdark))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("The family of tests"),
      selectInput("example", "Start from an example",
                  choices = c("(type my own)", names(example_sets)),
                  selected = names(example_sets)[2]),
      textAreaInput("pvals", "Raw p-values",
                    value = example_sets[[2]], rows = 6,
                    placeholder = "0.001, 0.02, 0.3 ..."),
      helpText("Separate with commas, spaces or new lines. Use a point for decimals (0.05)."),
      uiOutput("n_text"),

      h4("Decision level"),
      sliderInput("alpha", "α — the level at which the FWER is controlled",
                  min = 0.01, max = 0.25, value = 0.05, step = 0.01),

      h4("Show the arithmetic for"),
      selectInput("focus", NULL, choices = methods_all, selected = "holm"),
      downloadButton("save", "Download table", style = "width:100%;")
    ),

    mainPanel(
      width = 9,
      h4(textOutput("headline")),
      DTOutput("tab"),
      helpText("Rejected: adjusted p-value at most α, in bold dark blue. Not rejected:
                greyed out. Rows are sorted by raw p-value."),
      br(),
      fluidRow(
        column(7,
               h4(textOutput("focus_title")),
               div(class = "note", uiOutput("focus_rule")),
               br(),
               DTOutput("steps")),
        column(5,
               h4("Adjusted against raw p-values"),
               plotOutput("plot", height = "300px"),
               checkboxInput("logscale", "Log scale (spreads out the small p-values)", value = TRUE),
               uiOutput("count_text"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$example, {
    if (input$example %in% names(example_sets))
      updateTextAreaInput(session, "pvals", value = example_sets[[input$example]])
  }, ignoreInit = TRUE)

  ## once the text no longer matches the chosen example, the menu says so, and
  ## picking the example again reloads it
  observeEvent(input$pvals, {
    ex <- input$example
    if (ex %in% names(example_sets) && trimws(input$pvals) != trimws(example_sets[[ex]]))
      updateSelectInput(session, "example", selected = "(type my own)")
  }, ignoreInit = TRUE)

  parsed <- reactive(parse_input(input$pvals))

  ## results are shown only when every entry is a valid p-value: n drives every
  ## correction, so a silently dropped entry would change all of them
  pv <- reactive({
    x <- parsed()
    req(length(x$tok) > 0, !has_decimal_comma(input$pvals),
        length(x$not_number) == 0, length(x$out_of_range) == 0)
    x$val
  })

  output$n_text <- renderUI({
    x <- parsed()
    warn <- function(...) div(style = sprintf("color:%s; font-weight:bold; margin-top:4px;", crukpink), ...)
    show <- function(v) paste(head(v, 10), collapse = ", ")
    if (!length(x$tok)) return(helpText("No p-values yet."))
    if (has_decimal_comma(input$pvals))
      return(warn("These look like decimal commas (e.g. 0,05). A comma separates p-values,
                   so 0,05 is read as 0 and 5. Write decimals with a point, and separate
                   the p-values with commas, spaces or new lines: 0.05, 0.01, 0.2"))
    probs <- c(
      if (length(x$not_number))
        sprintf("not a number: %s", show(x$not_number)),
      if (length(x$out_of_range))
        sprintf("outside [0, 1]: %s", show(x$out_of_range)))
    if (length(probs))
      return(warn(sprintf("Not valid p-values (%s). Correct or remove them to see the results.",
                          paste(probs, collapse = "; "))))
    n <- length(x$tok)
    helpText(sprintf("n = %d %s", n, plural(n, "p-value", "p-values")))
  })

  res <- reactive({
    v <- pv()
    out <- data.frame(Hypothesis = sprintf("H%d", seq_along(v)), p = v)
    for (nm in names(methods_all)) out[[nm]] <- adjust_p(v, methods_all[[nm]])
    out[order(out$p), ]
  })

  output$headline <- renderText({
    n <- length(pv())
    sprintf("%d %s, α = %.2f", n, plural(n, "hypothesis", "hypotheses"), input$alpha)
  })

  output$tab <- renderDT({
    d <- res()
    cut <- input$alpha * (1 + 1e-9)
    dt <- datatable(d, rownames = FALSE,
                    options = list(dom = table_dom(d), ordering = FALSE, pageLength = rows_per_page),
                    colnames = c("Hypothesis", "Raw p", names(methods_all))) |>
      formatSignif(c("p", names(methods_all)), 4)
    for (nm in names(methods_all)) {
      dt <- dt |> formatStyle(nm,
              fontWeight = styleInterval(cut, c("bold", "normal")),
              color      = styleInterval(cut, c(crukdark, "#BBB")))
    }
    dt
  })

  output$focus_title <- renderText({
    nm <- names(methods_all)[match(input$focus, methods_all)]
    sprintf("%s, step by step", nm)
  })

  output$focus_rule <- renderUI({
    v <- pv()
    m <- input$focus
    tagList(
      p(tags$b("Adjusted p-value: "), step_formula(m, length(v))),
      p(tags$b("Equivalent rule on the raw p-values: "),
        if (m %in% single_step)
          "every p-value is compared with the same threshold (single-step)."
        else
          "work down the ordered p-values, rejecting while each is at most its own threshold, and stop at the first failure (step-down).")
    )
  })

  output$steps <- renderDT({
    v <- pv()
    m <- input$focus; n <- length(v)
    o <- order(v)
    adj <- adjust_p(v, m)[o]
    ps  <- v[o]
    rej <- rejects(adj, input$alpha)
    dec <- ifelse(rej, "reject", "retain")
    ## step-down: say where the procedure stopped, since a later p-value can be
    ## at most its own threshold and still be retained
    stop_at <- which(!rej)[1]
    if (!(m %in% single_step) && !is.na(stop_at)) {
      dec[stop_at] <- "retain: stop here"
      if (stop_at < n) dec[(stop_at + 1):n] <- sprintf("retain (stopped at step %d)", stop_at)
    }
    d <- data.frame(
      i           = seq_len(n),
      Hypothesis  = sprintf("H%d", o),
      p_ordered   = ps,
      Threshold   = vapply(seq_len(n), function(i) step_threshold(m, n, i, input$alpha), numeric(1)),
      Calculation = vapply(seq_len(n), function(i) step_calc(m, n, i, ps[i], adj[i]), character(1)),
      Adjusted    = adj,
      Decision    = dec
    )
    datatable(d, rownames = FALSE,
              options = list(dom = table_dom(d), ordering = FALSE, pageLength = rows_per_page),
              colnames = c("Step i", "Hypothesis", "Ordered p", "Threshold",
                           "How the adjusted p-value is obtained",
                           "Adjusted p", "Decision")) |>
      formatSignif(c("p_ordered", "Threshold", "Adjusted"), 4) |>
      formatStyle("Decision",
                  color      = styleEqual("reject", crukpink, default = "#999"),
                  fontWeight = styleEqual("reject", "bold", default = "normal"))
  })

  output$plot <- renderPlot({
    d <- res()
    long <- do.call(rbind, lapply(names(methods_all), function(nm)
      data.frame(p = d$p, adj = d[[nm]], Method = nm)))
    long$Method <- factor(long$Method, levels = names(methods_all))
    caption <- "dashed grey: no adjustment    dotted pink: the level α"
    if (input$logscale) {
      dropped <- any(long$p == 0 | long$adj == 0)
      long <- long[long$p > 0 & long$adj > 0, ]
      if (dropped) caption <- paste0(caption, "\np = 0 cannot be shown on a log scale")
    }
    req(nrow(long) > 0)
    g <- ggplot(long, aes(p, adj, colour = Method, shape = Method, linetype = Method)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey75", linetype = 2) +
      geom_hline(yintercept = input$alpha, colour = crukpink, linetype = 3, linewidth = 0.9) +
      geom_point(size = 2)
    if (nrow(long) > length(methods_all)) g <- g + geom_line(linewidth = 0.7)
    g <- g +
      scale_colour_manual(values = c(crukdark, cruklight, "#5B8C00", "#7B2D8E")) +
      scale_shape_manual(values = c(16, 17, 15, 4)) +
      scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed")) +
      labs(x = "raw p-value", y = "adjusted p-value", caption = caption) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right", legend.title = element_blank())
    if (input$logscale) g + scale_x_log10() + scale_y_log10()
    else g + coord_cartesian(ylim = c(0, 1))
  })

  output$count_text <- renderUI({
    d <- res()
    k <- vapply(names(methods_all), function(nm) sum(rejects(d[[nm]], input$alpha)), numeric(1))
    raw <- sum(rejects(d$p, input$alpha))
    tagList(
      p(tags$b(sprintf("Uncorrected: %d of %d rejected.", raw, nrow(d)))),
      tags$ul(lapply(names(methods_all), function(nm)
        tags$li(HTML(sprintf("%s: <b>%d</b>", nm, k[[nm]])))))
    )
  })

  ## ---- saving ------------------------------------------------------------
  ## the CSV carries alpha and each method's decisions, and starts with a UTF-8
  ## byte-order mark so that Excel shows "Šidák" correctly
  output$save <- downloadHandler(
    filename = function() sprintf("adjusted-p-%s.csv", format(Sys.time(), "%Y%m%d-%H%M%S")),
    content  = function(file) {
      d <- tryCatch(res(), error = function(e) NULL)
      out <- if (is.null(d)) {
        data.frame(Note = "No valid p-values were entered.")
      } else {
        o <- data.frame(Hypothesis = d$Hypothesis, `Raw p` = d$p, alpha = input$alpha,
                        check.names = FALSE)
        for (nm in names(methods_all)) {
          o[[paste(nm, "adjusted p")]] <- d[[nm]]
          o[[paste(nm, "decision")]]   <- ifelse(rejects(d[[nm]], input$alpha), "reject", "retain")
        }
        o
      }
      con <- file(file, open = "w", encoding = "UTF-8")
      on.exit(close(con))
      writeLines("﻿", con, sep = "")
      write.csv(out, con, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
