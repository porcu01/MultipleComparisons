library(shiny)
library(DT)
library(ggplot2)

crukdark  <- "#00007E"
cruklight <- "#009CEE"
crukpink  <- "#FF0087"

## the same colour and shape for each correction as in app 2
method_cols   <- c("Bonferroni" = crukdark, "Holm" = cruklight,
                   "Šidák" = "#5B8C00", "Holm–Šidák" = "#7B2D8E")
method_shapes <- c("Bonferroni" = 16, "Holm" = 17, "Šidák" = 15, "Holm–Šidák" = 18)

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
## k may be a single number (single-step) or one value per hypothesis (step-down),
## so recycle it before ifelse(), which otherwise returns only length(k) values
sidak_adj <- function(p, k) {
  k <- rep_len(k, length(p))
  ifelse(k == 1, p, -expm1(k * log1p(-p)))
}

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

## how each method builds its adjusted p-value, shown step by step (HTML)
step_formula <- function(method, n) {
  switch(method,
    bonferroni = sprintf("p̃<sub>i</sub> = min(1, %d × p<sub>i</sub>)", n),
    holm       = sprintf("p̃<sub>(i)</sub> = min(1, max<sub>j ≤ i</sub> (%d − j + 1) × p<sub>(j)</sub>)", n),
    sidak      = sprintf("p̃<sub>i</sub> = 1 − (1 − p<sub>i</sub>)<sup>%d</sup>", n),
    holm_sidak = sprintf("p̃<sub>(i)</sub> = max<sub>j ≤ i</sub> [1 − (1 − p<sub>(j)</sub>)<sup>%d − j + 1</sup>]", n)
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

## the tables and the step-by-step arithmetic all show 6 significant figures, so
## that an adjusted p-value and alpha (or a p-value and its threshold) never
## collapse onto the same rounded number when the decisions differ
sig_digits <- 6

## the arithmetic of a single step, written out with the actual numbers:
## the true product first, then any cap at 1, then any running maximum
step_calc <- function(method, n, i, p, adj, digits = 6) {
  fmt <- function(x) sprintf("%.*g", digits, x)
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
  ## "0,1" and "0,0" are legitimate p-value pairs, so only flag digits after the
  ## comma that cannot stand on their own as a p-value (0,05 and 0,5 but not 0,1)
  grepl("(^|[^0-9.eE+-])0,(0[0-9]|[2-9])[0-9]*($|[^0-9.])", txt)
}

## "0,1" and "0,0" are accepted as two p-values, but the reader is told so
has_ambiguous_comma <- function(txt)
  grepl("(^|[^0-9.eE+-])0,[01]($|[^0-9.])", txt)

plural <- function(k, one, many) if (k == 1) one else many

## both tables page at ten rows, so a family of any size stays on one screen;
## beyond that DataTables adds the row count and the page buttons
rows_per_page  <- 10
steps_per_page <- 10
table_dom <- function(d, per_page = rows_per_page)
  if (nrow(d) > per_page) "tip" else "t"

example_sets <- list(
  "Two hypotheses"                           = "0.0000000063, 0.048",
  ## at alpha = 0.05 the four corrections reject 2 (Bonferroni), 3 (Sidak),
  ## 5 (Holm) and 6 (Holm-Sidak, as many as without correction)
  "Ten tests"                                = "0.001, 0.003, 0.0051, 0.007, 0.008, 0.0101, 0.130, 0.340, 0.560, 0.870",
  ## chosen so that all four corrections separate at thresholds the slider reaches:
  ## Holm beats Bonferroni at 0.07, Holm-Sidak beats both at 0.18, Sidak beats
  ## Bonferroni at 0.22
  "Twenty p-values"                          = paste(sprintf("%.3f", sort(c(0.001, 0.002, 0.003, 0.004, 0.012,
                                                                           0.031, 0.049, 0.067, 0.091, 0.120,
                                                                           0.180, 0.210, 0.260, 0.330, 0.370,
                                                                           0.440, 0.510, 0.580, 0.630, 0.700))),
                                                     collapse = ", ")
)

## ---------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML(sprintf("
    .well { background:#F7FAFF; border-color:%s33; }
    h4 { color:%s; }
    .note { border-left:5px solid %s; padding:6px 12px; background:#F2F9FF; }
    .formula { font-size:17px; }
    .intro p { margin-bottom:8px; }
  ", crukdark, crukdark, cruklight)))),

  titlePanel(div(style = sprintf("color:%s", crukdark),
                 "From raw to adjusted p-values: Bonferroni, Holm, Šidák, Holm–Šidák"),
             windowTitle = "From raw to adjusted p-values"),

  div(class = "intro", style = "margin:14px 0 18px 0; color:#444; font-size:15px;",
      HTML(sprintf("
        <p>Enter the raw p-values of a <b>family of hypotheses</b>, the tests considered together. For each test the <b>adjusted</b>
        p-value under each correction is reported, and the arithmetic that produces it
        is shown step by step for the correction you choose. A hypothesis is rejected
        when its adjusted p-value is at most the level α, so each correction's
        decisions can be read from its own column.</p>
        <p>Each correction controls the <b style='color:%s'>family-wise error rate
        (FWER)</b>, the probability of <i>at least one</i> false positive anywhere in the
        family of hypotheses, at level α, under its own condition. <b>Bonferroni</b> and <b>Holm</b>:
        whatever the dependence among the p-values. <b>Šidák</b> and <b>Holm–Šidák</b>:
        when the p-values of the true null hypotheses are independent or, more
        generally, positively upper orthant dependent (PUOD), a sufficient condition:
        for any threshold c (e.g. 0.05), these p-values are at least as likely to be
        all non-significant together as they would be if they were independent.</p>
        <p>Holm rejects at least as much as Bonferroni, Šidák at least as much as
        Bonferroni, and Holm–Šidák at least as much as both Holm and Šidák. Šidák and
        Holm, however, are not ordered: depending on the p-values, either can reject
        more.</p>",
        crukdark))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("The family of hypotheses"),
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
      h4(textOutput("focus_title")),
      div(class = "note", uiOutput("focus_rule")),
      br(),
      DTOutput("steps"),
      br(),
      fluidRow(
        column(8,
               h4("Adjusted against raw p-values"),
               plotOutput("plot", height = "380px"),
               checkboxInput("logscale", "Log scale (spreads out the small p-values)", value = TRUE)),
        column(4,
               h4("How many are rejected"),
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
    tagList(
      helpText(sprintf("n = %d %s", n, plural(n, "p-value", "p-values"))),
      if (has_ambiguous_comma(input$pvals))
        div(style = sprintf("color:%s; margin-top:4px;", crukpink),
            "Note: 0,1 (or 0,0) is read as two p-values, 0 and 1. If you meant 0.1, write it with a point."))
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
    ## the decision is taken once, in R, by the same rejects() the step table,
    ## the counts and the CSV use. Deriving it again in the browser would drop
    ## the 1e-9 tolerance, because DT rounds the cut-off to 4 decimals on its
    ## way into the JavaScript, and a value sitting exactly on alpha would then
    ## be bold in one table and grey in the other.
    flags <- vapply(names(methods_all),
                    function(nm) paste0(".", nm), character(1))
    for (nm in names(methods_all))
      d[[paste0(".", nm)]] <- ifelse(rejects(d[[nm]], input$alpha), "y", "n")
    dt <- datatable(d, rownames = FALSE,
                    options = list(dom = table_dom(d), ordering = FALSE,
                                   pageLength = rows_per_page,
                                   columnDefs = list(list(visible = FALSE,
                                                          targets = unname(flags)))),
                    colnames = c("Hypothesis", "Raw p", names(methods_all),
                                 unname(flags))) |>
      formatSignif(c("p", names(methods_all)), sig_digits)
    for (nm in names(methods_all)) {
      dt <- dt |> formatStyle(nm, valueColumns = paste0(".", nm),
              fontWeight = styleEqual("y", "bold", default = "normal"),
              color      = styleEqual("y", crukdark, default = "#767676"))
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
      p(tags$b("Adjusted p-value: "), span(class = "formula", HTML(step_formula(m, length(v))))),
      if (length(v) == 1)
        p(tags$b("With a single hypothesis there is nothing to correct:"),
          " every adjusted p-value equals the raw one."),
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
    thr <- vapply(seq_len(n), function(i) step_threshold(m, n, i, input$alpha), numeric(1))
    if (!(m %in% single_step) && !is.na(stop_at)) {
      ## after the stop no comparison is made, so no threshold is shown
      if (stop_at < n) thr[(stop_at + 1):n] <- NA
      dec[stop_at] <- "retain: stop here"
      if (stop_at < n) dec[(stop_at + 1):n] <- sprintf("retain (stopped at step %d)", stop_at)
    }
    d <- data.frame(
      i           = seq_len(n),
      Hypothesis  = sprintf("H%d", o),
      p_ordered   = ps,
      ## formatted here, so that a step with no comparison shows a dash
      Threshold   = ifelse(is.na(thr), "—", sprintf("%.*g", sig_digits, thr)),
      Calculation = vapply(seq_len(n), function(i) step_calc(m, n, i, ps[i], adj[i], sig_digits),
                           character(1)),
      Adjusted    = adj,
      Decision    = dec
    )
    datatable(d, rownames = FALSE,
              options = list(dom = table_dom(d, steps_per_page), ordering = FALSE,
                             pageLength = steps_per_page),
              colnames = c("Step i", "Hypothesis", "Ordered p", "Threshold",
                           "How the adjusted p-value is obtained",
                           "Adjusted p", "Decision")) |>
      formatSignif(c("p_ordered", "Adjusted"), sig_digits) |>
      formatStyle("Threshold", textAlign = "right") |>
      formatStyle("Decision",
                  color      = styleEqual("reject", crukdark, default = "#767676"),
                  fontWeight = styleEqual("reject", "bold", default = "normal"))
  })

  output$plot <- renderPlot({
    d <- res()
    long <- do.call(rbind, lapply(names(methods_all), function(nm)
      data.frame(p = d$p, adj = d[[nm]], Method = nm)))
    long$Method <- factor(long$Method, levels = names(methods_all))
    caption <- "grey dashed line: no adjustment (adjusted = raw)\npink dotted line: the level α"
    if (input$logscale) {
      dropped <- any(long$p == 0 | long$adj == 0)
      long <- long[long$p > 0 & long$adj > 0, ]
      if (dropped) caption <- paste0(caption, "\np = 0 cannot be shown on a log scale")
    }
    validate(need(nrow(long) > 0,
                  "Every p-value is 0, which cannot be shown on a log scale: untick Log scale."))
    g <- ggplot(long, aes(p, adj, colour = Method, shape = Method, linetype = Method)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey75", linetype = "longdash") +
      geom_hline(yintercept = input$alpha, colour = crukpink, linetype = 3, linewidth = 0.9) +
      geom_point(size = 2)
    if (nrow(long) > length(methods_all)) g <- g + geom_line(linewidth = 0.7)
    g <- g +
      scale_colour_manual(values = method_cols) +
      scale_shape_manual(values = method_shapes) +
      scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
      labs(x = "raw p-value", y = "adjusted p-value", caption = caption) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            plot.caption = element_text(hjust = 0, colour = "grey35"))
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
        tags$li(HTML(sprintf("%s: <b>%d</b> of %d rejected", nm, k[[nm]], nrow(d))))))
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
