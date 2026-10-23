# Shiny apps for the Multiple Comparisons course

Four apps, meant to be used in order. They follow the course from the
arithmetic of adjusted p-values (app 1), to the FWER and power of the
general-purpose corrections (app 2), to the corrections built for comparing
group means (app 3), and finally to the FDR (app 4).

| App | File | What it is for |
|---|---|---|
| 1 | `1.Adjusted_p_calculus.R` | Calculating adjusted p-values by hand: Bonferroni, Holm, Šidák, Holm–Šidák |
| 2 | `2.Adjusted_p_FWER_power.R` | FWER and power of the same four corrections, by simulation |
| 3 | `3.Adjusted_p_Dunnett_Tukey_p_FWER_power.R` | Dunnett and Tukey (single-step and step-down), with alternatives for unequal variances and non-normal data |
| 4 | `4.FDR_simulation.R` | FWER against FDR in a gene expression experiment: Bonferroni, Holm, Benjamini–Hochberg, Benjamini–Yekutieli |

## Running them

Each app is a single file. None of them is called `app.R`, so `runApp()` has
to be given the file, not the folder:

```r
setwd("course-slides/Final_version/Shiny_apps")
shiny::runApp("1.Adjusted_p_calculus.R")
shiny::runApp("2.Adjusted_p_FWER_power.R")
shiny::runApp("3.Adjusted_p_Dunnett_Tukey_p_FWER_power.R")
shiny::runApp("4.FDR_simulation.R")
```

All four apps need `shiny`, `shinyWidgets`, `DT` and `ggplot2`. App 3 also
needs `mvtnorm`.

## Common layout (apps 2–4)

- **Left**: the parameters, and the **Run simulation** button.
- **Top**: a **Free exploration** tab, followed by five guided steps. Each step
  applies a set of parameters and tells you what to look at. Its explanation
  stays hidden behind **Show the explanation** until you have run the
  simulation and read the numbers for yourself.
- **Right**: the summary table, a plot of the error rate against power, and the
  explanation of the current step.

Across the simulated experiments, the apps report:

- **FWER**: the proportion of experiments with at least one false positive.
- **Power**: the mean proportion of real effects that are detected.
- **FDR** (app 4 only): the mean of `Q = V/R`, with `Q = 0` when `R = 0`. This
  is the Benjamini–Hochberg definition, not the positive FDR.

---

## 1. Adjusted p-values: the calculation

This app does no simulation. You type in a list of raw p-values, or pick one of
the examples (two hypotheses, ten genes, twenty p-values),
and choose α. The app then shows:

- a table of the adjusted p-values for **Bonferroni, Holm, Šidák and
  Holm–Šidák**. Hypotheses rejected at α are shown in bold.
- the arithmetic for a method you choose. Each step shows the ordered p-value,
  the threshold, the multiplication by `n` or `n − i + 1` (or `1 − (1 − p)^k`
  for Šidák), and, for the step-down methods, the running maximum that keeps
  the adjusted p-values monotone.
- a plot of raw against adjusted p-values, with the number of rejections for
  each method.

The app covers the FWER only. The FDR is in app 4.

## 2. FWER and power of the four corrections

A family of `n` tests, of which a chosen percentage are true nulls. No data are
simulated. Each test's two-sample t statistic, on `2(n_arm − 1)` degrees of
freedom, is drawn directly from its sampling distribution, with the same
one-factor dependence model as app 4 (see below). The parameters are:

- **The family**: `n` and the percentage of true nulls.
- **Dependence**: `|ρ|` (0 = independent), and the proportion of tests loading
  the opposite way.
- **Experiment**: effect size, samples per group, one- or two-sided test, α,
  and the number of simulated experiments.

The app compares no correction, Bonferroni, Holm, Šidák and Holm–Šidák.

Guided steps:

1. **The error rate grows with the family.** The uncorrected FWER follows
   `1 − (1 − α)^m0`.
2. **Holm beats Bonferroni, for free.** Holm has the same FWER guarantee and is
   never less powerful.
3. **What Šidák adds, and what it assumes.** The gain over Bonferroni is tiny,
   and it is justified by assuming the tests are independent.
4. **When the tests are not independent.** Positive dependence makes every
   correction conservative, and costs power.
5. **Many tests, few real effects.** This is where FWER control becomes costly,
   and it leads into app 4.

## 3. Dunnett, Tukey and their alternatives

Real data are simulated: `k` groups of `n_arm` observations, of which a chosen
number of groups carry a real effect. You can choose between two families of
comparisons:

- **Against a control** (Dunnett family): `k − 1` comparisons.
- **All pairs** (Tukey family): `k(k − 1)/2` comparisons.

Two settings change the data:

- **Largest / smallest SD**: the group SDs rise linearly from group 1 (the
  control) to group `k`. A ratio of 1 means equal variances.
- **Right-skewed data**: standardised exponential errors instead of normal
  ones.

| | Against a control | All pairs |
|---|---|---|
| General-purpose | Bonferroni, Holm, Šidák, Holm–Šidák | Bonferroni, Holm, Šidák, Holm–Šidák |
| Max-t, single-step | Dunnett | Tukey |
| Max-t, step-down | Holm–Dunnett | Holm–Tukey |
| Unequal variances | Welch + Holm (in the spirit of T2/T3) | Games–Howell |
| Non-normal data (ranks) | Steel | Steel–Dwass |

How the methods work:

- The pooled-variance methods use a common SD with `k(n_arm − 1)` degrees of
  freedom.
- The critical values are computed once for each configuration, not once per
  experiment. This keeps the app fast (300 experiments take under a second).
- **Dunnett** uses the multivariate t quantile, with correlation 1/2 between
  comparisons.
- **Tukey** uses the studentized range, `q(α, k, ν)/√2`.
- **Steel and Steel–Dwass** use the asymptotic normal distribution of the
  pairwise rank sums.

Because every group has the same size, these critical values are exact for
Dunnett, Tukey and the step-down Dunnett. **Holm–Tukey** uses the usual
approximation: at each step, the range is taken over the number of groups
still involved.

Under the global null with equal variances (20,000 experiments, `k = 4`,
`n_arm = 10`), the simulated FWER is 0.050 for Dunnett (theory: 0.0500), 0.049
for Tukey and 0.049 for Games–Howell. Bonferroni, Holm and Šidák are
conservative (0.039–0.044).

Guided steps:

1. **Using the design beats correcting blindly.**
2. **The power that correlation buys.**
3. **When the variances differ.**
4. **When the data are not normal.**
5. **Many groups.**

## 4. FWER and FDR in a gene expression experiment

One experiment consists of `N` genes, measured in two groups of `n_arm`
samples each. `m0` of the genes are true nulls. The remaining `m1 = N − m0` are
differentially expressed, with the chosen effect size.

No expression values are simulated. Each gene's two-sample t statistic is
drawn directly from its sampling distribution. Its numerator is

```
z_i = mu_i + s_i * sqrt(|rho|) * F + sqrt(1 - |rho|) * e_i
```

and it is divided by that gene's own variance estimate:
`t_i = z_i / sqrt(w_i)`, where `w_i ~ chi-squared(df)/df` and
`df = 2(n_arm − 1)`.

- The single latent factor `F` stands for whatever makes genes move together: a
  batch effect, a shared pathway, or a common control sample.
- `mu_i = 0` for a true null, and `effect * sqrt(n_arm/2)` otherwise.
- Because each gene carries its own variance estimate, the correlation between
  the t statistics is slightly below `|rho|`.

The dependence between genes can be set to one of three cases:

- **Independent**: `rho = 0`, so the genes share no common signal.
- **Positive (co-regulated)**: all loadings are `s_i = +1`, so every pair of
  genes has correlation `+|rho|`.
- **Negative (opposite directions)**: a randomly chosen share of the genes
  load `−1`. The two blocks then push against each other, with correlation
  `+|rho|` within a block and `−|rho|` between blocks. A two-sided test sees
  only `|t|`, so flipping the sign of a null gene's loading leaves the joint
  distribution of the null p-values exactly unchanged, and the negative
  structure would be invisible. The negative-dependence step therefore switches
  to a one-sided test, where the sign does reach the p-values.

You choose what to control: the FWER at level α (Bonferroni, Holm), or the FDR
at level q (BH, BY).

Guided steps:

1. **All genes null.** Under the global null, FDR = FWER exactly.
2. **A few real effects.** Bonferroni's power collapses while BH holds up, and
   BH's realised FDR sits near `π0 · q`, not `q`.
3. **More genes.** The uncorrected FWER saturates at 1 while the false
   positives pile up. The FDR does not move.
4. **Co-regulated genes.** The mean FDP is still controlled, but its spread
   grows. The FDR is a promise about the average experiment, not about yours.
5. **Genes pushing in opposite directions.** This is the case BH's proof does
   not cover, and it shows what BY costs in power.

App 4 also shows the experiment-by-experiment distribution of the false
discovery proportion, and the p-values from a single experiment.

---

## Saving results

Each app has a **Download** button, which sends the current table to your
browser as a timestamped CSV.

The simulation apps include the settings in every row of the CSV, so a saved
file can be read on its own.
