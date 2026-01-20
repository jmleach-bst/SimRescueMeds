
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `SimRescueMeds`

<!-- badges: start -->

<!-- badges: end -->

The goal of `SimRescueMeds` is to enable (relatively) simple simulation
studies that impose rescue medication (or really any intercurrent event)
on longitudinal data.

1.  Longitudinal data is simulated according a linear mixed model. This
    is the data one would observe if no intercurrent event, e.g., rescue
    therapy, occurs.
2.  Rescue medication/intercurrent events are simulated using models
    adapted from Thomadakis et al. (2019), which allows for reasonably
    straightforward generation of events that, if excluding post-event
    data, produce MCAR, MAR, or MNAR missingness on the longitudinal
    data using either exponential or Weibull models.
3.  We also intend to develop functions that can modify post-event data.
    Only one of these has been developed thus far,
    `rescue_effect_means()`; see documentation for details.

## Installation

You can install the development version of `SimRescueMeds` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jmleach-bst/SimRescueMeds")
```

## Simulation examples

First, load the package.

``` r
library(SimRescueMeds)
#> Loading required package: matrixcalc
#> Loading required package: MASS
```

Use `build_design_matrix()` to build the design matrix for the fixed
effects, in which you specify the number of individuals in each group
`N_t` and `N_c`, the number of time points `K`, the starting time, and
time scale. The example below scales such that the file time point has
$t = 1$, e.g., 1-year follow-up.

``` r
xmat_t <- build_design_matrix(
    N_t = 1,
    N_c = 0,
    K = 5,
    time_start = 0,
    time_scale = 1/4
  )
xmat_t
#>      intercept x time x_by_time
#> [1,]         1 1 0.00      0.00
#> [2,]         1 1 0.25      0.25
#> [3,]         1 1 0.50      0.50
#> [4,]         1 1 0.75      0.75
#> [5,]         1 1 1.00      1.00
```

Use `simulate_lmm_rct()` to generate data sets based on the linear mixed
model, which generates data from a multivariate Normal distribution.
Here, we specify a single-individual design matrix for each of 2
treatment groups and the number of individuals in each group. Use
`build_lmm_cov()` to more easily construct the covariance matrix as a
function of random effects (co)variances and within-individual random
error covariances.

``` r
set.seed(3002242)
df <- simulate_lmm_rct(
  N_t = 50,
  N_c = 50,
  xmat_t = xmat_t,
  xmat_c = build_design_matrix(
    N_t = 0,
    N_c = 1,
    K = 5,
    time_start = 0,
    time_scale = 1/4
  ),
  betas = c(
    beta0 = 8,
    beta1 = 0,
    beta2 = -1,
    beta3 = -1
  ),
  Sigma = build_lmm_cov(
    zmat = xmat_t[, c("intercept", "time")],
    re_sigma = c(1.15, 1.05),
    re_corr_mat = matrix(
      c(1, 0,
        0, 1),
      byrow = TRUE,
      nrow = 2
      ),
    ws_sigma = 1.25,
    ws_corr_mat = diag(5),
    print_intermediate = FALSE
    )
)

head(df, 10)
#>           y x time x_by_time id
#> 1  6.777270 1 0.00      0.00  1
#> 2  6.503446 1 0.25      0.25  1
#> 3  6.408070 1 0.50      0.50  1
#> 4  5.286074 1 0.75      0.75  1
#> 5  6.041579 1 1.00      1.00  1
#> 6  7.642811 1 0.00      0.00  2
#> 7  7.122279 1 0.25      0.25  2
#> 8  7.278219 1 0.50      0.50  2
#> 9  5.342140 1 0.75      0.75  2
#> 10 6.014225 1 1.00      1.00  2
```

Simulate rescue therapy based on a proportional hazards model in the
manner of Thomadakis et al. (2019), which allows the hazard of rescue
therapy within some interval $t_{i,j} \le t < t_{i,j+1}$ to be a
function of the last observed outcome measurement, $y_{i,j}$, and the
next measurement, $y_{i,j+1}$, the latter of which is unobserved in the
study since it occurs after rescue therapy. While simple in some
respects, this approach allows us to easily define MCAR
$(\theta_1 = \theta_2 = 0)$, MAR $(\theta_1 \ne 0,\,\theta_2 = 0)$, or
MNAR $(\theta_2 \ne 0)$.

$$
h_i(t) = \lambda\exp\left(y_{i,j}\theta_1 + y_{i,j+1}\theta_2\right), \quad t_{i,j} \le t < t_{i,j+1}.
$$

``` r
set.seed(707582)

rt <- sim_hazard_thomadakis_df(
  data = df,
  theta1 = 1/3,
  theta2 = 0,
  lambda = 0.05,
  dist = "exponential",
  print_intermediate = FALSE
)

head(rt)
#>   id        tte  eventtime event event_factor
#> 1  1 0.17237285 0.17237285     1            1
#> 2  2 5.54637407 1.00000000     0            0
#> 3  3 8.41025687 1.00000000     0            0
#> 4  4 0.03216887 0.03216887     1            1
#> 5  5 2.42347970 1.00000000     0            0
#> 6  6 0.52273847 0.52273847     1            1
```

Finally, for analysis we would exclude the events that occur post-event
using `exclude_post_tte()`.

``` r
post_rt <- exclude_post_tte(
  original_data = df,
  tte_data = rt
)

head(post_rt, 15)
#>    id        y x time x_by_time        tte  eventtime event event_factor
#> 1   1 6.777270 1 0.00      0.00 0.17237285 0.17237285     1            1
#> 2   2 7.642811 1 0.00      0.00 5.54637407 1.00000000     0            0
#> 3   2 7.122279 1 0.25      0.25 5.54637407 1.00000000     0            0
#> 4   2 7.278219 1 0.50      0.50 5.54637407 1.00000000     0            0
#> 5   2 5.342140 1 0.75      0.75 5.54637407 1.00000000     0            0
#> 6   2 6.014225 1 1.00      1.00 5.54637407 1.00000000     0            0
#> 7   3 6.083945 1 0.00      0.00 8.41025687 1.00000000     0            0
#> 8   3 7.900913 1 0.25      0.25 8.41025687 1.00000000     0            0
#> 9   3 7.085961 1 0.50      0.50 8.41025687 1.00000000     0            0
#> 10  3 4.586002 1 0.75      0.75 8.41025687 1.00000000     0            0
#> 11  3 2.852984 1 1.00      1.00 8.41025687 1.00000000     0            0
#> 12  4 8.579877 1 0.00      0.00 0.03216887 0.03216887     1            1
#> 13  5 6.900070 1 0.00      0.00 2.42347970 1.00000000     0            0
#> 14  5 7.325157 1 0.25      0.25 2.42347970 1.00000000     0            0
#> 15  5 7.661118 1 0.50      0.50 2.42347970 1.00000000     0            0
#>    last_time_preRT last_y_preRT
#> 1                0     6.777270
#> 2                1     6.014225
#> 3                1     6.014225
#> 4                1     6.014225
#> 5                1     6.014225
#> 6                1     6.014225
#> 7                1     2.852984
#> 8                1     2.852984
#> 9                1     2.852984
#> 10               1     2.852984
#> 11               1     2.852984
#> 12               0     8.579877
#> 13               1     7.294891
#> 14               1     7.294891
#> 15               1     7.294891
```

## Calculation of Monte Carlo Standard Errors (MCSE)

It is good practice to report some measure of precision associated with
the estimates from a simulation study, because in general the results
can change if a new seed, and therefore new data, are generated. We
largely follow Morris, White, and Crowther (2019) in using formulas for
these calculations, but we also include options for jackknife estimation
of MCSE, as described by Koehler, Brown, and Haneuse (2009) and
re-emphasized by Kelter (2024). These are general functions that are not
necessarily dependent on simulating interrcurrent events.

Consider a very simple example, in which we estimate the regression
parameter $\beta$ from a simple linear model:

$$
Y_i = \alpha + X_i\beta + \epsilon_i 
$$

where $X_i$ is a binary predictor taking value 0 or 1 and
$\epsilon_i \sim \mathcal{N}(0, \sigma^2$. This is straightforward and
quick to simulate. Note that we would usually save both the data and the
analysis, but for this demonstration it is unnecessary. Below we use
`mcse_estimands()`, which reports several estimates and MCSE for some
typical estimands of interest. In this case the data frame input only
consists of a column for estimates of $\hat{\beta}$ and a column for
estimates of $\widehat{SE}(\hat{\beta})$. We have arguments for
`estimand_name` and `se_name`, which tells the function which columns
are which, and we need to specify `true_value` in order to calculate,
e.g., bias. We see that for this relatively simple model, even
$M = 2500$ produces relatively small MCSE for every measure except
percent-bias. In general, we may except the MCSE to decrease with
increasing $M$.

``` r
# Set seed
set.seed(723862)

# Number of simulations
M <- 2500

# Sample size
N <- 150

# Parameters
alpha <- 0.5
beta <- 1.5
sigma <- 1.25

# Balanced groups
x <- c(
  rep(1, N/2),
  rep(0, N/2)
)

# Simulate data
beta_hat <- c()
se_hat <- c()
for (m in 1:M) {
  dfm <- data.frame(
    x = x,
    y = alpha + x*beta + rnorm(N, 0, sigma)
  )
  lm_m <- lm(
    formula = y ~ x,
    data = dfm 
  )
  beta_hat[m] <- coefficients(lm_m)["x"]
  se_hat[m] <- summary(lm_m)$coefficients["x", "Std. Error"]
}
df <- data.frame(beta_hat, se_hat)

knitr::kable(
  mcse_estimands(
    data = df,
    estimand_name = "beta_hat",
    se_name = "se_hat",
    include_bias_percent = TRUE,
    report_proportion = FALSE,
    true_value = beta
  )
)
```

| measure      |   estimate |      mcse |
|:-------------|-----------:|----------:|
| Percent-Bias | -0.0242720 | 0.2715892 |
| Bias         | -0.0003641 | 0.0040738 |
| EmpSE        |  0.2036919 | 0.0028812 |
| ModSE        |  0.2040271 | 0.0002385 |
| MSE          |  0.0414739 | 0.0012012 |

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Kelter2024" class="csl-entry">

Kelter, Riko. 2024. “The Bayesian Simulation Study (BASIS) Framework for
Simulation Studies in Statistical and Methodological Research.”
*Biometrical Journal* 66 (January).
<https://doi.org/10.1002/bimj.202200095>.

</div>

<div id="ref-Koehler2009" class="csl-entry">

Koehler, Elizabeth, Elizabeth Brown, and Sebastien J. P. A. Haneuse.
2009. “On the Assessment of Monte Carlo Error in Simulation-Based
Statistical Analyses.” *American Statistician* 63: 155–62.
<https://doi.org/10.1198/tast.2009.0030>.

</div>

<div id="ref-Morris2019" class="csl-entry">

Morris, Tim P., Ian R. White, and Michael J. Crowther. 2019. “Using
Simulation Studies to Evaluate Statistical Methods.” *Statistics in
Medicine* 38: 2074–2102. <https://doi.org/10.1002/sim.8086>.

</div>

<div id="ref-Thomadakis2019" class="csl-entry">

Thomadakis, Christos, Loukia Meligkotsidou, Nikos Pantazis, and Giota
Touloumi. 2019. “Longitudinal and Time-to-Drop-Out Joint Models Can Lead
to Seriously Biased Estimates When the Drop-Out Mechanism Is at Random.”
*Biometrics* 75 (March): 58–68. <https://doi.org/10.1111/biom.12986>.

</div>

</div>
