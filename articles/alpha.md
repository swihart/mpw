# Weibull-Weibull-Positive Stable

## PS($\alpha,\alpha,0$) Positive Stable Frailty Example

This is the $F_{\alpha}$ referenced in the manuscript:

- $F_{\alpha}(x) = 1 - \exp\left\lbrack - \exp\alpha\eta \right\rbrack$

This is the PS($\alpha,\alpha,0$) Positive Stable Frailty density:

- $f_{u_{i}}\left( u_{i}|\alpha \right) = - \frac{1}{\pi u_{i}}\sum_{k = 1}^{\infty}\frac{\Gamma(k\alpha + 1)}{k!}\left( - u_{i}^{- \alpha} \right)^{k}\sin(\alpha k\pi)$

We analyze data that looks similar to that of Figure 2 [Thomas et al
(2021) *Safety and Efficacy of the BNT162b2 mRNA Covid-19 Vaccine
through 6 Months*](https://pubmed.ncbi.nlm.nih.gov/34525277/).

We plot the data points, color-code the groups, and just connect the
dots with interpolated, dashed lines.

``` r
library(mpw)
## This data resembles Figure 2 in Thomas et al
## https://pubmed.ncbi.nlm.nih.gov/34525277/
time<- c(0,14,28,42,56,70,84,98,112,126,140,154,168,182,196)
F1  <- c(0,.18,.19,.22,.25,.27,.28,.34,.44,.50,.60,.72,.75,.81,.93)/100
F0  <- c(0,.29,.60,1,1.38,1.75,2.25,2.97,3.50,4.25,4.94,5.53,6.00,6.31,6.94)/100

## plot the data points, 
## with interpolated lines
plot  (time, F0, type="l", col="black", xlab="Days", ylab="CDF", lty=2)
points(time, F0, col="black", pch=15)
lines (time, F1, col="blue", lty=2)
points(time, F1, col="blue", pch=16)
legend("topleft", c("Placebo", "Vaccine"),
       col=c("black","blue"),
       lty=1,
       pch=c(15,16)
)
```

![](alpha_files/figure-html/thomas-data-1.png)

Set the `h_parm` and `frailty` distribution for this example. Beware
that H_PARM has different bounds depending on what value of FRAILTY is
selected.

``` r
H_PARM <- 0.72
FRAILTY <- "PS"
```

We first consider two knots which will allow us to model each curve with
three pieces. We set the first piece of each group to be equal so we can
get an HR of 1 (more in subsequent sections). This first piece can be
specified to be over a short time range where both groups would be
considered equal due to the *ramp-up time* of a vaccine. Below we set
the ramp-up time to end 5 days after first dose, the second piece to go
from then to 90 days after second dose, and the third piece from then
until the end of the study. We do this by specifying two knot values, as
below:

### Two knots

``` r
tvec.in <- c(5, 90+21)
tvec.in
#> [1]   5 111
```

#### Two step fitting procedure

##### 1. Fit for placebo group

``` r
fit_F_alpha <- function(x){mean((popavg_dist(time, 
                                            knots= tvec.in, 
                                            logk0=x[1], 
                                            g0=x[2], 
                                            delta_vec=x[-c(1,2)],
                                            h_parm = H_PARM,
                                            frailty=FRAILTY)  - F0)^2) }

init.vals <- c(log(2.4), -10, rep( 0,length(tvec.in)))
plac_fit_F_alpha <- 
  optim(init.vals, fit_F_alpha,
        method="Nelder-Mead",
        control=list(maxit=1e8))
print(plac_fit_F_alpha)
#> $par
#> [1]  -0.8590240 -11.5514302   1.6053115  -0.3700213
#> 
#> $value
#> [1] 1.574351e-06
#> 
#> $counts
#> function gradient 
#>      167       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL

logk0p = plac_fit_F_alpha$par[1]
g0p = plac_fit_F_alpha$par[2]
delta_vec_p =plac_fit_F_alpha$par[-1*c(1,2)]
```

##### 2. Fit for vaccine group

Now take `logk0p` and `g0p`, the parameters that control the first
(leftmost) piece for the placebo group and use them for `logk0v` and
`g0v` when estimating the pieces for the vaccine group. This forces the
first pieces to be the same for both groups (i.e. *ramp-up time*).

``` r
fit_F_alpha <- function(x){mean((popavg_dist(time, 
                                            knots= tvec.in, 
                                            logk0=logk0p, 
                                            g0=g0p, 
                                            delta_vec=x,
                                            h_parm = H_PARM,
                                            frailty=FRAILTY)  - F1)^2) }

init.vals <- rep( 0,length(tvec.in))
vacc_fit_F_alpha <- 
  optim(init.vals, fit_F_alpha,
        method="Nelder-Mead",
        control=list(maxit=1e8))
print(vacc_fit_F_alpha)
#> $par
#> [1] 0.6274834 0.9321016
#> 
#> $value
#> [1] 1.41648e-07
#> 
#> $counts
#> function gradient 
#>       75       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL

logk0v = logk0p
g0v = g0p
delta_vec_v =vacc_fit_F_alpha$par
```

#### Plot population avg CDF

- Quick check of fit

We use the knots and the original time in `time.dist` to do a quick
check of the fit. We put vertical dashed lines in the plot to denote
knot placement.

``` r
time.dist <- sort(c(seq(min(c(time, tvec.in)), 
                        max(c(time, tvec.in)), 0.1),
                    tvec.in)
                  )

plac.dist <- popavg_dist(      x = time.dist, 
                           knots = tvec.in,
                           logk0 = logk0p,
                              g0 = g0p,
                       delta_vec = delta_vec_p,
                          h_parm = H_PARM,
                         frailty = FRAILTY)

vacc.dist <- popavg_dist(      x = time.dist, 
                           knots = tvec.in,
                           logk0 = logk0v,
                              g0 = g0v,
                       delta_vec = delta_vec_v,
                          h_parm = H_PARM,
                         frailty = FRAILTY)
  
plot  (time.dist, plac.dist, type="l", col="black", xlab="Days", ylab="CDF")
points(time     , F0       , col="black", pch=15)
lines (time.dist, vacc.dist, col="blue")
points(time     , F1, col="blue", pch=16)
abline(v=tvec.in, lty=2, col="grey")
legend("topleft", c("Placebo", "Vaccine"),
       col=c("black","blue"),
       lty=1,
       pch=c(15,16)
)
axis(1, at=c(tvec.in[1]))
```

![](alpha_files/figure-html/alpha-02-knots-popavg-dist-1.png)

Note from Days 0 to 5 the curves are the same. Over this range the HR
will be 1. Otherwise the fit looks okay – some points are below their
line and some are above.

#### Plot population avg HR

- uses parameters from two step fitting procedure
- HR(0) = 1

We use the parameters of the previous fits to estimate the population
average hazard, which we can fit with
[`mpw::popavg_haz()`](https://swihart.github.io/mpw/reference/popavg_haz.md).
We adjust the time input so that the first element is just above 0.

``` r
time.haz <- c(1e-4, time.dist[-1])
plac.haz <- popavg_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0p,
                              g0 = g0p,
                       delta_vec = delta_vec_p,
                          h_parm = H_PARM,
                         frailty = FRAILTY)

vacc.haz <- popavg_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0v,
                              g0 = g0v,
                       delta_vec = delta_vec_v,
                          h_parm = H_PARM,
                         frailty = FRAILTY)

plot  (time.haz, vacc.haz/plac.haz, type="l", col="purple", xlab="Days",
       ylab="Hazard Ratio", lwd=2, ylim=c(-0.05,1.05))

abline(v=tvec.in, lty=2, col="grey")

axis(1, at=c(tvec.in[1]))
```

![](alpha_files/figure-html/alpha-02-knots-popavg-haz-1.png)

#### Add subject-specific HRs

- should flatten

The function
[`mpw::subjspec_haz()`](https://swihart.github.io/mpw/reference/subjspec_haz.md)
does not require the arguments `frailty` or `h_parm` because it takes
the inputs for the b,k, and delta pieces and models the subject-specific
(aka conditional) hazard. Doing this for the placebo group, then the
vaccine group and dividing them gives the subject-specific hazard ratio
(HR).

``` r
plac.haz.ss <- subjspec_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0p,
                              g0 = g0p,
                       delta_vec = delta_vec_p)

vacc.haz.ss <- subjspec_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0v,
                              g0 = g0v,
                       delta_vec = delta_vec_v)

plot  (time.haz, vacc.haz/plac.haz, type="l", col="purple", xlab="Days",
       ylab="Hazard Ratio", lwd=2, ylim=c(-0.05,1.05))

axis(1, at=c(tvec.in[1]))

lines  (time.haz, vacc.haz.ss/plac.haz.ss, type="l", col="gold", xlab="Days",
        lwd=2, lty=2)
legend("topright", c("Population HR", "Subject-specific HR"),
       col=c("purple","gold"),
       lty=c(1,2),
       lwd=2
)
```

![](alpha_files/figure-html/alpha-02-knots-subjspec-haz-1.png)

The plot above shows that for a positive-stable frailty with $\alpha$ =
0.72 the subject-specific HR is flatter and lower than the population
average HR.

### Sixteen knots

By increasing the knots, we get a “higher resolution” fit.

``` r
tvec.in <- c(5, seq(14,max(time)-14,length.out=14),max(time)-7)
tvec.in
#>  [1]   5.00000  14.00000  26.92308  39.84615  52.76923  65.69231  78.61538
#>  [8]  91.53846 104.46154 117.38462 130.30769 143.23077 156.15385 169.07692
#> [15] 182.00000 189.00000
```

#### Two step fitting procedure

##### 1. Fit for placebo group

``` r
fit_F_alpha <- function(x){mean((popavg_dist(time, 
                                            knots= tvec.in, 
                                            logk0=x[1], 
                                            g0=x[2], 
                                            delta_vec=x[-c(1,2)],
                                            h_parm = H_PARM,
                                            frailty=FRAILTY)  - F0)^2) }

init.vals <- c(log(2.4), -10, rep( 0,length(tvec.in)))
plac_fit_F_alpha <- 
  optim(init.vals, fit_F_alpha,
        method="Nelder-Mead",
        control=list(maxit=1e8))
print(plac_fit_F_alpha)
#> $par
#>  [1]  -0.19809428 -10.00290490  -0.46528979   1.50975780  -0.41049204
#>  [6]   0.13659323   0.10835991   0.26344501   0.18239172  -0.08269246
#> [11]  -0.07470972   0.19303094   0.05622029  -0.55575772  -0.99909140
#> [16]   0.53537197   0.38629201   0.88779139
#> 
#> $value
#> [1] 2.366405e-07
#> 
#> $counts
#> function gradient 
#>      909       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL

logk0p = plac_fit_F_alpha$par[1]
g0p = plac_fit_F_alpha$par[2]
delta_vec_p =plac_fit_F_alpha$par[-1*c(1,2)]
```

##### 2. Fit for vaccine group

Now take `logk0p` and `g0p`, the parameters that control the first
(leftmost) piece for the placebo group and use them for `logk0v` and
`g0v` when estimating the pieces for the vaccine group. This forces the
first pieces to be the same for both groups (i.e. *ramp-up time*).

``` r
fit_F_alpha <- function(x){mean((popavg_dist(time, 
                                            knots= tvec.in, 
                                            logk0=logk0p, 
                                            g0=g0p, 
                                            delta_vec=x,
                                            h_parm = H_PARM,
                                            frailty=FRAILTY)  - F1)^2) }

init.vals <- rep( 0,length(tvec.in))
vacc_fit_F_alpha <- 
  optim(init.vals, fit_F_alpha,
        method="Nelder-Mead",
        control=list(maxit=1e8))
print(vacc_fit_F_alpha)
#> $par
#>  [1] -0.80982557 -0.01045652  0.62457154 -0.33277056 -0.21486434  0.42390586
#>  [7]  1.31843225 -0.36908115  1.26785248 -0.74096065  0.23847180 -0.30173673
#> [13] -0.44639190 -0.40910808  1.19813727  0.38037584
#> 
#> $value
#> [1] 1.488447e-08
#> 
#> $counts
#> function gradient 
#>     2483       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL

logk0v = logk0p
g0v = g0p
delta_vec_v =vacc_fit_F_alpha$par
```

#### Plot population avg CDF

- Quick check of fit

We use the knots and the original time in `time.dist` to do a quick
check of the fit. We put vertical dashed lines in the plot to denote
knot placement.

``` r
time.dist <- sort(c(seq(min(c(time, tvec.in)), 
                        max(c(time, tvec.in)), 0.1),
                    tvec.in)
                  )

plac.dist <- popavg_dist(      x = time.dist, 
                           knots = tvec.in,
                           logk0 = logk0p,
                              g0 = g0p,
                       delta_vec = delta_vec_p,
                          h_parm = H_PARM,
                         frailty = FRAILTY)

vacc.dist <- popavg_dist(      x = time.dist, 
                           knots = tvec.in,
                           logk0 = logk0v,
                              g0 = g0v,
                       delta_vec = delta_vec_v,
                          h_parm = H_PARM,
                         frailty = FRAILTY)
  
plot  (time.dist, plac.dist, type="l", col="black", xlab="Days", ylab="CDF")
points(time     , F0       , col="black", pch=15)
lines (time.dist, vacc.dist, col="blue")
points(time     , F1, col="blue", pch=16)
abline(v=tvec.in, lty=2, col="grey")
legend("topleft", c("Placebo", "Vaccine"),
       col=c("black","blue"),
       lty=1,
       pch=c(15,16)
)
axis(1, at=c(tvec.in[1]))
```

![](alpha_files/figure-html/alpha-16-knots-popavg-dist-1.png)

Note from Days 0 to 5 the curves are the same. Comparing with the two
knot example, the fit is much improved.

#### Plot population avg HR

- uses parameters from two step fitting procedure
- HR(0) = 1

We use the parameters of the previous fits to estimate the population
average hazard, which we can fit with
[`mpw::popavg_haz()`](https://swihart.github.io/mpw/reference/popavg_haz.md).
We adjust the time input so that the first element is just above 0.

``` r
time.haz <- c(1e-4, time.dist[-1])
plac.haz <- popavg_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0p,
                              g0 = g0p,
                       delta_vec = delta_vec_p,
                          h_parm = H_PARM,
                         frailty = FRAILTY)

vacc.haz <- popavg_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0v,
                              g0 = g0v,
                       delta_vec = delta_vec_v,
                          h_parm = H_PARM,
                         frailty = FRAILTY)

plot  (time.haz, vacc.haz/plac.haz, type="l", col="purple", xlab="Days", ylab="Hazard Ratio", lwd=2, ylim=c(-0.05,1.05))

abline(v=tvec.in, lty=2, col="grey")

axis(1, at=c(tvec.in[1]))
```

![](alpha_files/figure-html/alpha-16-knots-popavg-haz-1.png)

#### Add subject-specific HRs

- should flatten

The function
[`mpw::subjspec_haz()`](https://swihart.github.io/mpw/reference/subjspec_haz.md)
does not require the arguments `frailty` or `h_parm` because it takes
the inputs for the b,k, and delta pieces and models the subject-specific
(aka conditional) hazard. Doing this for the placebo group, then the
vaccine group and dividing them gives the subject-specific hazard ratio
(HR).

``` r
plac.haz.ss <- subjspec_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0p,
                              g0 = g0p,
                       delta_vec = delta_vec_p)

vacc.haz.ss <- subjspec_haz(      x = time.haz, 
                           knots = tvec.in,
                           logk0 = logk0v,
                              g0 = g0v,
                       delta_vec = delta_vec_v)

plot  (time.haz, vacc.haz/plac.haz, type="l", col="purple", xlab="Days",
       ylab="Hazard Ratio", lwd=2, ylim=c(-0.05,1.05))

axis(1, at=c(tvec.in[1]))

lines  (time.haz, vacc.haz.ss/plac.haz.ss, type="l", col="gold", xlab="Days",
        lwd=2, lty=2)
legend("topright", c("Population HR", "Subject-specific HR"),
       col=c("purple","gold"),
       lty=c(1,2),
       lwd=2
)
```

![](alpha_files/figure-html/alpha-16-knots-subjspec-haz-1.png)

The plot above shows that for a positive-stable frailty with $\alpha$ =
0.72 the subject-specific HR is flatter and lower than the population
average HR.
