# Population average hazard function

Population average hazard function

## Usage

``` r
popavg_haz(
  x,
  knots,
  logk0 = -0.1,
  g0 = 0.1,
  delta_vec,
  h_parm = 1,
  frailty = "PS",
  xi = 0.99,
  n = 0.99,
  s = 1.01
)
```

## Arguments

- x:

  The domain of the function; the survival times.

- knots:

  A vector of the internal knots or change-points. Must be within the
  range (min(x), max(x)).

- logk0:

  The log of the shape parameter for the first piece.

- g0:

  The shape parameter.

- delta_vec:

  the vector of delta parameters, in order.

- h_parm:

  The heterogeneity parameter or homogeneity parameter. You decide what
  h stands for.

- frailty:

  character matching "PS" for positive stable (alpha case in paper),
  "GA" for gamma (beta), "IG" for Inverse Gaussian (lambda), "TP1" (rho)
  for two point with mean of 1, and "TPU" for two point unrestricted
  (omega).

- xi:

  a scalar less than 1 for the `frailty="TP1"` case.

- n:

  a scalar less than 1 for the `frailty="TPU"` case.

- s:

  a scalar more than 1 for the `frailty="TPU"` case.

## Value

Value

## Examples

``` r
rnorm(1)
#> [1] -0.005571287
```
