# Subject specific hazard function

Take the parameters estimated from `popavg_dist` and calculate subject
specific hazards of the Weibull form `b*k*x^(k-1)`.

## Usage

``` r
subjspec_haz(x, knots, logk0 = -0.1, g0 = 0.1, delta_vec)
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

## Value

Value

## Examples

``` r
rnorm(1)
#> [1] 1.148412
```
