#' Population average hazard function
#'
#' @param x The domain of the function; the survival times.
#' @param knots A vector of the internal knots or change-points. Must be within
#' the range (min(x), max(x)).
#' @param logk0 The log of the shape parameter for the first piece.
#' @param g0 The shape parameter.
#' @param delta_vec the vector of delta parameters, in order.
#' @param h_parm The heterogeneity parameter or homogeneity parameter.  You
#' decide what h stands for.
#' @param frailty character matching "PS" for positive stable (alpha case in
#' paper), "GA" for gamma (beta),
#' "IG" for Inverse Gaussian (lambda), "TP1" (rho) for two point with mean of 1,
#' and "TPU" for two point unrestricted (omega).
#' @param xi a scalar less than 1 for the \code{frailty="TP1"} case.
#' @param n  a scalar less than 1 for the \code{frailty="TPU"} case.
#' @param s  a scalar more than 1 for the \code{frailty="TPU"} case.
#' @return Value
#' @export
#'
#' @examples
#' rnorm(1)
popavg_haz <- function(x,
                        knots,
                        logk0=-.1,
                        g0=.1,
                        delta_vec,
                        h_parm=1,
                        frailty="PS",
                        xi=0.99,
                        n= 0.99,
                        s= 1.01){

  if( length(knots) != length(delta_vec)){
    warning("Double check the lengths of `knots` and `delta_vec` -- they should be equal.")
  }
  xi_prime = xi + (1-xi)/h_parm;

  (result <- sapply(x,
                    function(w){suppressWarnings(
                      LI <- max((1:length(knots))[w >  knots]));
                      lower.index <- ifelse(is.infinite(LI),0,LI)

                      k0 = exp(logk0);

                      k =  k0+sum(delta_vec[0: lower.index])

                      b =  exp(g0 -
                                 sum(delta_vec[0: lower.index]*
                                       log(knots[0: lower.index])))
                      switch(frailty,
                             PS = b^h_parm*(h_parm*k)*x^(h_parm*k-1),
                             GA = b*k*x^(k-1) / (1 + b*x^(k)/h_parm),
                             IG = 1 -
                               exp(-h_parm * (sqrt(1 + 2/h_parm * b*x^k)-1)) +
                               any(c(k0, k0+cumsum(delta_vec)) < 0)*1e6,
                             TP1 = 1 -
                               ( (1-h_parm)*exp(-xi      * b*x^k) +
                                   (h_parm)*exp(-xi_prime* b*x^k)
                               ) + any(c(k0, k0+cumsum(delta_vec)) < 0)*1e6,
                             TPU = 1 -
                               ( (1-h_parm)*exp(-n* b*x^k) +
                                   (h_parm)*exp(-s* b*x^k)
                               ) + any(c(k0, k0+cumsum(delta_vec)) < 0)*1e6


                      )


                    }
  )
  )

  ifelse(length(result)<=1, m<-result, m<-diag(result))
  m
}
