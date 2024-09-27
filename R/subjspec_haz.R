#' Subject specific hazard function
#'
#' Take the parameters estimated from \code{popavg_dist} and calculate
#' subject specific hazards of the Weibull form `b*k*x^(k-1)`.
#'
#' @param x The domain of the function; the survival times.
#' @param knots A vector of the internal knots or change-points. Must be within
#' the range (min(x), max(x)).
#' @param logk0 The log of the shape parameter for the first piece.
#' @param g0 The shape parameter.
#' @param delta_vec the vector of delta parameters, in order.
#' @return Value
#' @export
#'
#' @examples
#' rnorm(1)
subjspec_haz <- function(x,
                        knots,
                        logk0=-.1,
                        g0=.1,
                        delta_vec){

  if( length(knots) != length(delta_vec)){
    warning("Double check the lengths of `knots` and `delta_vec` -- they should be equal.")
  }

  (result <- sapply(x,
                    function(w){suppressWarnings(
                      LI <- max((1:length(knots))[w >  knots]));
                      lower.index <- ifelse(is.infinite(LI),0,LI)

                      k0 = exp(logk0);

                      k =  k0+sum(delta_vec[0: lower.index])

                      b =  exp(g0 -
                                 sum(delta_vec[0: lower.index]*
                                       log(knots[0: lower.index])))

                      b*k*x^(k-1)

                    }
  )
  )

  ifelse(length(result)<=1, m<-result, m<-diag(result))
  m
}
