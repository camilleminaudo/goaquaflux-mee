###############################################################################
## synthetic_incubations.R
##
## Generator of synthetic floating-chamber CH4 incubations with a known
## ground truth (diffusive rate, bubble times and bubble sizes). Sourced by
##   * compare_goAquaFlux_FluxSeparator_MSP.R (benchmark, Appendix S1),
##   * make_FigS1_synthetic_incubations.R     (Figure S1).
##
## Model (Appendix S1, Eqs. S3-S5), at a regular time step dt:
##
##   C(t) = C0 + s * t + sum_i b(t - t_i; M_i, A_i) + e(t)
##
##   b(u) = 0                                  for u < 0
##        = (M + A) * u / r                    for 0 <= u < r   (ramp)
##        = M + A * exp(-(u - r) / tau)        for u >= r       (decay)
##
## with M_i the settled step (gas added to the headspace once mixed),
## A_i = overshoot * M_i the transient excess seen by the analyser before the
## bubble gas is mixed, r the ramp time and tau the decay time constant.
## Bubble sizes are log-normal with arithmetic mean `bubble_mean`; onsets are
## uniform in [edge, duration - edge] with a minimum spacing; the noise e(t)
## is white Gaussian (or AR(1) if s$ar1 > 0).
###############################################################################


#' Concentration response to one bubble
#'
#' @param u Numeric vector; time since the bubble onset (s).
#' @param M Numeric; settled step (gas units).
#' @param A Numeric; transient overshoot at the peak (gas units).
#' @param tau Numeric; decay time constant of the overshoot (s).
#' @param ramp Numeric; time from onset to the transient peak (s).
#' @return Numeric vector of concentration increases, same length as `u`.
bubble_shape <- function(u, M, A, tau, ramp) {
  out <- numeric(length(u))
  if (ramp > 0) { r <- u >= 0 & u < ramp; out[r] <- (M + A) * u[r] / ramp }
  p <- u >= ramp
  out[p] <- M + A * exp(-(u[p] - ramp) / tau)
  out
}


#' Random bubble onsets with a minimum spacing (rejection sampling)
#'
#' @param n Integer; number of bubbles.
#' @param lower,upper Numeric; interval of possible onsets (s).
#' @param min_spacing Numeric; minimum time between two onsets (s).
#' @param max_try Integer; maximum number of draws.
#' @return Sorted numeric vector of onsets, rounded to the second.
place_bubbles <- function(n, lower, upper, min_spacing, max_try = 5000) {
  if (n == 0) return(numeric(0))
  for (k in seq_len(max_try)) {
    tb <- sort(round(runif(n, lower, upper)))
    if (n == 1 || all(diff(tb) >= min_spacing)) return(tb)
  }
  stop("Could not place ", n, " bubbles")
}


#' Simulate one synthetic incubation
#'
#' The random draws are made in a fixed order (onsets, sizes, noise), so that
#' a given seed always yields the same incubation.
#'
#' @param sc One-row data.frame (or list) with the scenario: noise_sd (ppb),
#'   slope (ppb s-1), bubble_mean (ppb), n_bubbles, overshoot (A / M).
#' @param s List of constants: duration, dt, C0, ramp, tau, edge,
#'   min_spacing, size_sdlog, ar1 (see CFG$syn in
#'   compare_goAquaFlux_FluxSeparator_MSP.R).
#' @return List with t (s), y (ppb, with noise), tb (true onsets, s) and
#'   M (true settled steps, ppb).
simulate_series <- function(sc, s) {
  t  <- seq(0, s$duration - s$dt, by = s$dt)
  tb <- place_bubbles(sc$n_bubbles, s$edge, s$duration - s$edge, s$min_spacing)
  M  <- if (sc$n_bubbles > 0) sc$bubble_mean *
    exp(rnorm(sc$n_bubbles, -s$size_sdlog^2 / 2, s$size_sdlog)) else numeric(0)
  y <- s$C0 + sc$slope * t
  for (i in seq_along(tb)) {
    y <- y + bubble_shape(t - tb[i], M[i], sc$overshoot * M[i], s$tau, s$ramp)
  }
  eps <- if (s$ar1 > 0) {
    as.numeric(arima.sim(list(ar = s$ar1), n = length(t),
                         sd = sc$noise_sd * sqrt(1 - s$ar1^2)))
  } else {
    rnorm(length(t), 0, sc$noise_sd)
  }
  list(t = t, y = y + eps, tb = tb, M = M)
}
