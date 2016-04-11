#' Calculate forward/reverse flux
#'
#' Convenience function to calculate forward and reverse flux from net flux and
#' reversibility.
#'
#' @param net the net flux
#' @param reversibility the reversibility of the reaction
#' @param direction whether to go forward ("+") or reverse ("-")
#' @param model_offset completely reversible reactions (\code{reversibility = 1})
#' do not have any net flux so if net is != 0 this is a problem that can be solved
#' at steady-state but not during a model. The \code{model_offset} parameter offsets
#' completely reversible reactions from 1 by this small number to make the reaction
#' computationally feasible.
#' @export
#' @note This function uses standard evaluation.
dir_flux <- function(net, reversibility, direction, model_offset = 1e-9){

  # calculate flux
  if (direction == "+") {
    dir_flux <- net / (1 - reversibility)
    dir_sign <- 1
  } else if (direction == "-") {
    dir_flux <- net / (1 - reversibility) * reversibility
    dir_sign <- -1
  } else
    stop("direction not recognized: ", direction, call. = FALSE)

  # check for problems (this takes some processing - any way to get around even more?)
  if (any(is.infinite(dir_flux))) {
    # this is from reversibility == 1 cases, add the model_offset
    reversibility[reversibility == 1 & net > 0] <- reversibility[reversibility == 1 & net > 0] - model_offset
    reversibility[reversibility == 1 & net < 0] <- reversibility[reversibility == 1 & net < 0] + model_offset
    if (direction == "+") dir_flux <- net / (1 - reversibility)
    else if (direction == "-") dir_flux <- net / (1 - reversibility) * reversibility
  }

  if (any(dir_flux < 0, na.rm = TRUE)) {
    stop("negative directional flux does not make sense")
  }

  return(dir_sign * dir_flux)
}

#' Fractionate a delta value
#'
#' Convenience function to fractionate a delta value using a fractionation factor
#' in alpha or epsilon notation. This function does not approximate but calculates accurate
#' fractionation.
#'
#' @param delta The delta value (has to be in permil if \code{permil = TRUE})
#' @param alpha The fractionation factor in alpha notation.
#' @param epsilon The fractionation factor in epsilon notation (alpha - 1)
#' (has to be in permil if \code{permil = TRUE}).
#' @param permil Whether the delta and epsilon values are in permil notation
#' (i.e. all multiplied by 1000) or not.
#' @param multiply Whether to multiply the delta-derived ratio by the fractionation factor or divide it.
#' Divide is \code{multiply=FALSE} (the default), which is the correct behavior for kinetic fractionation
#' factors defined as k_light/k_heavy and equilibrium fractionation factors defined as substrate/product.
#' Multiply is \code{multiply=TRUE}, which is the correct behavior for kinetic fractionation
#' factors defined as k_heavy/k_light and equilibrium fractionation factors defined as product/substrate.
#' @note This function uses standard evaluation.
#' @export
#' @return The resulting delta value (in permil notation if \code{permil = TRUE})
fractionate <- function(delta, alpha = NULL, epsilon = NULL, permil = TRUE, multiply = FALSE) {

  multiplier <- if (permil) 1000 else 1

  if (is.null(alpha) && is.null(epsilon))
    stop("cannot calculate fractionation, no fractionation factor provided",
         call. = FALSE)
  else if (!is.null(alpha) && !is.null(epsilon))
    stop("cannot calculate fractionation, ambiguous fractionation factor provided as both alpha and epsilon value",
         .call = FALSE)
  else if (!is.null(epsilon))
    alpha <- epsilon/multiplier + 1

  if (any(alpha <= 0) || any(alpha > 2))
    stop("unrealistically high (>2) or impossibly low (<0) fractionation factors. ",
         "Make sure to indicate properly whether parameters to fractionate are alpha or eps: alpha = ",
         alpha %>% signif(3) %>% paste(collapse = ", "), call. = FALSE)

  # multiplication is computationally faster so switching to this
  alpha <- if (!multiply) 1/alpha else alpha

  # calculate delta value
  # (note: rather than multiplier ((d/multipler + 1) * alpha - 1), this is multiplied
  # through by multiplier because it is computationally ~30% faster!
  return((delta + multiplier) * alpha - multiplier)
}
