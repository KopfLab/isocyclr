#' Fractionate a delta value
#'
#' Convenience function to fractionate a delta value using a fractionation factor
#' in alpha or epsilon notation. This function does not approximate calculate accurate
#' fractionation.
#'
#' @param delta the delta value (has to be in permil if \code{permil = TRUE})
#' @param alpha the raw fractionation factor
#' @param epsilon the fractionation factor in epsilon notation (alpha - 1)
#' (has to be in permil if \code{permil = TRUE})
#' @param permil whether the delta and epsilon values are in permil notation
#' (i.e. all multiplied by 1000) or not (default is TRUE)
#' @param multiply whether to multiply the resulting ratio by the fractionation factor or divide it.
#' Divide is \code{multiply=FALSE} (the default), which is the correct behavior for kinetic fractionation
#' factors defined as k_light/k_heavy and equilibrium fractionation factors defined as reactant/product.
#' Multiply is \code{multiply=TRUE}, which is the correct behavior for kinetic fractionatoin
#' factors defined as k_heave/k_light and equilibrium fractionation factors defined as product/reactant.
#' @note This function uses standard evaluation.
#' @return the resulting delta value (in permil notation if \code{permil = TRUE})
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

  # multiplication is computationally faster so switching to this
  alpha <- if (!multiply) 1/alpha else alpha

  # calculate delta value
  # (note: rather than multiplier ((d/multipler + 1) * alpha - 1), this is multiplied
  # through by multiplier because it is computationally ~30% faster!
  return((delta + multiplier) * alpha - multiplier)
}
