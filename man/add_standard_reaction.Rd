% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reaction.R
\name{add_standard_reaction}
\alias{add_standard_reaction}
\title{Add a standard single fractionation reaction}
\usage{
add_standard_reaction(
  ip,
  equation,
  name = default_rxn_name(ip),
  flux = NULL,
  reversibility,
  ...,
  permil = TRUE,
  eq_ratio = c("P/S", "S/P"),
  abscissa = NULL
)
}
\arguments{
\item{ip}{An \link{isopath} object to which this reaction should be added.}

\item{name}{Name of the reaction. Must be unique (i.e. will overwrite an existing reaction of the same name). If omitted (the default), will name the reaction based on the existing number of reactions in the \link{isopath}, e.g. "rxn4".}

\item{flux}{The net flux through the reaction - can be a number or expression (referencing variables and parameters in the system)}

\item{reversibility}{The reversibility of the reaction. If this is omitted, the reaction is treated as being completely irreversible if a kinetic fractionation factor (\code{alpha.<isotope>} or \code{eps.<isotope>}) is provided or as in full equilibrium if an equilibrium fractionation factor (\code{alpha.<isotope>.eq} or \code{eps.<isotope>.eq}) is provided. If anything else is provided (can be a number or expression; evaluation to 0 implies no reversibility, evaluation to 1 implies full reversibility), requires either a fractionation factor for the reverse reaction or an equilibrium fractionation factor to be defined in the dots (\code{...}) in addition to the forward kinetic effect for each isotope system involved in the reaction.}

\item{...}{Fractionation factors for the reaction - can be numbers or expressions (referencing other variables and parameters in the system), naming convention: \code{alpha.<isotope> = ...}, or \code{eps.<isotope> = ...} for the kinetic isotope effect of the reaction in alpha or epsilon (permil) notation (with definition alpha = k_light / k_heavy and eps = (alpha - 1) * 1000, i.e. normal isotope effects are alpha > 1 and eps > 0, inverse isotope effects are alpha < 1 and eps < 0), and \code{alpha.<isotope>eq ...} or \code{eps.<isotope> = ...} for equilibrium isotope effects. If only kinetic fractionation factors are provided the reaction is created as an irreversible reaction, if only equilibrium fractionation factors are provided, it's created as an equilibrium reaction. For anything in between, use the \code{reversibility} parameter, If reversibility is set, a fractionation factor for the reverse reaction must be provided either directly (\code{alpha.<isotope>.rev = } or \code{eps.<isotope>.rev = ...}) or in the form of the equilibrium fractionation factor (\code{alpha.<isotope>.eq = ...} or \code{eps.<isotope>.eq = ...}). All equilibrium fractionatino factors are interpreted as product/substrate or substrate/product depending on the \code{eq_ratio} parameter..}

\item{permil}{Whether the delta and epsilon values are in permil notation (i.e. all multiplied by 1000) or not. All must be in the same notation.}

\item{eq_ratio}{eq_ratio Whether equilibrium fractionation factors are provided as product / substrate ratios (the default) or substrate / product ratios.}

\item{abscissa}{The reaction loaction. This is used purely for graphical representation (see \link{generate_reaction_diagram}). By default, a decent value will be guessed from the existing reactions in the \link{isopath}. This works fairly well as long as all reactions after the first one are added to the isopath in an order such that at least one reactant or product is already part of a different reaction in the pathway.}
}
\description{
This adds a standard 1 to 1 (one reactant, one product) reaction with standard
fractionation factors to the isopath (for all isotopes that are part of both
reactant AND product). Supports reversibility for easy implementation of bi-
directional fluxes. This is really intended as a convenience function to make
it easy to set up Hayes-type reaction systems. Reaction systems that only have
this kind of reaction also have the added benefit of having an analytical
solution for steady state. See \link{calculate_steady_state()} to make use of
this functionalty.
}
\seealso{
Other reaction functions: 
\code{\link{add_custom_reaction}()}
}
\concept{reaction functions}
