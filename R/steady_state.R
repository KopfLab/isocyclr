#' to evaluate steady state for a chain of reactions, need to
#'  - dissect all the flux isotopic composition calls $reactions$isotopes$x to see that it they all make function calls to "reversible_reaction" or "irreversible_reaction" and dissect out their parameters to get at the reversibility and alphas
#'  - make sure that only the starting component is fixed (no other components may be!)
#'  - that all reactions are only A to B! no multi component shenanigans
#'  - that the fluxes at branching points can be interpreted as branching ratios
#'  - the dtemplate data frame is simpler?
#'  - continue here...

#' actually just keep track of what type of reactions are what and do the steady-state for the hayes type models
#' --> only allow steady-state if it's all hayes type reactions
#'


#' for hays type things
#' pass in the isotopic composition of the flux INTO the system, make it clear that there
#' is flux OUT of the system into an external reservoir that is not part of the pathway
#'
#' has to double check that all reactions are of type "simple" and that there are no converging
#' reactant paths in any of the isotopes (diverging is okay)
#'
#' has to check that all diverging pathways have defined fluxes to be able to take the ratio!!
#'
#' has to check that there is an obvious starting point in the reactio network (otherwise error)
#'
#' when putting together the expressions, do it recusrively from start point to walk out to all the
#' branch tipping points (randomly select one branch over the other as the final points in the
#' network)
#'
#' once the expressions are generated, should be able to just evaluate in the context of the parameter
#' data frame, i.e for all scenarios at once!!
#'

