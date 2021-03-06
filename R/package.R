#' @name isocyclr
#' @title Model isotope fractionation in biogeochemical cycles and pathways
#' @description This package facilitates the modeling of isotopic effects through biogeochemical pathways and cycles including non-steady state and steady-state solutions.
#' @author Sebastian Kopf
#' @docType package
#' @import dplyr lazyeval deSolve rootSolve ggplot2
#' @importFrom tidyr gather gather_ extract extract_ spread spread_ pivot_longer
#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom utils modifyList
#' @importFrom rlang !! !!!
NULL

# todo list:
# - implement component straight from formula (inferring isotopes), e.g. NO3
# - implement check on mass balance consistency in equations (so that the same number of N and O get's added/removed on both sides)
# - finds ways to improve performance / make everything run faster if possible
