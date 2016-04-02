#' Add a custom reaction
#'
#' Anything's allowed for these, just the flux and the isotopic composition of the flux are specified.
#'
#' @param ip An \link{isopath} object to which this reaction should be added.
#' @param name Name of the reaction. Must be unique (i.e. will overwrite an existing reaction of the same name). If omitted (the default), will name the reaction based on the existing number of reactions in the \link{isopath}, e.g. "rxn4".
#' @param eq Valid equation in chemical format a * A + b * B + ... == x * X + y * Y + ... Note that all components must already be part of the \link{isopath} (add components using \link{add_component}).
#' @param flux The net flux through the reaction - can be a number or expression (referencing variables and parameters in the system)
#' @param ... isotopic composition of the flux transferred in the reaction - can be a number or expression (referencing variables and parameters in the system), naming convention: flux.[<component>.]<isotope> = ... (component can be omitted if the isotopic composition of the flux is the same for each pool)
#' @param abscissa The reaction loaction. This is used purely for graphical representation (see \link{generate_reaction_diagram}). By default, a decent value will be guessed from the existing reactions in the \link{isopath}. This works fairly well as long as all reactions after the first one are added to the isopath in an order such that at least one reactant or product is already part of a different reaction in the pathway.
#' @return The \link{isopath} object with the new reaction added.
#' @family reaction functions
#' @examples
#' isopath() %>% add_component(c("A", "B")) %>% add_custom_reaction(A == 3 * B)
#' @note Under the hood, this uses the standard evaluation function \link{add_reaction_}, which can be used directly if needed.
#' @export
add_custom_reaction <- function(ip, eq, name = default_rxn_name(ip), flux = NULL, ..., abscissa = NULL) {
  add_reaction_(ip, deparse(substitute(eq)), name = name, flux = lazy(flux),
                isotopes = lazy_dots(...), class = "custom", abscissa = abscissa)
}

#' Add a simple single fractionation reaction
#'
#' This adds a reaction
#'
#' - check later that fluxes all add up correctly
#' - check later that the starting point is non-variable - is that actually required? prob not!
#'
#' @param reversibility The reversibility of the reaction. If this is omitted, the reaction is treated as being irreversible. If anything else is provided (can be a number or expression evaluating to 0 to 1), requires either a fractionation factor for the reverse reaction or an equilibrium fractionation factor to be part of the dots (\code{...}) for each isotope system involved in the reaction.
#' @param ... Fractionation factors for the reaction - can be numbers or expressions (referencing other variables and parameters in the system), naming convention: \code{alpha.<isotope> = ...}, or \code{eps.<isotope> = ...} for the kinetic isotope effect of the reaction in alpha or epsilon (permil) notation (with definition alpha = k_light / k_heavy and eps = (alpha - 1) * 1000, i.e. normal isotope effects are alpha > 1 and eps > 0, inverse isotope effects are alpha < 1 and eps < 0). If reaction is reversible (\code{reversible = TRUE}), a fractionation factor for the reverse reaction must be provided either directly (\code{alpha.<isotope>.rev = } or \code{eps.<isotope>.rev = ...}) or in the form of an equilibrium fractionation factor (\code{alpha.<isotope>.eq = ...} or \code{eps.<isotope>.eq = ...}) with definition alpha.eq = R_substrate / R_product = alpha / alpha.rev).
#' @export
#' @inheritParams add_custom_reaction
#' @family reaction functions
add_simple_reaction <- function(ip, eq, name = default_rxn_name(ip), flux = NULL, reversibility, ..., abscissa = NULL) {

  # name / equation flexibility
  if (missing(eq)) {
    eq <- deparse(substitute(name))
    name <- NULL
  } else {
    eq <- deparse(substitute(eq))
  }

  # reaction components
  components <- parse_reaction_equation(eq) %>% as.list() %>% as_data_frame()
  if (length(missing_comp <- setdiff(names(components), names(ip$components))) > 0)
    stop("missing component definition(s), make sure to add this with add_component() first: ",
         missing_comp %>% paste(collapse = ", "), call. = FALSE)

  # reaction isotopes summary
  rxn_isotopes <-
    left_join(
      components %>% gather(component, comp_stoic),
      ip %>% get_component_matrix() %>% gather(isotope, iso_stoic, -component, -variable),
      by = "component") %>%
    filter(!is.na(iso_stoic)) %>%
    group_by(isotope) %>%
    mutate(n_reactant = sum(comp_stoic < 0),
           n_product = sum(comp_stoic > 0))

  # check that for each isotope it's a simple 1 to 1 transformation
  problem_isotopes <- rxn_isotopes  %>%
    filter(n_reactant > 1 | n_product > 1)
  if (nrow(problem_isotopes) > 0) {
    stop("simple reactions can only be 1 to 1 transformations for each isotope system, the following isotope system(s) have multiple reactants or products: ",
         with(problem_isotopes, paste0(isotope, " (reactants: ", n_reactant, ", products: ", n_product, ")")) %>%
           paste(collapse = ", "), call. = F)
  }

  # process only those isotopes that have exactly one reactant and one product
  rxn_isotopes <- rxn_isotopes %>%
    filter(n_reactant == 1 & n_product == 1) %>%
    select(isotope, component, comp_stoic) %>%
    spread(comp_stoic, component)

  # fractionation factors from dots
  ff_dots <- lazy_dots(...)

  # construct allowed pattern for the dots for the relevat isotopes
  ff_dots_pattern <- paste0("(alpha|eps).", rxn_isotopes$isotope) %>% setNames(rxn_isotopes$isotope)
  if (!missing(reversibility)) {
    ff_dots_pattern <- c(
      ff_dots_pattern,
      paste0("(alpha|eps).", isotopes, ".(rev|eq)") %>% setNames(paste0(rxn_isotopes$isotope, ".rev")))
  }

  # find matching dot names for each isotope
  ff_dots_names <- sapply(ff_dots_pattern, function(i){
    ff_name <- gsub(".", "\\.", paste0("^", i, "$"), fixed = T) %>% grep(names(ff_dots), value = T)
    if (length(ff_name) == 0) return (NA)
    else if (length(ff_name) > 1) stop("shouldn't be possible to happen")
    return (ff_name)
  })

  # stop if any of the required ff dots are missing
  ff_dots_missing <- ff_dots_pattern[is.na(ff_dots_names)]
  if (length(ff_dots_missing) > 0) {
    stop("not all required fractionation factors available, missing: ",
         ff_dots_missing %>% paste(collapse = ", "), call. = F)
  }

  # loop through each isotope to add reaction
  for (i in 1:nrow(rxn_isotopes)) {
    with(rxn_isotopes[i,],{
      message("isotope: ", isotope, " reactant: ", `-1`, " product: ", `1`)
    })
  }

  return(invisible(ip))
}

#' get default reaction name
default_rxn_name <- function(ip) {
  paste0("rxn", length(ip$reactions) + 1)
}

#' add reaction with standard evaluation
#'
#' @param flux lazy object for later evaluation
#' @param isotopes named list with lazy objects for later evaluation to calculate the isotope composition of the flux through this reaction
#' @param abscissa the reaction loaction (purely for sorting in a digram), a value will be inferred from the system if NULL
#' @param class the class of the reaction for internal processing purposes (e.g. steady state calculations)
#' @param args additional arguments (usually expressions)
#' @export
#' @note standard evaluation
add_reaction_ <- function(ip, eq, name, flux, isotopes = list(),
                          abscissa = NULL, class = "generic", args = list()) {

  if (!is(ip, "isopath")) stop ("reaction can only be added to an isopath", call. = FALSE)

  # default abscissa
  components <- parse_reaction_equation(eq)
  if (is.null(abscissa))
    abscissa <- default_abscissa(ip, components)

  ip$reactions[[name]] <-
    structure(
      list(
        name = name,
        abscissa = abscissa,
        components = components,
        flux = flux,
        isotopes = list(),
        args = args
      ), class = class)

  # check isotope names for the flux. prefix
  missing_prefix <- names(isotopes)[!grepl("flux\\.", names(isotopes))]
  if (length(missing_prefix) > 0)
    stop("missing prefix for isotope flux, please prefix each flux isotopic composition with 'flux.[<component>.].<isotope> = ...': ", missing_prefix %>% paste(collapse = ", "), call. = FALSE)
  names(isotopes) <- gsub("flux\\.", "", names(isotopes))

  # parse isotopes list for compound names (isotope.component)
  isotope.components <- c()
  for (iso in names(isotopes)) {
    parts <- strsplit(iso, ".", fixed = T)[[1]]
    if ( length(parts) == 1) {
      ip$reactions[[name]]$isotopes[[parts[1]]] <- isotopes[[iso]]
    } else if (length(parts) == 2) {
      isotope.components <- c(isotope.components, parts[1])
      ip$reactions[[name]]$isotopes[[parts[2]]][[parts[1]]] <- isotopes[[iso]]
    } else stop("cannot process compound isotope.component name: ", iso, call. = FALSE)
  }

  # check for missing components
  missing_comp <- c(
    setdiff(ip$reactions[[name]]$components %>% names(), names(ip$components)),
    setdiff(isotope.components, names(ip$components))
  )
  if (length(missing_comp) > 0)
    stop("missing component definition(s), make sure to add this with add_component() first: ",
         missing_comp %>% paste(collapse = ", "), call. = FALSE)

  # check for missing isotopes
  missing_isos <- setdiff(ip$reactions[[name]]$isotopes %>% names(), names(ip$isotopes))
  if (length(missing_isos) > 0)
    stop("missing isotope definition(s), make sure to add this with add_isotope() first: ",
         missing_isos %>% paste(collapse = ", "), call. = FALSE)

  # store info
  ip <- ip %>% store_info()

  return(invisible(ip))
}

#' calculate default abscissa
#' depending on exising reactions, calculate
#' where the new reaction should go
default_abscissa <- function(ip, new_components) {

  # if no other prior reactions yet
  if (length(ip$reactions) == 0)
    return(1)

  # combine new components and existing reactions
  rxns <-
    left_join(
      ip$info$reaction_component_matrix,
      new_components %>% as.list() %>% as_data_frame() %>% gather(new_component, new_comp_stoic),
      by = c("component" = "new_component")
    ) %>%
    select(abscissa, comp_stoic, new_comp_stoic) %>%
    distinct() %>%
    filter(!is.na(new_comp_stoic)) %>%
    filter(abscissa == min(abscissa)) %>%
    filter(comp_stoic == max(comp_stoic))

  # components not found
  if (nrow(rxns) == 0)
    return( max(ip$info$reaction_component_matrix$abscissa) + 1 )

  # figure out where to position reaction from existing component abscissa
  new_ab <- with(rxns[1,], {
    if (new_comp_stoic < 0) return (abscissa + 1)
    else if (new_comp_stoic > 0) return (abscissa)
    else stop ("should never happen")
  })
  return(new_ab)
}

#' get the reaction matrix for an isopath
#' @export
get_reaction_matrix <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get reaction matrix from an isopath")

  lapply(ip$reactions, function(i) {
    c(list(reaction = i$name, abscissa = i$abscissa), as.list(i$components)) %>% as_data_frame()
  }) %>%
    bind_rows()
}

#' get the combined reaction_component_matrix
#' @export
get_reaction_component_matrix <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only calculate reaction component matrix from an isopath")

  cps <- ip %>% get_component_matrix()
  rxn <- ip %>% get_reaction_matrix()

  if (nrow(cps) == 0 || nrow (rxn) == 0) return(data_frame())

  left_join(
    rxn %>% gather(component, comp_stoic, -reaction, -abscissa),
    cps %>% gather(isotope, iso_stoic, -component, -variable),
    by = "component"
  ) %>%
    # make sure to remove empty entires
    filter(!is.na(iso_stoic), !is.na(comp_stoic)) %>%
    # calculate abscissa for each component
    mutate(abscissa = abscissa - ifelse(comp_stoic < 0, 1, 0 )) %>%
    # arrange by isotope, then component abscissa order
    arrange(isotope, reaction, comp_stoic > 0, component) %>%
    # order columns
    select(isotope, reaction, component, comp_stoic, variable, abscissa)
}

#' parse reaction component
#' pulls out the stoichiometry and component
#' @param x the stoichiometric component, e.g. "A", "2 * A", etc.
#' @note standard evaluation
parse_reaction_component <- function(x) {
  m <- regexec("^\\s*(\\d+\\.?\\d*)?( \\* )?(\\w+)\\s*$", x)
  parts <- regmatches(x, m)[[1]]
  if ( length(parts) == 0)
    stop("cannot parse component: ", x, call. = F)
  value <- if (parts[2] == "") 1 else as.numeric(parts[2])
  return(value %>% setNames(parts[4]))
}

#' parse reaction equation into its components
#' @note standard evaluation
parse_reaction_equation <- function(eq) {
  sides <- strsplit(eq, "==", fixed = T)[[1]]
  err <- paste("please write equation in format 'a * A + b * B + ... == x * X + y * Y + ...':", eq)
  if (length(sides) != 2) stop(err, call. = FALSE) # missing two sides
  left <- strsplit(sides[1], "+", fixed = T)[[1]]
  right <- strsplit(sides[2], "+", fixed = T)[[1]]
  if (length(left) == 0 || length(right) == 0) stop(err, call. = FALSE)
  c(
    lapply(left, parse_reaction_component) %>% unlist() * -1,
    lapply(right, parse_reaction_component) %>% unlist() * 1
  )
}
