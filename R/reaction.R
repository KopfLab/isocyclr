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
