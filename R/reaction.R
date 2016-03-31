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
      as.list(new_components) %>% as_data_frame() %>%
        gather(new_component, new_comp_stoic),
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

  # figure out where to position from existing components
  new_ab <- with(rxns[1,], {
    if (comp_stoic < 0 && new_comp_stoic < 0) return(abscissa) # consumed in both reactions
    else if (comp_stoic < 0 && new_comp_stoic > 0) return(abscissa - 1) # produced in new reaction but consumed in prior
    else if (comp_stoic > 0 && new_comp_stoic < 0) return(abscissa + 1) # produced in old but consumed in new reaction
    else if (comp_stoic > 0 && new_comp_stoic > 0) return (abscissa) # produced in both reactions
    else stop ("should never happen")
  })
  return(new_ab)
}

#' Create a reaction for a pathway
#' @param eq valid equation with format a * A + b * B + ... == x * X + y * Y + ...
#' @param flux the net flux through the reaction - can be a number or expression (referencing variables and parameters in the system)
#' @param ... isotopic composition of the flux transferred in the reaction - can be a number or expression (referencing variables and parameters in the system), naming convention: flux.[<component>.]<isotope> = ... (component can be omitted if the isotopic composition of the flux is the same for each pool)
#' @param abscissa the reaction loaction (purely for sorting in a digram), a value will be inferred from the system if NULL
#' @export
add_reaction <- function(ip, name, eq, flux = NULL, ..., abscissa = NULL) {
  add_reaction_(ip, name, deparse(substitute(eq)),lazy(flux, env = parent.frame()), isotopes = lazy_dots(...), abscissa = abscissa)
}

#' add reaction with standard evaluation
#'
#' @param abscissa the reaction loaction (purely for sorting in a digram), a value will be inferred from the system if NULL
#' @param flux lazy object for later evaluation
#' @param isotopes named list with lazy objects for later evaluation
#' @export
#' @note standard evaluation
add_reaction_ <- function(ip, name, eq, flux, isotopes = list(), abscissa = NULL) {
  if (!is(ip, "isopath")) stop ("reaction can only be added to an isopath", call. = FALSE)

  components <- parse_reaction_equation(eq)
  if (is.null(abscissa)) {
    abscissa <- default_abscissa(ip, components)
  }

  ip$reactions[[name]] <-
    list(
      name = name,
      abscissa = abscissa,
      components = components,
      flux = flux,
      isotopes = list()
    )

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

  return(ip)
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
    filter(!is.na(iso_stoic), !is.na(comp_stoic))
}

#' parse reaction component
#' pulls out the stoichiometry and component
#' @param x the stoichiometric component, e.g. "A", "2 * A", etc.
#' @note standard evaluation
parse_reaction_component <- function(x) {
  m <- regexec("^\\s*(\\d*)( \\* )?(\\w+)\\s*$", x)
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
  if (length(sides) != 2) stop(err) # missing two sides
  left <- strsplit(sides[1], "+", fixed = T)[[1]]
  right <- strsplit(sides[2], "+", fixed = T)[[1]]
  if (length(left) == 0 || length(right) == 0) stop(err)
  c(
    lapply(left, parse_reaction_component) %>% unlist() * -1,
    lapply(right, parse_reaction_component) %>% unlist() * 1
  )
}
