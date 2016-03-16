#' Create a reaction for a pathway
#' @param eq valid equation with format a * A + b * B + ... == x * X + y * Y + ...
#' @param flux the net flux through the reaction - can be a number or expression (referencing variables and parameters in the system)
#' @param ... isotopic composition of the flux transferred in the reaction - can be a number or expression (referencing variables and parameters in the system), naming convention: flux.[<component>.]<isotope> = ... (component can be omitted if the isotopic composition of the flux is the same for each pool)
#' @export
add_reaction <- function(ip, name, eq, flux, ..., nr = new_nr()) {
  new_nr <- function() length(ip$reactions) + 1
  add_reaction_(ip, name, deparse(substitute(eq)), nr, lazy(flux, env = parent.frame()), isotopes = lazy_dots(...))
}

#' add reaction with standard evaluation
#'
#' @param nr the reaction number (for sorting it in a diagram)
#' @param flux lazy object for later evaluation
#' @param isotopes named list with lazy objects for later evaluation
#' @export
#' @note standard evaluation
add_reaction_ <- function(ip, name, eq, nr, flux, isotopes = list()) {
  if (!is(ip, "isopath")) stop ("reaction can only be added to an isopath", call. = FALSE)

  new_nr <- function() {
    length(ip$reactions) + 1
  }

  ip$reactions[[name]] <-
    list(
      name = name,
      nr = nr,
      components = parse_reaction_equation(eq),
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

  return(ip)
}

#' get the reaction matrix for an isopath
#' @export
get_reaction_matrix <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get reaction matrix from an isopath")
  lapply(ip$reactions, function(i) {
    c(list(reaction = i$name, rxn_nr = i$nr), as.list(i$components)) %>% as_data_frame()
  }) %>%
    bind_rows()
}

#' get the combined reaction_component_matrix
#' @export
get_reaction_component_matrix <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get reaction component matrix from an isopath")
  ip %>%
    get_reaction_matrix() %>%
    gather(component, comp_stoic, -reaction, -rxn_nr) %>%
    left_join(
      ip %>% get_component_matrix() %>%
        gather(isotope, iso_stoic, -component),
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
    lapply(left, parse_reaction_component) %>% unlist() * 1,
    lapply(right, parse_reaction_component) %>% unlist() * -1
  )
}
