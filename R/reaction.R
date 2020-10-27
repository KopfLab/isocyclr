
custom_rxn <- function(equation, flux = NULL, ...) {
  rxn_(deparse(substitute(equation)),
                flux = lazy(flux), isotopes = lazy_dots(...),
                class = "custom")
}

rxn_ <- function(equation, flux = NULL, isotopes = list(), class = "generic", args = list()) {

  # default abscissa
  components <- parse_reaction_equation(equation)

  rxn <-
    structure(
      list(
        equation = equation,
        components = parse_reaction_equation(equation),
        flux = flux,
        isotopes = list(),
        args = args
      ), class = c(class, "isorxn"))

  # check isotope names for the flux. prefix
  missing_prefix <- names(isotopes)[!grepl("flux\\.", names(isotopes))]
  if (length(missing_prefix) > 0)
    stop("missing prefix for isotope flux, please prefix each flux isotopic composition with 'flux.[<component>.].<isotope> = ...': ", missing_prefix %>% paste(collapse = ", "), call. = FALSE)
  names(isotopes) <- gsub("flux\\.", "", names(isotopes))

  # parse isotopes list for compound names (isotope.component)
  isotope.components <- c()
  for (iso in names(isotopes)) {
    parts <- strsplit(iso, ".", fixed = TRUE)[[1]]
    if ( length(parts) == 1) {
      rxn$isotopes[[parts[1]]] <- isotopes[[iso]]
    } else if (length(parts) == 2) {
      isotope.components <- c(isotope.components, parts[1])
      rxn$isotopes[[parts[2]]][[parts[1]]] <- isotopes[[iso]]
    } else stop("cannot process compound isotope.component name: ", iso, call. = FALSE)
  }

  return(rxn)
}


#' @export
print.isorxn <- function(x, ...) {
  sprintf(
    paste0(
      "Reaction (type '%s'):\n",
      " - Equation: %s\n",
      " - Flux: %s\n",
      " - Isotopic composition of flux:"
    ), class(x)[1], x$equation, deparse(x$flux$expr)) %>% cat()
  for (iso in names(x$isotopes)) {
    if (class(x$isotopes[[iso]]) == "lazy") {
      sprintf("\n   - %s (for all components): %s", iso,
              deparse(x$isotopes[[iso]]$expr)) %>% cat()
    } else {
      for (comp in names(x$isotopes[[iso]])) {
        sprintf("\n   - %s (for component %s): %s", iso, comp,
                deparse(x$isotopes[[iso]][[comp]]$expr)) %>% cat()
      }
    }
  }
}


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
#' @note Under the hood, this uses the standard evaluation function \link{add_reaction_}, which can be used directly if needed.
#' @export
add_custom_reaction <- function(ip, equation, name = default_rxn_name(ip), flux = NULL, ..., abscissa = NULL) {
  add_reaction_(ip, deparse(substitute(equation)), name = name, flux = lazy(flux),
                isotopes = lazy_dots(...), class = "custom", abscissa = abscissa)
}

#' Add a standard single fractionation reaction
#'
#' This adds a standard 1 to 1 (one reactant, one product) reaction with standard
#' fractionation factors to the isopath (for all isotopes that are part of both
#' reactant AND product). Supports reversibility for easy implementation of bi-
#' directional fluxes. This is really intended as a convenience function to make
#' it easy to set up Hayes-type reaction systems. Reaction systems that only have
#' this kind of reaction also have the added benefit of having an analytical
#' solution for steady state. See \link{calculate_steady_state()} to make use of
#' this functionalty.
#'
#' @param reversibility The reversibility of the reaction. If this is omitted, the reaction is treated as being completely irreversible if a kinetic fractionation factor (\code{alpha.<isotope>} or \code{eps.<isotope>}) is provided or as in full equilibrium if an equilibrium fractionation factor (\code{alpha.<isotope>.eq} or \code{eps.<isotope>.eq}) is provided. If anything else is provided (can be a number or expression; evaluation to 0 implies no reversibility, evaluation to 1 implies full reversibility), requires either a fractionation factor for the reverse reaction or an equilibrium fractionation factor to be defined in the dots (\code{...}) in addition to the forward kinetic effect for each isotope system involved in the reaction.
#' @param ... Fractionation factors for the reaction - can be numbers or expressions (referencing other variables and parameters in the system), naming convention: \code{alpha.<isotope> = ...}, or \code{eps.<isotope> = ...} for the kinetic isotope effect of the reaction in alpha or epsilon (permil) notation (with definition alpha = k_light / k_heavy and eps = (alpha - 1) * 1000, i.e. normal isotope effects are alpha > 1 and eps > 0, inverse isotope effects are alpha < 1 and eps < 0), and \code{alpha.<isotope>eq ...} or \code{eps.<isotope> = ...} for equilibrium isotope effects. If only kinetic fractionation factors are provided the reaction is created as an irreversible reaction, if only equilibrium fractionation factors are provided, it's created as an equilibrium reaction. For anything in between, use the \code{reversibility} parameter, If reversibility is set, a fractionation factor for the reverse reaction must be provided either directly (\code{alpha.<isotope>.rev = } or \code{eps.<isotope>.rev = ...}) or in the form of the equilibrium fractionation factor (\code{alpha.<isotope>.eq = ...} or \code{eps.<isotope>.eq = ...}). All equilibrium fractionatino factors are interpreted as product/substrate or substrate/product depending on the \code{eq_ratio} parameter..
#' @param permil Whether the delta and epsilon values are in permil notation (i.e. all multiplied by 1000) or not. All must be in the same notation.
#' @param eq_ratio eq_ratio Whether equilibrium fractionation factors are provided as product / substrate ratios (the default) or substrate / product ratios.
#' @export
#' @inheritParams add_custom_reaction
#' @family reaction functions
add_standard_reaction <- function(ip, equation, name = default_rxn_name(ip), flux = NULL, reversibility, ...,
                                permil = TRUE, eq_ratio = c("P/S", "S/P"), abscissa = NULL){

  # equation
  eq <- deparse(substitute(equation))

  # reaction components
  components <-
    parse_reaction_equation(eq) %>%
    as.list() %>% as_tibble() %>%
    gather(component, comp_stoic)

  # make sure it's a 1 to 1 reaction
  rxn_stoic <- components %>% summarize(n_reactant = sum(comp_stoic < 0), n_product = sum(comp_stoic > 0))
  if (rxn_stoic$n_reactant != 1 || rxn_stoic$n_product != 1) {
    stop("standard reactions can only be 1 to 1 transformations: ", eq,
         " (reactants: ", rxn_stoic$n_reactant, ", products: ", rxn_stoic$n_product, ")", call. = FALSE)
  }

  # check for missing components
  if (length(missing_comp <- setdiff(components$component, names(ip$components))) > 0)
    stop("missing component definition(s), make sure to add this with add_component() first: ",
         missing_comp %>% paste(collapse = ", "), call. = FALSE)

  # reaction isotopes summary
  rxn_isotopes <-
    right_join(
      ip %>% get_component_matrix() %>% gather(isotope, iso_stoic, -component, -variable),
      components, by = "component") %>%
    filter(!is.na(iso_stoic)) %>%
    group_by(isotope) %>%
    mutate(n_reactant = sum(comp_stoic < 0),
           n_product = sum(comp_stoic > 0)) %>%
    filter(n_reactant == 1 & n_product == 1) %>%
    select(isotope, component, comp_stoic) %>%
    spread(comp_stoic, component)

  # fractionation factors from dots
  ff_dots <- lazy_dots(...)

  # construct allowed pattern for the dots for the relevat isotopes
  if (missing(reversibility)) {
    # irreversible or equilibrium reaction
    ff_dots_pattern <- paste0("(alpha|eps).", rxn_isotopes$isotope, "(.eq)?") %>% setNames(rxn_isotopes$isotope)
  } else  {
    # reversible reaction
    ff_dots_pattern <- c(
      paste0("(alpha|eps).", rxn_isotopes$isotope) %>% setNames(rxn_isotopes$isotope),
      paste0("(alpha|eps).", rxn_isotopes$isotope, ".(rev|eq)") %>% setNames(paste0(rxn_isotopes$isotope, ".rev")))
  }

  # find matching dot names for each isotope
  ff_dots_names <- sapply(ff_dots_pattern, function(i){
    ff_name <- gsub(".", "\\.", paste0("^", i, "$"), fixed = T) %>% grep(names(ff_dots), value = TRUE)
    if (length(ff_name) == 0) return (NA)
    else if (length(ff_name) > 1) stop("shouldn't be possible to happen")
    return (ff_name)
  })

  # stop if any of the required ff dots are missing
  ff_dots_missing <- ff_dots_pattern[is.na(ff_dots_names)]
  if (length(ff_dots_missing) > 0) {
    stop("not all required fractionation factors available, missing: ",
         ff_dots_missing %>% paste(collapse = ", "), call. = FALSE)
  }

  # stop if equilibrum ff given but equilibrium_ratio not provided
  if (any(grepl("eq", ff_dots_names)) && (missing(eq_ratio) || !eq_ratio %in% c("P/S", "S/P"))) {
    stop("equilibrium fractionation factor(s) provided but unclear whether defined as product/substrate or substrate/product. Please specify parameter eq_ratio='P/S' or eq_ratio='S/P'.", call. = FALSE)
  } else if (!any(grepl("eq", ff_dots_names)) && !missing(eq_ratio)) {
    warning("The parameter `eq_ratio` was defined but no equilibrium fractionation factor(s) provided. Did you specify `eq_ratio` by accident or did you mean to supply .eq frationation factors but did not?", call. = F)
  }

  # determine type for the reaction
  if (missing(reversibility)) {
    if (all(grepl("eq", ff_dots_names))) type <- "EQ" # equilibrium
    else if (!any(grepl("eq", ff_dots_names))) type <- "IR" # irreverisble
    else stop("ambiguous reaction type (irreversible vs. equilibrium), mixed kinetic and equilibrium fractionation factors: ",
              ff_dots_names %>% paste(collapse = ", "), call. = F)
  } else
    type <- "RV" # reversible

  # construct flux calls
  if (type == "IR") {
    # irreversible
    mass_flux <- lazy(flux)
  } else if (type == "EQ") {
    # equilibrium
    mass_flux <- interp(lazy(dir_flux(flux_, rev = 1, dir = "+")), flux_ = lazy(flux))
    mass_flux_rev <- interp(lazy(dir_flux(flux_, rev = 1, dir = "-")), flux_ = lazy(flux))
  } else if (type == "RV") {
    # reversible reaction
    mass_flux <-
      interp(lazy(dir_flux(flux_, rev = rev_, dir = "+")),
             flux_ = lazy(flux), rev_ = substitute(reversibility))
    mass_flux_rev <-
      interp(lazy(dir_flux(flux_, rev = rev_, dir = "-")),
             flux_ = lazy(flux), rev_ = substitute(reversibility))
  }

  # construct bare flux isotope calls (simplify by only including parameters when
  # they do not equate to the default of the fractionate function)
  permil_text <- if(!permil) ", p = FALSE" else ""
  multi_text <- if(type == "EQ" && eq_ratio == "P/S") ", m = TRUE" else ""
  alpha_text <- sprintf("fractionate(d_, a = ff_%s%s)", permil_text, multi_text)
  eps_text <- sprintf("fractionate(d_, eps = ff_%s%s)", permil_text, multi_text)
  frac_alpha <- interp(lazy(x), x = parse(text = alpha_text, keep.source = F, n = NULL)[[1]])
  frac_eps <- interp(lazy(x), x = parse(text = eps_text, keep.source = F, n = NULL)[[1]])

  flux_isotope <- list()
  flux_isotope_rev <- list()

  # construct calls for each isotope
  for (i in 1:nrow(rxn_isotopes)) {

    iso <- rxn_isotopes[i,]$isotope
    reactant <- rxn_isotopes[i,]$`-1`
    product <- rxn_isotopes[i,]$`1`

    prefix <- gsub(ff_dots_pattern[iso], "\\1", ff_dots_names[iso])
    func <- if (prefix == "eps") frac_eps else frac_alpha
    ff_fwd <- ff_dots[[ff_dots_names[iso]]]
    flux_isotope[[paste0("flux.", iso)]] <-
      interp(func,d_ = as.name(paste0(reactant, ".", iso)), ff_ = ff_fwd)

    if (type == "EQ") {
      # in equilibrium scenario, the reverse flux is not fractionated (all
      # fractionation occurs on the forward flux)
      flux_isotope_rev[[paste0("flux.", iso)]] <-
        interp(lazy(d_), d_ = as.name(paste0(product, ".", iso)))
    } else if (type == "RV") {
      # reversible (more complex retrun flux)
      iso_rev <- paste0(iso, ".rev")
      prefix_rev <- gsub(ff_dots_pattern[iso_rev], "\\1", ff_dots_names[iso_rev])

      # figure out if equilibrium factor needs to be converted to forward
      suffix <- gsub(ff_dots_pattern[iso_rev], "\\2", ff_dots_names[iso_rev])
      if (suffix == "rev") {
        func <- if (prefix_rev == "eps") frac_eps else frac_alpha
        ff <- ff_dots[[ff_dots_names[iso_rev]]]
      } else if (suffix == "eq") {
        ff_fwd_alpha <- ff_fwd
        if (prefix == "eps") {
          # have to convert from eps to alpha
          ff_fwd_alpha <- interp(lazy(m_ * eps_ + 1), m_ = if (permil) 0.001 else 1, eps_ = ff_fwd_alpha)
        }
        ff_eq_alpha <- ff_dots[[ff_dots_names[iso_rev]]]
        if (prefix_rev == "eps") {
          # have to convert from eps to alpha
          ff_eq_alpha <- interp(lazy(m_ * eps_ + 1), m_ = if (permil) 0.001 else 1, eps_ = ff_eq_alpha)
        }

        # overall
        func <- frac_alpha
        if (eq_ratio == "S/P")
          ff <- interp(lazy(fwd/eq), fwd = ff_fwd_alpha, eq = ff_eq_alpha)
        else
          ff <- interp(lazy(fwd*eq), fwd = ff_fwd_alpha, eq = ff_eq_alpha)
      } else stop("should never happen")

      # add to flux list
      flux_isotope_rev[[paste0("flux.", iso)]] <-
        interp(func,d_ = as.name(paste0(product, ".", iso)), ff_ = ff)
    }

  }

  # add forward reaction
  ip <- add_reaction_(
    ip, eq, name = if(type != "IR") paste(name,"(forward)") else name,
    flux = mass_flux, isotopes = flux_isotope, class = "standard", abscissa = abscissa)

  if (type != "IR") {
    ip <- add_reaction_(
      ip, eq, name = paste(name,"(reverse)"),
      flux = mass_flux_rev, isotopes = flux_isotope_rev, class = "standard", abscissa = abscissa)
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
add_reaction_ <- function(ip, equation, name, flux, isotopes = list(),
                          abscissa = NULL, class = "generic", args = list()) {

  if (!is(ip, "isopath")) stop ("reaction can only be added to an isopath", call. = FALSE)

  # default abscissa
  components <- parse_reaction_equation(equation)
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
    parts <- strsplit(iso, ".", fixed = TRUE)[[1]]
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
      ip %>% get_reaction_component_matrix(),
      new_components %>% as.list() %>% as_tibble() %>% gather(new_component, new_comp_stoic),
      by = c("component" = "new_component")
    ) %>%
    select(abscissa, comp_stoic, new_comp_stoic) %>%
    distinct() %>%
    filter(!is.na(new_comp_stoic)) %>%
    filter(abscissa == min(abscissa)) %>%
    filter(comp_stoic == max(comp_stoic)) %>%
    arrange(new_comp_stoic)

  # components not found
  if (nrow(rxns) == 0)
    return( max(ip$get_reaction_component_matrix()$abscissa) + 1 )

  # figure out where to position reaction from existing component abscissa
  new_ab <- with(rxns[1,], {
    if (new_comp_stoic < 0) return (abscissa + 1)
    else if (new_comp_stoic > 0) return (abscissa)
    else stop ("should never happen")
  })
  return(new_ab)
}

#' parse reaction component
#' pulls out the stoichiometry and component
#' @param x the stoichiometric component, e.g. "A", "2 * A", etc.
#' @note standard evaluation
parse_reaction_component <- function(x) {
  m <- regexec("^\\s*(\\d+\\.?\\d*)?( \\* )?(\\w+)\\s*$", x)
  parts <- regmatches(x, m)[[1]]
  if ( length(parts) == 0)
    stop("cannot parse component: ", x, call. = FALSE)
  value <- if (parts[2] == "") 1 else as.numeric(parts[2])
  return(value %>% setNames(parts[4]))
}

#' parse reaction equation into its components
#' @note standard evaluation
parse_reaction_equation <- function(eq) {
  sides <- strsplit(eq, "==", fixed = TRUE)[[1]]
  err <- paste("please write equation in format 'a * A + b * B + ... == x * X + y * Y + ...':", eq)
  if (length(sides) != 2) stop(err, call. = FALSE) # missing two sides
  left <- strsplit(sides[1], "+", fixed = TRUE)[[1]]
  right <- strsplit(sides[2], "+", fixed = TRUE)[[1]]
  if (length(left) == 0 || length(right) == 0) stop(err, call. = FALSE)
  c(
    lapply(left, parse_reaction_component) %>% unlist() * -1,
    lapply(right, parse_reaction_component) %>% unlist() * 1
  )
}
