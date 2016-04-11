context("Standard reactions")

test_that("Testing standard reactions", {
  sys <- isopath() %>%
    add_isotope(c("carbon", "nitrogen")) %>%
    add_component("A", carbon, nitrogen) %>%
    add_component("B", carbon, nitrogen) %>%
    add_component("C", carbon)
  expect_error( sys %>% add_standard_reaction(A == D), "missing component definition")
  expect_error( sys %>% add_standard_reaction(A == B + C), "standard reactions can only be 1 to 1 transformations")

  # missing fractionatoin factor for carbon systen
  expect_error( sys %>% add_standard_reaction(A == C), "not all required fractionation factors available")
  expect_is( sys %>% add_standard_reaction(A == C, alpha.carbon = 1), "isopath")

  # missing equilibrium definitions
  expect_error(sys %>% add_standard_reaction(A == C, alpha.carbon.eq = 1), "equilibrium .* unclear")
  expect_error(sys %>% add_standard_reaction(A == C, alpha.carbon.eq = 1, eq_ratio = "err"), "equilibrium .* unclear")
  expect_is( sys %>% add_standard_reaction(A == C, alpha.carbon.eq = 1, eq_ratio = "P/S"), "isopath")

  # now both carbon and nitrogen are involved in the reaction, missing nitrogen
  expect_error( sys %>% add_standard_reaction(A == B, alpha.carbon = 1), "not all required fractionation factors available")
  expect_is( sys %>% add_standard_reaction(A == B, alpha.carbon = 1, eps.nitrogen = 1), "isopath")

  # ambiguous reaction type with carbon ff supplied as kinetic and nitrogen ff as equilibrium
  expect_error(sys %>% add_standard_reaction(A == B, alpha.carbon = 1, alpha.nitrogen.eq = 1, eq_ratio = "P/S"),
               "ambiguous reaction type")
  expect_warning(sys %>% add_standard_reaction(A == B, alpha.carbon = 1, alpha.nitrogen = 1, eq_ratio = "P/S"),
                                             "defined but no equilibrium fractionation factor")
  expect_is( sys %>% add_standard_reaction(A == B, alpha.carbon.eq = 1, eps.nitrogen.eq = 1, eq_ratio = "P/S"), "isopath")

  # once reversability is activated, need to actually supply the reverse too, missing rev. carbon or rev. nitrogen
  expect_error( sys %>% add_standard_reaction(A == C, eps.carbon = 1, rev = 0.5),
                "not all required fractionation factors available")
  expect_error( sys %>% add_standard_reaction(A == B, eps.carbon = 1, alpha.carbon.re = 1, nitrogen.alpha = 0.2, rev = 0.5),
                "not all required fractionation factors available")
  expect_is(sys %>% add_standard_reaction(A == C, eps.carbon = 1, alpha.carbon.rev = x, rev = 0.5), "isopath")

  # test irreversible system
  expect_equal(sys %>% add_standard_reaction(A == B, alpha.carbon = cff, eps.nitrogen = nff, flux = my_flux) %>%
                 get_ode_matrix(),
               data_frame(
                 x = c("A", "A.carbon", "A.nitrogen", "B", "B.carbon", "B.nitrogen"),
                 value = c("A", "A.carbon", "A.nitrogen", "B", "B.carbon", "B.nitrogen"),
                 `dx/dt` = c("-1 * my_flux",
                             "-1 * my_flux/A * (fractionate(A.carbon, a = cff) - A.carbon)",
                             "-1 * my_flux/A * (fractionate(A.nitrogen, eps = nff) - A.nitrogen)",
                             "1 * my_flux",
                             "1 * my_flux/B * (fractionate(A.carbon, a = cff) - B.carbon)",
                             "1 * my_flux/B * (fractionate(A.nitrogen, eps = nff) - B.nitrogen)")
               ))

  # test equilibrium system
  expect_equal(sys %>% add_standard_reaction(A == C, alpha.carbon.eq = ceq, flux = my_flux, eq_ratio = "S/P") %>%
                 get_ode_matrix(),
               data_frame(
                 x = c("A", "A.carbon", "A.nitrogen", "C", "C.carbon"),
                 value = c("A", "A.carbon", "A.nitrogen", "C", "C.carbon"),
                 `dx/dt` = c(
                   "-1 * dir_flux(my_flux, rev = 1, dir = \"+\") + -1 * dir_flux(my_flux, rev = 1, dir = \"-\")",
                   "-1 * dir_flux(my_flux, rev = 1, dir = \"+\")/A * (fractionate(A.carbon, a = ceq) - A.carbon) + -1 * dir_flux(my_flux, rev = 1, dir = \"-\")/A * (C.carbon - A.carbon)",
                   "-1 * dir_flux(my_flux, rev = 1, dir = \"+\")/A * (A.nitrogen - A.nitrogen) + -1 * dir_flux(my_flux, rev = 1, dir = \"-\")/A * (A.nitrogen - A.nitrogen)",
                   "1 * dir_flux(my_flux, rev = 1, dir = \"+\") + 1 * dir_flux(my_flux, rev = 1, dir = \"-\")",
                   "1 * dir_flux(my_flux, rev = 1, dir = \"+\")/C * (fractionate(A.carbon, a = ceq) - C.carbon) + 1 * dir_flux(my_flux, rev = 1, dir = \"-\")/C * (C.carbon - C.carbon)")))
  expect_equal((sys %>% add_standard_reaction(A == C, alpha.carbon.eq = ceq, flux = my_flux, eq_ratio = "P/S") %>%
                  get_ode_matrix())$`dx/dt`[2],
               "-1 * dir_flux(my_flux, rev = 1, dir = \"+\")/A * (fractionate(A.carbon, a = ceq, m = TRUE) - A.carbon) + -1 * dir_flux(my_flux, rev = 1, dir = \"-\")/A * (C.carbon - A.carbon)")


  # test reversible system
  expect_equal(sys %>% add_standard_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.rev = cff.rev,
                                           flux = my_flux, reversibility = my_rev) %>%
                 get_flux_matrix(),
               data_frame(reaction = c("rxn1 (forward)", "rxn1 (reverse)"),
                          flux = c("dir_flux(my_flux, rev = my_rev, dir = \"+\")", "dir_flux(my_flux, rev = my_rev, dir = \"-\")")))
  expect_equal(sys %>% add_standard_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.rev = cff.rev,
                                           flux = my_flux, reversibility = my_rev) %>%
                 get_flux_isotope_matrix(),
               data_frame(
                 reaction = rep(c("rxn1 (forward)", "rxn1 (reverse)"), each = 2),
                 isotope = "carbon",
                 component = rep(c("A", "C"), times=2),
                 flux_isotope = rep(c("fractionate(A.carbon, a = cff.fwd)",
                                      "fractionate(C.carbon, a = cff.rev)"), each = 2)
               ))

  # test reversible that uses equilibrium isotope effect
  expect_equal((sys %>% add_standard_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.eq = cff.eq, eq_ratio = "S/P",
                                           flux = my_flux, reversibility = my_rev) %>%
                 get_flux_isotope_matrix())$flux_isotope[4],
               "fractionate(C.carbon, a = cff.fwd/cff.eq)")
  expect_equal((sys %>% add_standard_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.eq = cff.eq, eq_ratio = "P/S",
                                            flux = my_flux, reversibility = my_rev) %>%
                  get_flux_isotope_matrix())$flux_isotope[4],
               "fractionate(C.carbon, a = cff.fwd * cff.eq)")
  expect_equal((sys %>% add_standard_reaction(A == C, alpha.carbon = cff.fwd, eps.carbon.eq = eff.eq, eq_ratio = "S/P",
                                            flux = my_flux, reversibility = my_rev) %>%
                  get_flux_isotope_matrix())$flux_isotope[4],
               "fractionate(C.carbon, a = cff.fwd/(0.001 * eff.eq + 1))")
  expect_equal((sys %>% add_standard_reaction(A == C, alpha.carbon = cff.fwd, eps.carbon.eq = eff.eq, eq_ratio = "P/S",
                                            flux = my_flux, reversibility = my_rev) %>%
                  get_flux_isotope_matrix())$flux_isotope[4],
               "fractionate(C.carbon, a = cff.fwd * (0.001 * eff.eq + 1))")
})
