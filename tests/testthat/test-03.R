context("Steady-state")

test_that("Testing simple reactions steady state", {
  sys <- isopath() %>%
    add_isotope(c("carbon", "nitrogen")) %>%
    add_component("A", carbon, nitrogen) %>%
    add_component("B", carbon, nitrogen) %>%
    add_component("C", carbon)
  expect_error( sys %>% add_simple_reaction(A == D), "missing component definition")
  expect_error( sys %>% add_simple_reaction(A == B + C), "simple reactions can only be 1 to 1 transformations")

  # missing fractionatoin factor for carbon systen
  expect_error( sys %>% add_simple_reaction(A == C), "not all required fractionation factors available")
  expect_is( sys %>% add_simple_reaction(A == C, alpha.carbon = 1), "isopath")

  # now both carbon and nitrogen are involved in the reaction, missing nitrogen
  expect_error( sys %>% add_simple_reaction(A == B, alpha.carbon = 1), "not all required fractionation factors available")
  expect_is( sys %>% add_simple_reaction(A == B, alpha.carbon = 1, eps.nitrogen = 1), "isopath")

  # once reversability is activated, need to actually supply the reverse too, missing rev. carbon or rev. nitrogen
  expect_error( sys %>% add_simple_reaction(A == C, eps.carbon = 1, rev = 0.5),
                "not all required fractionation factors available")
  expect_error( sys %>% add_simple_reaction(A == B, eps.carbon = 1, alpha.carbon.re = 1, nitrogen.alpha = 0.2, rev = 0.5),
                "not all required fractionation factors available")
  expect_is(sys %>% add_simple_reaction(A == C, eps.carbon = 1, alpha.carbon.rev = x, rev = 0.5), "isopath")

  # test irreversible system
  expect_equal(sys %>% add_simple_reaction(A == B, alpha.carbon = cff, eps.nitrogen = nff, flux = my_flux) %>%
                 get_flux_matrix(),
               data_frame(reaction = "rxn1", flux = "my_flux"))
  expect_equal(sys %>% add_simple_reaction(A == B, alpha.carbon = cff, eps.nitrogen = nff, flux = my_flux) %>%
                 get_flux_isotope_matrix(),
               data_frame(
                 reaction = "rxn1",
                 isotope = rep(c("carbon", "nitrogen"), each=2),
                 component = rep(c("A", "B"), times=2),
                 flux_isotope = rep(c("fractionate(A.carbon, a = cff)",
                                      "fractionate(A.nitrogen, eps = nff)"), each = 2))
  )

  # test reversible system
  expect_equal(sys %>% add_simple_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.rev = cff.rev,
                                           flux = my_flux, reversibility = my_rev) %>%
                 get_flux_matrix(),
               data_frame(reaction = c("rxn1 (forward)", "rxn1 (reverse)"),
                          flux = c("flux(my_flux, rev = my_rev, dir = \"+\")", "flux(my_flux, rev = my_rev, dir = \"-\")")))
  expect_equal(sys %>% add_simple_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.rev = cff.rev,
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
  expect_equal((sys %>% add_simple_reaction(A == C, alpha.carbon = cff.fwd, alpha.carbon.eq = cff.eq,
                                           flux = my_flux, reversibility = my_rev) %>%
                 get_flux_isotope_matrix())$flux_isotope[4],
               "fractionate(C.carbon, a = cff.fwd/cff.eq)")
  expect_equal((sys %>% add_simple_reaction(A == C, alpha.carbon = cff.fwd, eps.carbon.eq = eff.eq,
                                            flux = my_flux, reversibility = my_rev) %>%
                  get_flux_isotope_matrix())$flux_isotope[4],
               "fractionate(C.carbon, a = cff.fwd/(0.001 * eff.eq + 1))")
})
