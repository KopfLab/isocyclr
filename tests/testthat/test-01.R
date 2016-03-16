context("Struture")

test_that("Generating isopath works", {
  expect_is(isopath(), "isopath")
})

test_that("Adding isotope works", {
  expect_error(add_isotope("not correct"), "can only be added to an isopath")
  sys <- isopath()
  expect_error( add_isotope(sys, "!@#$"), "only alphanumeric")
  expect_equal( add_isotope(sys, "C")$isotopes %>% names(), "C")
})

test_that("Adding components works", {
  expect_error(add_component("not correct"), "can only be added to an isopath")
  sys <- isopath() %>% add_isotope("N")
  expect_error( add_component(sys, "!@#$"), "only alphanumeric")
  expect_equal( add_component(sys, "A")$components %>% names, "A")
  expect_error( add_component(sys, "A", x * N), "cannot parse")
  expect_error( add_component(sys, "A", `%$#`), "cannot parse")
  expect_equal( add_component(sys, "A", N)$components$A$isotopes, c(N = 1) )
  expect_equal( add_component(sys, "A", 2 * N)$components$A$isotopes, c(N = 2) )
  expect_error( add_component(sys, "A", 2 * N, C), "missing isotope definition" )
  sys <- sys %>% add_isotope("C")
  expect_equal( add_component(sys, "A", 2 * N, C)$components$A$isotopes, c(N = 2, C = 1) )
})

test_that("Adding reaction equations works", {
  expect_error(add_reaction("not correct"), "can only be added to an isopath")
  sys <- isopath() %>%
    add_isotope("N") %>%
    add_component("A", N) %>%
    add_component("B", N)
  expect_error( add_reaction(sys, "rxn1", A), "please write equation in format")
  expect_equal( add_reaction(sys, "rxn1", A == B)$reactions %>% names, "rxn1")
  expect_error( add_reaction(sys, "rxn1", A == C), "missing component definition")
  sys <- sys %>% add_component("C") %>% add_component("D")
  expect_equal( add_reaction(sys, "rxn1", A + 5*B == C + 2*D)$reactions$rxn1$components,
                c(A = 1, B = 5, C = -1, D = -2))
})

test_that("Isopath structure matrices work", {
  expect_error(get_reaction_matrix("not correct"), "can only get .* from an isopath")
  expect_error(get_component_matrix("not correct"), "can only get .* from an isopath")

  sys <- isopath() %>%
    add_isotope("C") %>%
    add_isotope("N") %>%
    add_component("X", 2 * C, N) %>%
    add_component("Y", C, variable = FALSE) %>%
    add_component("Z") %>%
    add_component("W") %>%
    add_reaction("rxn1", X == 3 * Y) %>%
    add_reaction("rxn2", Y + 2 * Z == W)

  expect_equal(sys %>% get_component_matrix(),
    data_frame(component = c("X", "Y", "Z", "W"),
               variable = c(T, F, T, T),
               C = c(2, 1, NA, NA),
               N = c(1, NA, NA, NA),
               unspecified = c(NA, NA, 1, 1))
  )

  expect_equal(sys %>% get_reaction_matrix(),
               data_frame(reaction = c("rxn1", "rxn2"),
                          rxn_nr = c(1,2),
                          X = c(1, NA), Y = c(-3, 1),
                          Z = c(NA, 2), W = c(NA, -1))
  )

  expect_equal(sys %>% get_reaction_component_matrix(),
               data_frame(
                 reaction = rep(c("rxn1", "rxn2"), each = 3),
                 rxn_nr = rep(c(1, 2), each = 3),
                 component = c("X", "X", "Y", "Y", "Z", "W"),
                 comp_stoic = c(1, 1, -3, 1, 2, -1),
                 variable = c(T, T, F, F, T, T),
                 isotope = c("C", "N", "C", "C", "unspecified", "unspecified"),
                 iso_stoic = c(2, 1, 1, 1, 1, 1)
               ))

})

test_that("Adding reaction flux and isotopes works", {
  sys <- isopath() %>%
    add_isotope("C") %>%
    add_component("X", C) %>%
    add_component("Y", C)
  expect_error( add_reaction(sys, "rxn1", X == Y, C = 1), "missing prefix for isotope flux")
  expect_error( add_reaction(sys, "rxn1", X == Y, flux.N = 1), "missing isotope definition")
  expect_error( add_reaction(sys, "rxn1", X == Y, flux.C.abc = 1), "missing component definition")
  expect_equal( add_reaction(sys, "rxn1", X == Y, flux.C = 1)$reactions$rxn1$isotopes$C$expr, 1)
  expect_equal( add_reaction(sys, "rxn1", X == Y, flux.C = X + Y)$reactions$rxn1$isotopes$C$expr %>% deparse(), "X + Y")
  expect_equal( add_reaction(sys, "rxn1", X == Y, flux.X.C = 5, flux.Y.C = 3)$reactions$rxn1$isotopes$C %>%
                  names(), c("X", "Y"))
  expect_equal( add_reaction(sys, "rxn1", X == Y, flux.X.C = pi(), flux.Y.C = 3)$reactions$rxn1$isotopes$C$X$expr %>%
                  deparse(), "pi()")
  expect_equal( get_variables(sys), c("X", "X.C", "Y", "Y.C") )
  expect_equal( sys %>% add_component("X", C, variable = F) %>% get_variables(), c("Y", "Y.C") )
  expect_error( set_parameters(sys, data_frame(a = 5)), "parameters required for minimal parameter set missing")
  expect_equal( set_parameters(sys, data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4))$parameters,
                data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4))
})
