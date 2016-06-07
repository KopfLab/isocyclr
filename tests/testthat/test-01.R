context("Structure")

test_that("Generating isopath works", {
  expect_is(isopath(), "isopath")
})

test_that("Adding isotope works", {
  expect_error(add_isotope("not correct"), "can only be added to an isopath")
  sys <- isopath()
  expect_error( add_isotope(sys, "!@#$"), "only alphanumeric")
  expect_equal( add_isotope(sys, "C")$isotopes %>% names(), "C")
  expect_equal( add_isotope(sys, c("C", "N"))$isotopes %>% names(), c("C", "N"))
})

test_that("Adding components works", {
  expect_error(add_component("not correct", "A"), "can only be added to an isopath")
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
  expect_error(add_custom_reaction("not correct"), "can only be added to an isopath")
  sys <- isopath() %>%
    add_isotope("N") %>%
    add_component("A", N) %>%
    add_component("B", N)
  expect_error( add_custom_reaction(sys, A), "please write equation in format")
  expect_equal( add_custom_reaction(sys, A == B)$reactions %>% names(), "rxn1")
  expect_equal( add_custom_reaction(sys, A == B, "my_rxn")$reactions %>% names(), "my_rxn")
  expect_error( add_custom_reaction(sys, A == C), "missing component definition")
  sys <- sys %>% add_component("C") %>% add_component("D")
  expect_equal( add_custom_reaction(sys, A + 5*B == C + 2*D)$reactions$rxn1$components,
                c(A = -1, B = -5, C = 1, D = 2))
})

test_that("Isopath structure matrices work", {
  expect_error(get_reaction_matrix2("not correct"), "can only get .* from an isopath")
  expect_error(get_component_matrix("not correct"), "can only get .* from an isopath")

  sys <- isopath() %>%
    add_isotope("C") %>%
    add_isotope("N") %>%
    add_component("X", 2 * C, N) %>%
    add_component("Y", C, variable = FALSE) %>%
    add_component("Z") %>%
    add_component("W") %>%
    add_custom_reaction(X == 3 * Y) %>%
    add_custom_reaction(Y + 2 * Z == W)

  # component matrix
  expect_equal(sys %>% get_component_matrix(),
    data_frame(component = c("X", "Y", "Z", "W"),
               variable = c(T, F, T, T),
               C = c(2, 1, NA, NA),
               N = c(1, NA, NA, NA),
               unspecified = c(NA, NA, 1, 1))
  )

  # reaction matrix
  expect_equal(sys %>% get_reaction_matrix2(),
               data_frame(reaction = c("rxn1", "rxn2"),
                          abscissa = c(1,2),
                          flux = c("NULL", "NULL"),
                          X = c(-1, NA), Y = c(3, -1),
                          Z = c(NA, -2), W = c(NA, 1))
  )

  expect_equal(sys %>% get_reaction_component_matrix2(),
               data_frame(
                 component = c("W", "X", "Y", "Y", "Z"),
                 abscissa = c(2, 0, 1, 1, 1),
                 variable = c(TRUE, TRUE, FALSE, FALSE, TRUE),
                 reaction = c("rxn2", "rxn1", "rxn1", "rxn2", "rxn2"),
                 comp_stoic = c(1, -1, 3, -1, -2),
                 flux = c("NULL", "NULL", "NULL", "NULL", "NULL"),
                 pool_size = c("W", "X", "Y", "Y", "Z"),
                 `dx/dt` = c("1 * NULL", "-1 * NULL", "0", "0", "-2 * NULL")
               ))


  expect_equal(
    isopath() %>%
      add_component(LETTERS[1:6]) %>%
      add_custom_reaction(A == D) %>%
      add_custom_reaction(A == C) %>% # old reactant = new reactant
      add_custom_reaction(D == F) %>% # old reactant = new product
      add_custom_reaction(E == A) %>% # new reactant = old product
      add_custom_reaction(B == D) %>% # new product  = old product
      get_reaction_component_matrix2() %>%
      select(abscissa, component, comp_stoic),
    data_frame(
      abscissa = c(-1, 0, 0, 0, 0, 1, 1, 1, 1, 2),
      component = c("E", "A", "A", "A", "B", "D", "C", "D", "D", "F"),
      comp_stoic = c(-1, 1, -1, -1, -1, 1, 1, 1, -1, 1))

  )

})

test_that("Adding reaction flux and isotopes works", {
  sys <- isopath() %>%
    add_isotope("C") %>%
    add_component("X", C) %>%
    add_component("Y", C)
  expect_error( add_custom_reaction(sys, X == Y, C = 1), "missing prefix for isotope flux")
  expect_error( add_custom_reaction(sys, X == Y, flux.N = 1), "missing isotope definition")
  expect_error( add_custom_reaction(sys, X == Y, flux.C.abc = 1), "missing component definition")
  expect_equal( add_custom_reaction(sys, X == Y, flux.C = 1)$reactions$rxn1$isotopes$C$expr, 1)
  expect_equal( add_custom_reaction(sys, X == Y, flux.C = X + Y)$reactions$rxn1$isotopes$C$expr %>% deparse(), "X + Y")
  expect_equal( add_custom_reaction(sys, X == Y, flux.X.C = 5, flux.Y.C = 3)$reactions$rxn1$isotopes$C %>%
                  names(), c("X", "Y"))
  expect_equal( add_custom_reaction(sys, X == Y, flux.X.C = pi(), flux.Y.C = 3)$reactions$rxn1$isotopes$C$X$expr %>%
                  deparse(), "pi()")
})

test_that("Adding parameters works", {

  expect_error(set_parameters("incorrect"), "can only be .* for an isopath")

  sys <- isopath() %>%
    add_isotope("C") %>%
    add_component("X", C) %>%
    add_component("Y", C)
  expect_equal( {sys2 <- set_parameters(sys, data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4)); sys2$parameters},
                data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4))
  expect_equal( set_parameters(sys, X = 1, X.C = 2, Y = 3, Y.C = 4)$parameters,
                data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4))
  expect_equal( set_parameters(sys2, X = 100, new = 0.1)$parameters,
                data_frame(X = 100, X.C = 2, Y = 3, Y.C = 4, new = 0.1))

  # expand parameters
  expect_equal(
    expand_parameters(sys2, Z = c(1,2))$parameters,
    data_frame(X = c(1, 1), X.C = c(2, 2), Y = c(3, 3), Y.C = c(4, 4), Z = c(1, 2)))
  expect_equal(
    (sys2 %>% expand_parameters(Z = c(1,2)) %>% expand_parameters(new = 1:3))$parameters,
    data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4, Z = c(1, 1, 1, 2, 2, 2), new = c(1:3, 1:3)))

  # expanding parameters
  expect_error(expand_parameters("correct"), "can only be .* for an isopath")
  expect_error(expand_parameters(sys, x=1), "no parameters set yet for this isopath")


})

test_that("Evaluation works", {
  #expect_error(get_reaction_component_matrix2(NULL), "can only get flux matrix from an isopath")
  #expect_error(get_reaction_isotope_matrix2(NULL), "can only get flux matrix from an isopath")
  #expect_error(get_reaction_ode_matrix2(NULL), "can only get flux matrix from an isopath")

  sys <- isopath() %>%
    add_isotope("C") %>% add_isotope("N") %>%
    add_component("X", C, N) %>% add_component("Y", C, N) %>%
    add_custom_reaction(X == 2 * Y, flux = dm, flux.N = dN, flux.X.C = X.dC, flux.Y.C = Y.dC) %>%
    set_parameters(X = 1, X.C = 1, X.N = 1, Y = 1, Y.C = 1, Y.N = 1)

  # flux matrix
  # expect_equal(sys  %>% get_flux_matrix(), data_frame(reaction = "rxn1", flux = "dm"))
  # expect_equal(sys  %>% get_flux_matrix(eval = TRUE, param = data_frame(dm = 3)), data_frame(reaction = "rxn1", flux = 3))
  # expect_equal(sys %>% set_parameters(dm = 3) %>%
  #                get_flux_matrix(eval = TRUE), data_frame(reaction = "rxn1", flux = 3))

  # flux component summary
  # FIXME update documentation

  # flux isotopes matrix
  # expected <- data_frame(reaction = "rxn1", isotope = c("N", "N", "C", "C"),
  #                        component = c("X", "Y", "X", "Y"),
  #                        flux_isotope = c("dN", "dN", "X.dC", "Y.dC"))
  # expect_equal(sys  %>% get_flux_isotope_matrix(), expected)
  # expect_error(sys  %>% get_flux_isotope_matrix(eval = TRUE))
  # expect_equal(sys  %>%
  #                get_flux_isotope_matrix(eval = TRUE,
  #                                        param = data_frame(dN = 0.1, X.dC = 0.4, Y.dC = 0.6)),
  #              expected %>% mutate(flux_isotope = c(0.1, 0.1, 0.4, 0.6)))
  # expect_equal(sys %>%
  #                set_parameters(dN = 0.1, X.dC = 0.4, Y.dC = 0.6) %>%
  #                get_flux_isotope_matrix(eval = TRUE),
  #              expected %>% mutate(flux_isotope = c(0.1, 0.1, 0.4, 0.6)))
  #
  # expect_equal(sys  %>% set_parameters(dm = 2, dN = 0.1, X.dC = 0.4, Y.dC = 0.6) %>%
  #                get_ode_matrix(eval = T),
  #              data_frame(x = c("X", "X.C", "X.N", "Y", "Y.C", "Y.N"),
  #                         value = c(1, 1, 1, 1, 1, 1),
  #                         `dx/dt` = c(-2.0, 1.2, 1.8, 4.0, -1.6, -3.6)))

})
