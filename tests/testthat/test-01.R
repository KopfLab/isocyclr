context("Struture")

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
  expect_error( add_custom_reaction(sys, "rxn1", A), "please write equation in format")
  expect_equal( add_custom_reaction(sys, "rxn1", A == B)$reactions %>% names(), "rxn1")
  expect_equal( add_custom_reaction(sys, A == B)$reactions %>% names(), "rxn1")
  expect_error( add_custom_reaction(sys, "rxn1", A == C), "missing component definition")
  sys <- sys %>% add_component("C") %>% add_component("D")
  expect_equal( add_custom_reaction(sys, "rxn1", A + 5*B == C + 2*D)$reactions$rxn1$components,
                c(A = -1, B = -5, C = 1, D = 2))
})

test_that("System info cache works", {
  sys <- isopath() %>%
    add_isotope("N") %>%
    add_component("A", N) %>%
    add_component("B", N, variable = FALSE)
  expect_equal(sys$info$reaction_component_matrix, sys %>% get_reaction_component_matrix())
  expect_equal(sys$info$variable_reaction_component_matrix, sys %>% get_reaction_component_matrix() %>% filter(variable == T))
  expect_equal(sys$info$variables, sys %>% get_variables())
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
    add_custom_reaction("rxn1", X == 3 * Y) %>%
    add_custom_reaction("rxn2", Y + 2 * Z == W)

  expect_equal(sys %>% get_component_matrix(),
    data_frame(component = c("X", "Y", "Z", "W"),
               variable = c(T, F, T, T),
               C = c(2, 1, NA, NA),
               N = c(1, NA, NA, NA),
               unspecified = c(NA, NA, 1, 1))
  )

  expect_equal(sys %>% get_reaction_matrix(),
               data_frame(reaction = c("rxn1", "rxn2"),
                          abscissa = c(1,2),
                          X = c(-1, NA), Y = c(3, -1),
                          Z = c(NA, -2), W = c(NA, 1))
  )

  expect_equal(sys %>% get_reaction_component_matrix(),
               data_frame(
                 reaction = paste0("rxn", c(1, 1, 2, 1, 2, 2)),
                 abscissa = c(0, 1, 1, 0, 1, 2),
                 component = c("X", "Y", "Y", "X", "Z", "W"),
                 comp_stoic = c(-1, 3, -1, -1, -2, 1),
                 variable = c(T, F, F, T, T, T),
                 isotope = c("C", "C", "C", "N", "unspecified", "unspecified")
               ))


  expect_equal(
    isopath() %>%
      add_component(LETTERS[1:6]) %>%
      add_custom_reaction(A == D) %>%
      add_custom_reaction(A == C) %>% # old reactant = new reactant
      add_custom_reaction(D == F) %>% # old reactant = new product
      add_custom_reaction(E == A) %>% # new reactant = old product
      add_custom_reaction(B == D) %>% # new product  = old product
      get_reaction_component_matrix() %>%
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
  expect_error( add_custom_reaction(sys, "rxn1", X == Y, C = 1), "missing prefix for isotope flux")
  expect_error( add_custom_reaction(sys, "rxn1", X == Y, flux.N = 1), "missing isotope definition")
  expect_error( add_custom_reaction(sys, "rxn1", X == Y, flux.C.abc = 1), "missing component definition")
  expect_equal( add_custom_reaction(sys, "rxn1", X == Y, flux.C = 1)$reactions$rxn1$isotopes$C$expr, 1)
  expect_equal( add_custom_reaction(sys, "rxn1", X == Y, flux.C = X + Y)$reactions$rxn1$isotopes$C$expr %>% deparse(), "X + Y")
  expect_equal( add_custom_reaction(sys, "rxn1", X == Y, flux.X.C = 5, flux.Y.C = 3)$reactions$rxn1$isotopes$C %>%
                  names(), c("X", "Y"))
  expect_equal( add_custom_reaction(sys, "rxn1", X == Y, flux.X.C = pi(), flux.Y.C = 3)$reactions$rxn1$isotopes$C$X$expr %>%
                  deparse(), "pi()")
})

test_that("Adding parameters works", {
  sys <- isopath() %>%
    add_isotope("C") %>%
    add_component("X", C) %>%
    add_component("Y", C)
  expect_equal( get_variables(sys), c("X", "X.C", "Y", "Y.C") )
  expect_equal( sys %>% add_component("X", C, variable = F) %>% get_variables(), c("Y", "Y.C") )
  expect_equal( {sys2 <- set_parameters(sys, data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4)); sys2$parameters},
                data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4))
  expect_equal( set_parameters(sys, X = 1, X.C = 2, Y = 3, Y.C = 4)$parameters,
                data_frame(X = 1, X.C = 2, Y = 3, Y.C = 4))
  expect_equal( set_parameters(sys2, X = 100, new = 0.1)$parameters,
                data_frame(X = 100, X.C = 2, Y = 3, Y.C = 4, new = 0.1))
  # multi line
  sys3 <- set_parameters(sys2, X = 100, new = 0.1)
  expect_equal(
    set_parameters(sys3, Z = c(1,2))$parameters,
    data_frame(X = c(100, 100), X.C = c(2, 2), Y = c(3, 3), Y.C = c(4, 4),
               new = c(0.1, 0.1), Z = c(1, 2)))
  expect_equal(
    (sys3 %>% set_parameters(Z = c(1,2)) %>% set_parameters(X = 5))$parameters,
    data_frame(X = c(5, 5), X.C = c(2, 2), Y = c(3, 3), Y.C = c(4, 4),
               new = c(0.1, 0.1), Z = c(1, 2)))

})

test_that("Evaluation works", {
  expect_error(get_flux_matrix(NULL), "can only get flux matrix from an isopath")
  expect_error(get_flux_isotope_matrix(NULL), "can only get flux isotope matrix from an isopath")
  expect_error(get_component_change_summary(NULL), "can only get component change for an isopath")
  expect_error(get_isotope_change_summary(NULL), "can only get isotope change for an isopath")

  sys <- isopath() %>%
    add_isotope("C") %>% add_isotope("N") %>%
    add_component("X", C, N) %>% add_component("Y", C, N) %>%
    add_custom_reaction("rxn1", X == 2 * Y, flux = dm, flux.N = dN, flux.X.C = X.dC, flux.Y.C = Y.dC) %>%
    set_parameters(X = 1, X.C = 1, X.N = 1, Y = 1, Y.C = 1, Y.N = 1)

  # flux matrix
  expect_equal(sys  %>% get_flux_matrix(), data_frame(reaction = "rxn1", flux = "dm"))
  expect_error(sys  %>% get_flux_matrix(eval = T))
  expect_equal(sys  %>% get_flux_matrix(eval = T, param = data_frame(dm = 3)), data_frame(reaction = "rxn1", flux = 3))
  expect_equal(sys %>% set_parameters(dm = 3) %>%
                 get_flux_matrix(eval = T), data_frame(reaction = "rxn1", flux = 3))

  # flux component summary
  expect_equal(sys %>% set_parameters(dm = 2) %>%
                 get_component_change_summary(),
               data_frame(component = c("X", "Y"), pool_size = c(1, 1), `dx/dt` = c(-2, 4)))

  # flux isotopes matrix
  expected <- data_frame(reaction = "rxn1", isotope = c("N", "N", "C", "C"),
                         component = c("X", "Y", "X", "Y"),
                         flux_isotope = c("dN", "dN", "X.dC", "Y.dC"))
  expect_equal(sys  %>% get_flux_isotope_matrix(), expected)
  expect_error(sys  %>% get_flux_isotope_matrix(eval = T))
  expect_equal(sys  %>%
                 get_flux_isotope_matrix(eval = T,
                                         param = data_frame(dN = 0.1, X.dC = 0.4, Y.dC = 0.6)),
               expected %>% mutate(flux_isotope = c(0.1, 0.1, 0.4, 0.6)))
  expect_equal(sys %>%
                 set_parameters(dN = 0.1, X.dC = 0.4, Y.dC = 0.6) %>%
                 get_flux_isotope_matrix(eval = T),
               expected %>% mutate(flux_isotope = c(0.1, 0.1, 0.4, 0.6)))

  expect_equal(sys  %>% set_parameters(dm = 2, dN = 0.1, X.dC = 0.4, Y.dC = 0.6) %>%
                 get_isotope_change_summary(),
               data_frame(isotope = c("N", "N", "C", "C"),
                          component = c("X", "Y", "X", "Y"),
                          pool_isotope = c(1, 1, 1, 1),
                          `dx/dt` = c(1.8, -3.6, 1.2, -1.6)))

})
