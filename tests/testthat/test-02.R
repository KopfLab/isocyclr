context("Modelling")

test_that("Calculations work", {

  # test fractionate function
  expect_error(fractionate(5), "no fractionation factor provided")
  expect_error(fractionate(alpha = 0.99, eps = -10), "ambiguous fractionation factor")
  expect_equal(fractionate(d = 10, eps = -10, permil = TRUE), 1000*fractionate(d = 0.01, eps = -0.01, permil = FALSE))
  expect_equal(fractionate(d = 10, a = 0.99, mult = TRUE), fractionate(d = 10, a = 1/0.99, mult = FALSE))
  expect_equal(fractionate(d = 10, a = 0.99), fractionate(d = 10, eps = -10))
  expect_equal(fractionate(d = 0.01, a = 0.99, p = FALSE), fractionate(d = 0.01, eps = -0.01, p = FALSE))
  expect_equal(fractionate(d = 0.01, a = 0.99, p = FALSE), (0.01+1)/0.99 - 1)
  expect_equal(fractionate(d = 5, eps = 10, p = TRUE, mult = TRUE), ((5/1000+1) * (10/1000+1) - 1) * 1000)

  # test flux function
  expect_error(dir_flux(dir = "test"), "direction not recognized")
  expect_error(dir_flux(5, rev = 2, dir = "+"), "negative directional flux does not make sense")
  expect_equal(dir_flux(5, rev = 1, dir = "+", model_offset = 1e-12) + dir_flux(5, rev = 1, dir = "-", model_offset = 1e-12), 5)
  expect_equal(dir_flux(10, rev = 0, dir = "+"), 10)
  expect_equal(dir_flux(10, rev = 0, dir = "-"), 0)
  expect_equal(dir_flux(10, rev = 0.5, dir = "+"), 20)
  expect_equal(dir_flux(10, rev = 0.5, dir = "-"), -10)
})


test_that("Running model works", {

  expect_is({
    sys <- isopath() %>%
    add_isotope("C") %>% add_isotope("N") %>%
    add_component("X", C, N) %>% add_component("Y", C, N) %>%
    add_custom_reaction(X == 2.5 * Y, name = "my_rxn", flux.N = dN, flux.X.C = X.dC) %>%
    set_parameters(X.C = 1, X.N = 1, Y = 1, Y.C = 1, Y.N = 1)
    sys}, "isopath")

  expect_error(run_model("error"), "can only run model for an isopath")
  expect_error(run_model(sys), "time_steps is required")
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* object 'X' not found")
  expect_is({sys <- sys %>% set_parameters(transform(sys$parameters, X = c(1, 1))); sys}, "isopath")
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* there seem to be multiple identical run scenarios")
  expect_is({sys <- sys %>% set_parameters(transform(sys$parameters, X = c(1, 2))); sys}, "isopath")
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* object .* not found")
  expect_is({sys <- sys %>% set_parameters(dm = 0.1, dN = -5, X.dC = 10); sys}, "isopath")
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* some derivatives could not be computed")
  expect_is({sys <- sys %>%
    add_custom_reaction(X == 2.5 * Y, name = "my_rxn", flux = dm, flux.N = dN, flux.X.C = X.dC, flux.Y.C = Y.dC) %>%
    set_parameters(Y.dC = 5); sys}, "isopath")
  expect_message(sys %>% run_model(2), "Running model for 2 scenario(s)*")
  expect_message(sys %>% set_parameters(dm = 0.5) %>% run_model(2), "encountered the following error .* depleted .* pools")
  expect_equal(
    sys %>% set_parameters(dm = c(0.2, 0.4))  %>%  run_model(2) %>% signif(3),
    data_frame(
      X.dC = 10, Y.dC = 5, dN = -5,
      dm = c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4),
      time = c(0, 1, 2, 0, 1, 2),
      X = c(1, 0.8, 0.6, 2, 1.6, 1.2),
      X.C = c(1, -1.25, -5, 1, -1.25, -5),
      X.N = c(1, 2.5, 5, 1, 2.5, 5),
      Y = c(1, 1.5, 2, 1, 2, 3),
      Y.C = c(1, 2.33, 3, 1, 3, 3.67),
      Y.N = c(1, -1, -2, 1, -2, -3)))

})


test_that("Running steady-state works", {
  # @TODO
  expect_is({
    sys <- isopath() %>%
    add_isotope("C") %>% add_isotope("N") %>%
    add_component("X", C, N) %>% add_component("Y", C, N) %>%
    add_custom_reaction(X == 2.5 * Y, name = "my_rxn", flux = dm, flux.N = dN, flux.X.C = X.dC) %>%
    set_parameters(X.C = 1, X.N = 1, Y = 1, Y.C = 1, Y.N = 1, X.dC = 0)
    sys}, "isopath")

  expect_error(run_steady_state(NULL), "can only run model for an isopath")
  expect_error(run_steady_state(sys), "encountered the following error during pre-check .* object 'X' not found")
  expect_error( sys %>% set_parameters(transform(sys$parameters, X = c(1, 1))) %>% run_steady_state(),
                "encountered the following error during pre-check .* there seem to be multiple identical run scenarios")
  expect_error(run_steady_state(sys), "encountered the following error during pre-check .* object .* not found")
  expect_is({sys <- sys %>% set_parameters(dm = 0.1, dN = -5, X.dC = 10, X = 1); sys}, "isopath")
  expect_message(tryCatch(sys %>% set_parameters(dm = -1) %>% run_steady_state(), error = function(e){}), "depleted .* pool")
  expect_error(capture.output(sys %>% set_parameters(dm = -1) %>% run_steady_state()),
               "None of the scenarios could be run to steady-state")
  expect_message(tryCatch(sys %>% add_component("X", variable = F, C, N) %>%
                            set_parameters(dm = 1e6) %>% run_steady_state(), error = function(e){}), "overflowing pool")

})
