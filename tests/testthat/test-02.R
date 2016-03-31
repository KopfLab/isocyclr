context("Modelling")

test_that("Running model works", {

  sys <- isopath() %>%
    add_isotope("C") %>% add_isotope("N") %>%
    add_component("X", C, N) %>% add_component("Y", C, N) %>%
    add_custom_reaction("rxn1", X == 2.5 * Y, flux = dm, flux.N = dN, flux.X.C = X.dC) %>%
    set_parameters(X.C = 1, X.N = 1, Y = 1, Y.C = 1, Y.N = 1)

  expect_error(run_model(NULL), "can only run model for an isopath")
  expect_error(run_model(sys), "time_steps is required")
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* parameters required for minimal parameter set missing")
  sys <- sys %>% set_parameters(X = c(1, 1))
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* there seem to be multiple identical run scenarios")
  sys <- sys %>% set_parameters(X = c(1, 2))
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* object .* not found")
  sys <- sys %>% set_parameters(dm = 0.1, dN = -5, X.dC = 10)
  expect_error(run_model(sys, 5), "encountered the following error during pre-check .* missing isotope flux")
  sys <- sys %>%
    add_custom_reaction("rxn1", X == 2.5 * Y, flux = dm, flux.N = dN, flux.X.C = X.dC, flux.Y.C = Y.dC) %>%
    set_parameters(Y.dC = 5)
  expect_message(sys %>% run_model(2), "Running model for 2 scenarios")
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