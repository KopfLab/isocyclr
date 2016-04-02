context("Steady-state")

test_that("Testing simple reactions steady state", {
  sys <- isopath() %>%
    add_isotope(c("carbon", "nitrogen")) %>%
    add_component("A", carbon, nitrogen) %>%
    add_component("B", carbon, nitrogen) %>%
    add_component("C", carbon)
  #expect_error( sys %>% add_simple_reaction(A == B + C), "simple reactions can only be 1 to 1 transformations .* carbon")
  #expect_error( sys %>% add_simple_reaction(A == B), "simple reactions can only be 1 to 1 transformations .* nitrogen")
})
