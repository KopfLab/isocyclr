

#' plot reaction diagram
#'
#' @export
generate_reaction_diagram <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only generate diagram for isopath")

  # determin x locations for components
  general_x <-
    ip %>% get_reaction_component_matrix() %>%
    mutate(x = rxn_nr - ifelse(comp_stoic > 0, 1, 0 ))

  # determine y locations for components
  components_xy <-
    general_x %>% select(isotope, component, x) %>%
    distinct() %>% group_by(isotope, x) %>%
    mutate(y = seq(-length(x)+1, length(x)-1, length.out = length(x))) %>%
    ungroup()

  # determin reaction lines
  rxns_xy <-
    general_x %>% left_join(components_xy %>% select(-x), by = c("component", "isotope")) %>%
    mutate(
      xstart = ifelse(comp_stoic > 0, x, x - 0.5),
      xend = ifelse(comp_stoic > 0, x + 0.5, x),
      ystart = ifelse(comp_stoic > 0, y, 0),
      yend = ifelse(comp_stoic > 0, 0, y)
    )

  # plot everything
  ggplot() +
    geom_segment(
      data = rxns_xy, size = 1, alpha = 0.5,
      #arrow = arrow(length = unit(0.03, "npc")),
      map = aes(x=xstart, y=ystart,
                xend=xend, yend=yend,
                color = reaction)) +
    geom_label(
      data = components_xy, hjust = 0.5,
      map = aes(x, y, label = component, fill = component)) +
    facet_grid(isotope~.) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    scale_x_continuous("", expand = c(0, 0.5)) +
    scale_y_continuous("", expand = c(0, 0.5))
}
