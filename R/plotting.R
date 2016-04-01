

#' plot reaction diagram
#'
#' @export
generate_reaction_diagram <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only generate diagram for isopath", call. = FALSE)

  if (length(ip$reactions) == 0 || length(ip$components) == 0)
    stop("need at least two components and one reaction in order to plot the reaction diagram", call. = FALSE)

  # determin x locations for components
  general_x <-
    ip %>% get_reaction_component_matrix() %>%
    mutate(x = abscissa)

  # determine y locations for components
  components_xy <-
    general_x %>%
    select(isotope, component, x) %>%
    distinct() %>% arrange(component) %>%
    group_by(isotope, x) %>%
    mutate(y = seq(-length(x)+1, length(x)-1, length.out = length(x))) %>%
    ungroup()

  # reaction lines (offset for parallel reaction lines)
  rxn_components_xy <-
    left_join(
      general_x %>% select(reaction, component, isotope, comp_stoic, x),
      components_xy,
      by = c("component", "isotope", "x"))
  rxns_xy <-
    inner_join(
      rxn_components_xy %>% filter(comp_stoic > 0) %>% rename(xstart = x, ystart = y),
      rxn_components_xy %>% filter(comp_stoic < 0) %>% rename(xend = x, yend = y),
      by = c("reaction", "isotope")
    ) %>%
    group_by(xstart, ystart, xend, yend) %>% arrange(reaction) %>%
    # offset in y direction depending on multi-line reactions
    mutate(y_offset = 0.05 * seq(-n()+1, n()-1, length.out = n())) %>%
    ungroup() %>%
    mutate(
      ystart = ystart + y_offset,
      yend = yend + y_offset
    )

  # plot everything
  ggplot() +
    geom_segment(
      data = rxns_xy, size = 1.5, alpha = 0.8,
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
      axis.ticks = element_blank(),
      legend.box = "horizontal"
    ) +
    scale_x_continuous("", expand = c(0, 0.5)) +
    scale_y_continuous("", expand = c(0, 0.5))
}
