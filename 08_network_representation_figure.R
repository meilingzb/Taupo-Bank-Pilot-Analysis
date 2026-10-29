# 08 Figure 4: alternative network representations
#
# Purpose:
# Compare three representations of the Taupo Bank reference-scenario network:
# (a) the full multilayer network, (b) depth strata treated independently, and
# (c) a depth-integrated network with one node per seamount.
#
# Mathematical graph representation of the Taupo Bank pilot system

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(patchwork)
library(grid)

# Inputs

analysis_file <- file.path(
  "Outputs",
  "connectivity_analysis.rds"
)

if (!file.exists(analysis_file)) {
  stop(
    "Could not find Outputs/connectivity_analysis.rds. ",
    "Run 04_network_analysis.R first."
  )
}

connectivity_analysis <- readRDS(
  analysis_file
)

reference_pld_days <- 60
reference_precompetency_days <- 21

depth_levels <- c(
  "200-700 m",
  "700-1500 m",
  "1500-2500 m",
  ">2500 m"
)

depth_display_labels <- c(
  "200-700 m" = "200–700 m",
  "700-1500 m" = "700–1500 m",
  "1500-2500 m" = "1500–2500 m",
  ">2500 m" = "≥2500 m"
)

site_order <- c(
  "T22",
  "T21",
  "T20",
  "L11"
)

# These fixed positions reproduce the proposal figure geometry.
site_x <- c(
  1.00,
  2.55,
  4.10,
  5.65
)

output_dir <- file.path(
  "Outputs",
  "Proposal_figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

edge_scale_max <- 0.10

# Reference-scenario connectivity

depth_resolved <- connectivity_analysis$connectivity %>%
  filter(
    analysis_period == "pooled",
    period_type == "pooled",
    pld_days == reference_pld_days,
    precompetency_days == reference_precompetency_days,
    source_node == "T22",
    connectivity > 0
  ) %>%
  mutate(
    source_node = as.character(source_node),
    destination_node = as.character(destination_node),
    source_depth_stratum = as.character(source_depth_stratum),
    destination_depth_stratum = as.character(destination_depth_stratum)
  )

# The depth-integrated table is already source-normalised in 04 by averaging
# equally across source depth strata represented at each source seamount.
depth_integrated <- connectivity_analysis$collapsed_connectivity %>%
  filter(
    analysis_period == "pooled",
    period_type == "pooled",
    pld_days == reference_pld_days,
    precompetency_days == reference_precompetency_days,
    source_node == "T22",
    connectivity > 0
  ) %>%
  mutate(
    source_node = as.character(source_node),
    destination_node = as.character(destination_node)
  )

# Nodes and edges: depth-resolved representations

depth_node_table <- tidyr::expand_grid(
  site = site_order,
  depth = depth_levels
) %>%
  mutate(
    name = paste(
      site,
      depth,
      sep = "|"
    ),
    x = site_x[
      match(
        site,
        site_order
      )
    ],
    y = 5 -
      match(
        depth,
        depth_levels
      )
  )

depth_edge_table <- depth_resolved %>%
  transmute(
    from = paste(
      source_node,
      source_depth_stratum,
      sep = "|"
    ),
    to = paste(
      destination_node,
      destination_depth_stratum,
      sep = "|"
    ),
    connectivity = connectivity,
    connection_type = if_else(
      source_depth_stratum ==
        destination_depth_stratum,
      "Within-depth",
      "Cross-depth"
    ),
    same_site_cross_depth =
      source_node ==
        destination_node &
        source_depth_stratum !=
        destination_depth_stratum
  )

graph_multilayer <- tbl_graph(
  nodes = depth_node_table,
  edges = depth_edge_table,
  node_key = "name",
  directed = TRUE
)

graph_independent <- tbl_graph(
  nodes = depth_node_table,
  edges = depth_edge_table %>%
    filter(
      connection_type ==
        "Within-depth"
    ),
  node_key = "name",
  directed = TRUE
)

# Manual same-site cross-depth arrows

# These are plotted directly between the T22 site-depth nodes so the
# dashed cross-depth arrows do not collide with the open node circles.

depth_y_lookup <- c(
  "200-700 m" = 4,
  "700-1500 m" = 3,
  "1500-2500 m" = 2,
  ">2500 m" = 1
)

same_site_cross_segments <- depth_resolved %>%
  filter(
    source_node ==
      destination_node,
    source_depth_stratum !=
      destination_depth_stratum
  ) %>%
  transmute(
    source_depth_stratum,
    destination_depth_stratum,
    connectivity,
    y_source = unname(
      depth_y_lookup[
        source_depth_stratum
      ]
    ),
    y_destination = unname(
      depth_y_lookup[
        destination_depth_stratum
      ]
    )
  ) %>%
  mutate(
    x = site_x[1],
    xend = site_x[1],
    direction =
      sign(
        y_destination -
          y_source
      ),
    node_gap = 0.12,
    y =
      y_source +
      node_gap * direction,
    y_terminal =
      y_destination -
      node_gap * direction,
    arrowhead_length = 0.12,
    y_shaft_end =
      y_terminal -
      arrowhead_length * direction,
    plot_linewidth =
      0.42 +
      pmin(
        connectivity /
          edge_scale_max,
        1
      ) *
      (
        2.05 -
          0.42
      )
  )

# Nodes and edges: depth-integrated representation

integrated_node_table <- tibble(
  name = site_order,
  x = site_x,
  y = 1
)

integrated_edge_table <- depth_integrated %>%
  transmute(
    from = source_node,
    to = destination_node,
    connectivity = connectivity
  )

graph_integrated <- tbl_graph(
  nodes = integrated_node_table,
  edges = integrated_edge_table,
  node_key = "name",
  directed = TRUE
)

# Background layers

depth_bands <- tibble(
  ymin = c(
    3.54,
    2.54,
    1.54,
    0.54
  ),
  ymax = c(
    4.46,
    3.46,
    2.46,
    1.46
  ),
  band_fill = c(
    "grey98",
    "grey94",
    "grey98",
    "grey94"
  )
)

depth_guides <- tibble(
  y = c(
    4,
    3,
    2,
    1
  )
)

network_x_min <- 0.50
network_x_max <- 6.05
panel_ab_x_max <- 8.05

# Theme

theme_network <- function(
    base_size = 10.5) {

  theme_void(
    base_size = base_size
  ) +
    theme(
      text =
        element_text(
          colour = "black"
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 10.8,
          hjust = 0,
          margin = margin(
            b = 6
          )
        ),
      axis.text.x =
        element_text(
          colour = "black",
          size = 9,
          margin = margin(
            t = 5
          )
        ),
      axis.text.y =
        element_text(
          colour = "black",
          size = 8.8,
          margin = margin(
            r = 8
          )
        ),
      axis.title.x =
        element_text(
          colour = "black",
          face = "bold",
          size = 9.3,
          margin = margin(
            t = 7
          )
        ),
      plot.margin =
        margin(
          6,
          8,
          6,
          8
        )
    )
}

# Base depth-resolved panel

depth_panel_base <- function(
    graph_object,
    panel_title) {

  ggraph(
    graph_object,
    layout = "manual",
    x = x,
    y = y
  ) +

    geom_rect(
      data = depth_bands,
      aes(
        xmin = network_x_min,
        xmax = network_x_max,
        ymin = ymin,
        ymax = ymax,
        fill = band_fill
      ),
      inherit.aes = FALSE,
      colour = NA,
      show.legend = FALSE
    ) +

    scale_fill_identity(
      guide = "none"
    ) +

    geom_segment(
      data = depth_guides,
      aes(
        x = 0.68,
        xend = 5.92,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      colour = "grey79",
      linewidth = 0.34
    ) +

    scale_x_continuous(
      breaks = site_x,
      labels = site_order,
      limits = c(
        network_x_min,
        panel_ab_x_max
      ),
      expand = c(
        0,
        0
      )
    ) +

    scale_y_continuous(
      breaks = c(
        1,
        2,
        3,
        4
      ),
      labels =
        rev(
          unname(
            depth_display_labels[depth_levels]
          )
        ),
      limits = c(
        0.52,
        4.48
      ),
      expand = c(
        0,
        0
      )
    ) +

    labs(
      title = panel_title,
      x = NULL,
      y = NULL
    ) +

    coord_fixed(
      ratio = 1,
      clip = "off"
    ) +

    theme_network()
}

# Edge scales

edge_width_scale <- scale_edge_width_continuous(
  range = c(
    0.42,
    2.05
  ),
  limits = c(
    0,
    edge_scale_max
  ),
  guide = "none"
)

edge_type_scale <- scale_edge_linetype_manual(
  values = c(
    "Cross-depth" = "22",
    "Within-depth" = "solid"
  ),
  guide = "none"
)

# Panel a: multilayer network

panel_a <- depth_panel_base(
  graph_multilayer,
  "a) Multilayer network"
) +

  geom_edge_loop(
    aes(
      edge_width = connectivity,
      span = 300,
      direction = 180,
      strength = 1.25
    ),
    n = 160,
    colour = "black",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.25,
        "mm"
      ),
      type = "closed"
    ),
    end_cap = ggraph::circle(
      2.3,
      "mm"
    ),
    show.legend = FALSE
  ) +

  geom_edge_arc(
    aes(
      filter =
        from !=
        to &
        !same_site_cross_depth,
      edge_width = connectivity,
      edge_linetype =
        connection_type
    ),
    strength = 0.10,
    fold = TRUE,
    colour = "black",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.25,
        "mm"
      ),
      type = "closed"
    ),
    start_cap = ggraph::circle(
      2.1,
      "mm"
    ),
    end_cap = ggraph::circle(
      2.3,
      "mm"
    ),
    show.legend = FALSE
  ) +

  geom_segment(
    data =
      same_site_cross_segments,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = y_shaft_end,
      linewidth = plot_linewidth
    ),
    inherit.aes = FALSE,
    colour = "black",
    linetype = "longdash",
    lineend = "butt",
    show.legend = FALSE
  ) +

  geom_segment(
    data =
      same_site_cross_segments,
    aes(
      x = x,
      xend = xend,
      y = y_shaft_end,
      yend = y_terminal,
      linewidth = plot_linewidth
    ),
    inherit.aes = FALSE,
    colour = "black",
    linetype = "solid",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.45,
        "mm"
      ),
      angle = 24,
      type = "closed"
    ),
    show.legend = FALSE
  ) +

  geom_node_point(
    shape = 21,
    size = 3.2,
    stroke = 0.85,
    fill = "white",
    colour = "black"
  ) +

  scale_linewidth_identity(
    guide = "none"
  ) +

  edge_width_scale +
  edge_type_scale +

  annotate(
    "text",
    x = 6.43,
    y = 4.13,
    label =
      "Connectivity probability",
    hjust = 0,
    fontface = "bold",
    size = 3.25
  ) +

  geom_segment(
    aes(
      x = 6.45,
      xend = 7.00,
      y = 3.62,
      yend = 3.62
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth =
      0.42 +
      0.50 *
      (
        2.05 -
          0.42
      ),
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.15,
        "mm"
      ),
      type = "closed"
    )
  ) +

  annotate(
    "text",
    x = 7.15,
    y = 3.62,
    label = "5.0%",
    hjust = 0,
    size = 2.95
  ) +

  geom_segment(
    aes(
      x = 6.45,
      xend = 7.00,
      y = 3.05,
      yend = 3.05
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 2.05,
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.15,
        "mm"
      ),
      type = "closed"
    )
  ) +

  annotate(
    "text",
    x = 7.15,
    y = 3.05,
    label = "10.0%",
    hjust = 0,
    size = 2.95
  )

# Panel b: depth strata analysed independently

panel_b <- depth_panel_base(
  graph_independent,
  "b) Depth strata analysed independently"
) +

  geom_edge_loop(
    aes(
      edge_width = connectivity,
      span = 300,
      direction = 180,
      strength = 1.25
    ),
    n = 160,
    colour = "black",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.25,
        "mm"
      ),
      type = "closed"
    ),
    end_cap = ggraph::circle(
      2.3,
      "mm"
    ),
    show.legend = FALSE
  ) +

  geom_edge_arc(
    aes(
      filter =
        from !=
        to,
      edge_width = connectivity
    ),
    strength = 0.10,
    fold = TRUE,
    colour = "black",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.25,
        "mm"
      ),
      type = "closed"
    ),
    start_cap = ggraph::circle(
      2.1,
      "mm"
    ),
    end_cap = ggraph::circle(
      2.3,
      "mm"
    ),
    show.legend = FALSE
  ) +

  geom_node_point(
    shape = 21,
    size = 3.2,
    stroke = 0.85,
    fill = "white",
    colour = "black"
  ) +

  scale_edge_width_continuous(
    range = c(
      0.42,
      2.05
    ),
    limits = c(
      0,
      edge_scale_max
    ),
    guide = "none"
  ) +

  annotate(
    "text",
    x = 6.43,
    y = 4.13,
    label = "Connection type",
    hjust = 0,
    fontface = "bold",
    size = 3.25
  ) +

  geom_segment(
    aes(
      x = 6.45,
      xend = 7.00,
      y = 3.55,
      yend = 3.55
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.75,
    linetype = "longdash",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.75,
        "mm"
      ),
      angle = 24,
      type = "open"
    )
  ) +

  annotate(
    "text",
    x = 7.15,
    y = 3.55,
    label = "Cross-depth",
    hjust = 0,
    size = 2.95
  ) +

  geom_segment(
    aes(
      x = 6.45,
      xend = 7.00,
      y = 2.95,
      yend = 2.95
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.75,
    linetype = "solid",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.15,
        "mm"
      ),
      type = "closed"
    )
  ) +

  annotate(
    "text",
    x = 7.15,
    y = 2.95,
    label = "Within-depth",
    hjust = 0,
    size = 2.95
  )

# Panel c: depth-integrated network

panel_c <- ggraph(
  graph_integrated,
  layout = "manual",
  x = x,
  y = y
) +

  geom_rect(
    aes(
      xmin = network_x_min,
      xmax = network_x_max,
      ymin = 0.72,
      ymax = 1.28
    ),
    inherit.aes = FALSE,
    fill = "grey97",
    colour = NA
  ) +

  geom_segment(
    aes(
      x = 0.68,
      xend = 5.92,
      y = 1,
      yend = 1
    ),
    inherit.aes = FALSE,
    colour = "grey79",
    linewidth = 0.34
  ) +

  geom_edge_loop(
    aes(
      edge_width = connectivity,
      span = 300,
      direction = 180,
      strength = 1.25
    ),
    n = 160,
    colour = "black",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.25,
        "mm"
      ),
      type = "closed"
    ),
    end_cap = ggraph::circle(
      2.3,
      "mm"
    ),
    show.legend = FALSE
  ) +

  geom_edge_arc(
    aes(
      filter =
        from !=
        to,
      edge_width = connectivity
    ),
    strength = 0.10,
    fold = TRUE,
    colour = "black",
    lineend = "round",
    arrow = arrow(
      length = unit(
        1.25,
        "mm"
      ),
      type = "closed"
    ),
    start_cap = ggraph::circle(
      2.1,
      "mm"
    ),
    end_cap = ggraph::circle(
      2.3,
      "mm"
    ),
    show.legend = FALSE
  ) +

  geom_node_point(
    shape = 21,
    size = 3.2,
    stroke = 0.85,
    fill = "white",
    colour = "black"
  ) +

  scale_x_continuous(
    breaks = site_x,
    labels = site_order,
    limits = c(
      network_x_min,
      network_x_max
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(
    limits = c(
      0.70,
      1.30
    ),
    breaks = NULL,
    expand = c(
      0,
      0
    )
  ) +

  labs(
    title =
      "c) Depth-integrated network",
    x = "Seamount",
    y = NULL
  ) +

  coord_fixed(
    ratio = 1,
    clip = "off"
  ) +

  scale_edge_width_continuous(
    range = c(
      0.42,
      2.05
    ),
    limits = c(
      0,
      edge_scale_max
    ),
    guide = "none"
  ) +

  theme_network() +

  theme(
    axis.text.y =
      element_blank(),
    axis.title.y =
      element_blank()
  )

# Assemble

figure_4 <- panel_a /
  panel_b /
  panel_c +
  plot_layout(
    heights = c(
      1.00,
      1.00,
      0.42
    )
  )

# Save

figure_4_file <- file.path(
  output_dir,
  "Figure_4_taupo_network.pdf"
)

ggsave(
  filename = figure_4_file,
  plot = figure_4,
  device = cairo_pdf,
  width = 10.0,
  height = 7.6,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_4_taupo_network.png"
  ),
  plot = figure_4,
  width = 10.0,
  height = 7.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

