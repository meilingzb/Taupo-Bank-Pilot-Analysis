# 07 Figure 3: geographic connectivity in the reference scenario
#
# Purpose:
# Map depth-resolved potential connectivity from Taupo Bank for the reference
# biological scenario (PLD = 60 d; pre-competency = 21 d).
#
# Visual encoding:
# - one panel per Taupo source-depth stratum;
# - line colour indicates destination depth;
# - line width indicates source-normalised connectivity probability; and
# - loops around T22 represent same-seamount retention.

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(terra)
library(scales)
library(grid)
library(patchwork)

# Load project inputs

source(
  "00_setup_and_inputs.R"
)

# Output directory

figure_output_dir <- file.path(
  "Outputs",
  "Proposal_figures"
)

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Load connectivity analysis

analysis_file <- file.path(
  "Outputs",
  "connectivity_analysis.rds"
)

if (!file.exists(analysis_file)) {

  stop(
    paste0(
      "Could not find ",
      analysis_file,
      ". Current working directory: ",
      getwd()
    )
  )
}

connectivity_analysis <- readRDS(
  analysis_file
)

connectivity <- connectivity_analysis$connectivity

# Reference biological scenario

reference_pld_days <- 60

reference_precompetency_days <- 21

# Depth strata

depth_levels <- c(
  "200-700 m",
  "700-1500 m",
  "1500-2500 m",
  ">2500 m"
)

# Display labels

depth_labels <- c(
  "200–700 m",
  "700–1500 m",
  "1500–2500 m",
  "≥2500 m"
)

names(
  depth_labels
) <- depth_levels

# Destination-depth colours

depth_colours <- c(
  "#0072B2",
  "#009E73",
  "#E69F00",
  "#7A3E9D"
)

names(
  depth_colours
) <- depth_levels

# Helper: curvature for external connections

get_curve <- function(depth) {

  depth <- as.character(
    depth
  )

  case_when(
    depth == "200-700 m" ~ -0.20,
    depth == "700-1500 m" ~ -0.07,
    depth == "1500-2500 m" ~ 0.07,
    depth == ">2500 m" ~ 0.20,
    TRUE ~ 0
  )
}

# Helper: radius for same-site retention loops

get_loop_radius <- function(depth) {

  depth <- as.character(
    depth
  )

  case_when(
    depth == "200-700 m" ~ 0.08,
    depth == "700-1500 m" ~ 0.11,
    depth == "1500-2500 m" ~ 0.14,
    depth == ">2500 m" ~ 0.17,
    TRUE ~ 0.10
  )
}

# Select pooled reference-scenario connectivity

reference_connectivity <- connectivity %>%

  filter(
    analysis_period == "pooled",
    source_node == "T22",
    pld_days == reference_pld_days,
    precompetency_days == reference_precompetency_days,
    connectivity > 0
  )

# Some versions of the analysis output include period_type.
# Where present, retain pooled records only.

if (
  "period_type" %in%
  names(
    reference_connectivity
  )
) {

  reference_connectivity <- reference_connectivity %>%

    filter(
      period_type == "pooled"
    )
}

# Standardise depth ordering

reference_connectivity <- reference_connectivity %>%

  mutate(
    source_depth_stratum = factor(
      source_depth_stratum,
      levels = depth_levels
    ),
    destination_depth_stratum = factor(
      destination_depth_stratum,
      levels = depth_levels
    )
  )

# Validate that the reference biological scenario contains non-zero connections.

if (
  nrow(
    reference_connectivity
  ) == 0
) {

  stop(
    paste0(
      "No records were found for the ",
      reference_pld_days,
      "-day PLD / ",
      reference_precompetency_days,
      "-day pre-competency reference scenario."
    )
  )
}

# Prepare site coordinates

plot_nodes <- candidate_nodes %>%

  distinct(
    node_id,
    .keep_all = TRUE
  ) %>%

  mutate(
    plot_lon = lon,
    plot_lat = lat
  )

# Separate external connectivity from same-site retention

external_edges <- reference_connectivity %>%

  filter(
    destination_node != source_node
  ) %>%

  left_join(
    plot_nodes %>%
      transmute(
        source_node = node_id,
        source_lon = plot_lon,
        source_lat = plot_lat
      ),
    by = "source_node"
  ) %>%

  left_join(
    plot_nodes %>%
      transmute(
        destination_node = node_id,
        destination_lon = plot_lon,
        destination_lat = plot_lat
      ),
    by = "destination_node"
  )

retention_edges <- reference_connectivity %>%

  filter(
    destination_node == source_node
  )

# Taupo Bank coordinates

taupo_point <- plot_nodes %>%

  filter(
    node_id == "T22"
  )

if (
  nrow(
    taupo_point
  ) == 0
) {

  stop(
    "Taupo Bank T22 was not found in candidate_nodes."
  )
}

taupo_lon <- taupo_point %>%

  pull(
    plot_lon
  ) %>%

  first()

taupo_lat <- taupo_point %>%

  pull(
    plot_lat
  ) %>%

  first()

# Create same-site retention loops
#
# The latitude radius is increased slightly so loops appear approximately
# circular under the final non-fixed plotting aspect.
#
# This affects only the plotted loop shape; it does not alter connectivity.

make_retention_loops <- function(edges) {

  if (
    nrow(
      edges
    ) == 0
  ) {

    return(
      tibble()
    )
  }

  loop_list <- lapply(
    seq_len(
      nrow(
        edges
      )
    ),
    function(i) {
      edge <- edges %>%
        slice(
          i
        )
      source_depth <- edge %>%
        pull(
          source_depth_stratum
        ) %>%
        as.character()
      destination_depth <- edge %>%
        pull(
          destination_depth_stratum
        ) %>%
        as.character()
      connectivity_value <- edge %>%
        pull(
          connectivity
        )
      longitude_radius <- get_loop_radius(
        destination_depth
      )
      latitude_radius <- 1.50 *
        longitude_radius
      # Leave a small gap so the arrowhead shows direction.
      theta <- seq(
        pi / 4,
        2 * pi + pi / 4 - 0.25,
        length.out = 90
      )
      tibble(
        loop_id = paste(
          source_depth,
          destination_depth,
          sep = "_"
        ),
        source_depth_stratum = factor(
          source_depth,
          levels = depth_levels
        ),
        destination_depth_stratum = factor(
          destination_depth,
          levels = depth_levels
        ),
        connectivity = connectivity_value,
        x = taupo_lon +
          longitude_radius *
          cos(theta),
        y = taupo_lat +
          latitude_radius *
          sin(theta)
      )
    }
  )

  bind_rows(
    loop_list
  )
}

retention_loops <- make_retention_loops(
  retention_edges
)

# Determine geographic extent

connected_ids <- unique(
  c(
    "T22",
    external_edges %>%
      pull(
        destination_node
      ) %>%
      as.character()
  )
)

extent_nodes <- plot_nodes %>%

  filter(
    node_id %in%
      connected_ids
  )

longitude_range <- extent_nodes %>%

  pull(
    plot_lon
  ) %>%

  range(
    na.rm = TRUE
  )

latitude_range <- extent_nodes %>%

  pull(
    plot_lat
  ) %>%

  range(
    na.rm = TRUE
  )

# Map margins are deliberately compact so the connected sites occupy most of
# each panel when the figure is inserted at proposal-page width.

longitude_min <- min(
  longitude_range
) - 0.40

longitude_max <- max(
  longitude_range
) + 0.40

latitude_min <- min(
  latitude_range
) - 0.45

latitude_max <- max(
  latitude_range
) + 0.45

# Nearby candidate sites

background_nodes <- plot_nodes %>%

  filter(
    plot_lon >= longitude_min,
    plot_lon <= longitude_max,
    plot_lat >= latitude_min,
    plot_lat <= latitude_max
  )

# Connected nodes

connected_nodes <- plot_nodes %>%

  filter(
    node_id %in%
      connected_ids
  )

# GEBCO bathymetry

gebco_crop <- terra::crop(
  gebco,
  terra::ext(
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max
  )
)

# Slight aggregation improves plotting speed while retaining
# broad seabed morphology.

gebco_plot <- terra::aggregate(
  gebco_crop,
  fact = 3,
  fun = mean,
  na.rm = TRUE
)

gebco_df <- as.data.frame(
  gebco_plot,
  xy = TRUE,
  na.rm = FALSE
)

names(
  gebco_df
) <- c(
  "x",
  "y",
  "elevation"
)

# Bathymetric contours
#
# Draw 250 m contours throughout the complete mapped depth range.

minimum_bathymetry <- floor(
  min(
    gebco_df$elevation,
    na.rm = TRUE
  ) / 250
) * 250

background_contours <- seq(
  from = minimum_bathymetry,
  to = -250,
  by = 250
)

# Begin geographic figure

figure_3_map <- ggplot() +

  # GEBCO bathymetry

  geom_contour(
    data = gebco_df,
    aes(
      x = x,
      y = y,
      z = elevation
    ),
    breaks = background_contours,
    colour = "grey87",
    linewidth = 0.20
  ) +

  # Nearby candidate nodes

  geom_point(
    data = background_nodes,
    aes(
      x = plot_lon,
      y = plot_lat
    ),
    shape = 21,
    size = 1.4,
    stroke = 0.35,
    fill = "white",
    colour = "grey72"
  )

# External connectivity

for (
  depth_current in depth_levels
) {

  edge_subset <- external_edges %>%

    filter(
      as.character(
        destination_depth_stratum
      ) == depth_current
    )

  if (
    nrow(
      edge_subset
    ) > 0
  ) {

    figure_3_map <- figure_3_map +

      geom_curve(
        data = edge_subset,
        aes(
          x = source_lon,
          y = source_lat,
          xend = destination_lon,
          yend = destination_lat,
          linewidth = connectivity,
          colour = destination_depth_stratum
        ),
        curvature = get_curve(
          depth_current
        ),
        alpha = 0.92,
        lineend = "round",
        arrow = grid::arrow(
          length = grid::unit(
            0.12,
            "cm"
          ),
          type = "closed"
        ),
        show.legend = FALSE
      )
  }
}

# Same-site retention

if (
  nrow(
    retention_loops
  ) > 0
) {

  figure_3_map <- figure_3_map +

    geom_path(
      data = retention_loops,
      aes(
        x = x,
        y = y,
        group = loop_id,
        linewidth = connectivity,
        colour = destination_depth_stratum
      ),
      alpha = 0.95,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.10,
          "cm"
        ),
        type = "closed"
      ),
      show.legend = FALSE
    )
}

# Nodes, labels and facets

figure_3_map <- figure_3_map +

  # External destinations

  geom_point(
    data = connected_nodes %>%
      filter(
        node_id != "T22"
      ),
    aes(
      x = plot_lon,
      y = plot_lat
    ),
    shape = 21,
    size = 2.3,
    stroke = 0.6,
    fill = "white",
    colour = "black"
  ) +

  # Taupo source

  geom_point(
    data = taupo_point,
    aes(
      x = plot_lon,
      y = plot_lat
    ),
    shape = 21,
    size = 3.4,
    stroke = 0.7,
    fill = "black",
    colour = "black"
  ) +

  # External labels

  geom_text_repel(
    data = connected_nodes %>%
      filter(
        node_id != "T22"
      ),
    aes(
      x = plot_lon,
      y = plot_lat,
      label = node_id
    ),
    size = 3,
    box.padding = 0.45,
    point.padding = 0.35,
    force = 1.8,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.colour = "grey45",
    seed = 1,
    max.overlaps = Inf
  ) +

  # Taupo label

  geom_text(
    data = taupo_point,
    aes(
      x = plot_lon,
      y = plot_lat,
      label = node_id
    ),
    nudge_x = 0.15,
    nudge_y = -0.15,
    size = 3
  ) +

  # Source-depth panels

  facet_wrap(
    ~ source_depth_stratum,
    nrow = 1,
    drop = FALSE,
    labeller = labeller(
      source_depth_stratum =
        depth_labels
    )
  ) +

  # Destination-depth colours

  scale_colour_manual(
    values = depth_colours,
    breaks = depth_levels,
    labels = depth_labels,
    limits = depth_levels,
    drop = FALSE,
    guide = "none"
  ) +

  # Connectivity line widths

  scale_linewidth_continuous(
    range = c(
      0.35,
      2.8
    ),
    limits = c(
      0,
      0.10
    ),
    breaks = c(
      0.01,
      0.05,
      0.10
    ),
    guide = "none"
  ) +

  # Normal Cartesian coordinates
  #
  # This lets the final figure dimensions control the panel shape.

  coord_cartesian(
    xlim = c(
      longitude_min,
      longitude_max
    ),
    ylim = c(
      latitude_min,
      latitude_max
    ),
    expand = FALSE,
    clip = "off"
  ) +

  labs(
    x = "Longitude",
    y = "Latitude"
  ) +

  # Style matched to Figure 2

  theme_bw(
    base_size = 10.5
  ) +

  theme(
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      size = 13,
      margin = margin(
        t = 1,
        b = 6
      )
    ),
    axis.title = element_text(
      face = "plain",
      colour = "black",
      size = 11
    ),
    axis.text = element_text(
      colour = "black",
      size = 9
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.40
    ),
    axis.ticks.length = unit(
      0.11,
      "cm"
    ),
    panel.spacing.x = grid::unit(
      0.35,
      "lines"
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 7,
      r = 4,
      b = 0,
      l = 4
    )
  )

# Compact two-row legend
#
# Both rows are contained within one plotting panel so they appear
# as one wrapped legend rather than two separate legend blocks.
#
# The vertical positions are deliberately close:
#   row 1 = 1.10
#   row 2 = 0.90

legend_panel <- ggplot() +

  # ---------------------------------------------------------------
# Row 1: connectivity
# ---------------------------------------------------------------

annotate(
  "text",
  x = 0.15,
  y = 1.10,
  label = "Connectivity",
  hjust = 0,
  fontface = "bold",
  size = 2.9
) +

  # 1%

  geom_segment(
    aes(
      x = 2.00,
      xend = 2.38,
      y = 1.10,
      yend = 1.10
    ),
    colour = "black",
    linewidth = 0.45,
    lineend = "round",
    arrow = arrow(
      length = unit(
        0.95,
        "mm"
      ),
      type = "closed"
    )
  ) +

  annotate(
    "text",
    x = 2.52,
    y = 1.10,
    label = "1.0%",
    hjust = 0,
    size = 2.6
  ) +

  # 5%

  geom_segment(
    aes(
      x = 3.45,
      xend = 3.87,
      y = 1.10,
      yend = 1.10
    ),
    colour = "black",
    linewidth = 1.25,
    lineend = "round",
    arrow = arrow(
      length = unit(
        0.95,
        "mm"
      ),
      type = "closed"
    )
  ) +

  annotate(
    "text",
    x = 4.01,
    y = 1.10,
    label = "5.0%",
    hjust = 0,
    size = 2.6
  ) +

  # 10%

  geom_segment(
    aes(
      x = 5.00,
      xend = 5.47,
      y = 1.10,
      yend = 1.10
    ),
    colour = "black",
    linewidth = 2.10,
    lineend = "round",
    arrow = arrow(
      length = unit(
        0.95,
        "mm"
      ),
      type = "closed"
    )
  ) +

  annotate(
    "text",
    x = 5.62,
    y = 1.10,
    label = "10.0%",
    hjust = 0,
    size = 2.6
  ) +

  # ---------------------------------------------------------------
# Row 2: destination depth and retention
# ---------------------------------------------------------------

annotate(
  "text",
  x = 0.15,
  y = 0.90,
  label = "Destination depth",
  hjust = 0,
  fontface = "bold",
  size = 2.9
) +

  # 200–700 m

  geom_segment(
    aes(
      x = 2.55,
      xend = 2.98,
      y = 0.90,
      yend = 0.90
    ),
    colour = "#0072B2",
    linewidth = 2.25,
    lineend = "round"
  ) +

  annotate(
    "text",
    x = 3.10,
    y = 0.90,
    label = "200–700 m",
    hjust = 0,
    size = 2.55
  ) +

  # 700–1500 m

  geom_segment(
    aes(
      x = 4.65,
      xend = 5.08,
      y = 0.90,
      yend = 0.90
    ),
    colour = "#009E73",
    linewidth = 2.25,
    lineend = "round"
  ) +

  annotate(
    "text",
    x = 5.20,
    y = 0.90,
    label = "700–1500 m",
    hjust = 0,
    size = 2.55
  ) +

  # 1500–2500 m

  geom_segment(
    aes(
      x = 6.92,
      xend = 7.35,
      y = 0.90,
      yend = 0.90
    ),
    colour = "#E69F00",
    linewidth = 2.25,
    lineend = "round"
  ) +

  annotate(
    "text",
    x = 7.47,
    y = 0.90,
    label = "1500–2500 m",
    hjust = 0,
    size = 2.55
  ) +

  # ≥2500 m

  geom_segment(
    aes(
      x = 9.42,
      xend = 9.85,
      y = 0.90,
      yend = 0.90
    ),
    colour = "#7A3E9D",
    linewidth = 2.25,
    lineend = "round"
  ) +

  annotate(
    "text",
    x = 9.97,
    y = 0.90,
    label = "≥2500 m",
    hjust = 0,
    size = 2.55
  ) +

  # Same-site retention

  geom_point(
    aes(
      x = 11.55,
      y = 0.90
    ),
    shape = 21,
    size = 4.6,
    stroke = 0.9,
    fill = "white",
    colour = "#009E73"
  ) +

  geom_point(
    aes(
      x = 11.55,
      y = 0.90
    ),
    shape = 16,
    size = 1.55,
    colour = "black"
  ) +

  annotate(
    "text",
    x = 11.90,
    y = 0.90,
    label = "Same-site retention",
    hjust = 0,
    size = 2.55
  ) +

  # Compact legend viewport

  coord_cartesian(
    xlim = c(
      0,
      14.4
    ),
    ylim = c(
      0.78,
      1.22
    ),
    clip = "off"
  ) +

  theme_void() +

  theme(
    plot.margin = margin(
      0,
      0,
      0,
      0
    )
  )

# Combine map and compact wrapped legend
#
# The legend receives enough height for two close rows, but remains
# visually subordinate to the geographic panels.

figure_3 <- figure_3_map /
  legend_panel +
  plot_layout(
    heights = c(
      12,
      1.45
    )
  )

# Save the figure as a vector PDF
#
# A less panoramic shape means the figure can occupy more vertical
# space when inserted at page width in the proposal.

figure_3_file <- file.path(
  figure_output_dir,
  "Figure_3_reference_geographic_connectivity_60dPLD_21dPC.pdf"
)

ggsave(
  filename = figure_3_file,
  plot = figure_3,
  device = cairo_pdf,
  width = 9.0,
  height = 6.2,
  units = "in",
  bg = "white"
)

