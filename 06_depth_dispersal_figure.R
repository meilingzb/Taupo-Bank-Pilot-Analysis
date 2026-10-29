# 06 Figure 2: depth-dependent dispersal from Taupo Bank
#
# Purpose:
# Show how the maximum horizontal displacement reached by each pilot particle
# varies with release depth and release month. The figure is a transport
# diagnostic: it describes physical dispersal from Taupo Bank rather than
# habitat encounters or connectivity.
#
# Calculation:
# For each particle, great-circle distance is calculated from the release
# position to every trajectory position. The maximum of those distances is then
# retained. Black points and lines show the mean maximum displacement at each
# release depth, while shaded bands identify the four analysis depth strata.

library(dplyr)
library(ggplot2)
library(tibble)
library(grid)

# File paths

particle_file <- file.path(
  "Outputs",
  "particle_trajectories.rds"
)

summary_file <- file.path(
  "Outputs",
  "Taupo_pilot",
  "taupo_particle_maximum_horizontal_displacement.csv"
)

figure_output_dir <- file.path(
  "Outputs",
  "Proposal_figures"
)

# Create output directory if required

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Release-month lookup

# The three release IDs correspond to the three monthly pilot releases.

release_month_lookup <- tibble(
  release_id = c(
    1,
    2,
    3
  ),
  release_month = factor(
    c(
      "June",
      "July",
      "August"
    ),
    levels = c(
      "June",
      "July",
      "August"
    )
  )
)

# Haversine distance

# Calculates great-circle horizontal distance between two geographic
# coordinates in kilometres.

haversine_km <- function(
    lon1,
    lat1,
    lon2,
    lat2) {

  earth_radius_km <- 6371.0088

  to_radians <- pi / 180

  lon1_rad <- lon1 * to_radians
  lat1_rad <- lat1 * to_radians

  lon2_rad <- lon2 * to_radians
  lat2_rad <- lat2 * to_radians

  delta_lon <- lon2_rad - lon1_rad
  delta_lat <- lat2_rad - lat1_rad

  a <- sin(delta_lat / 2)^2 +
    cos(lat1_rad) *
    cos(lat2_rad) *
    sin(delta_lon / 2)^2

  # Numerical protection against very small floating-point errors.

  a <- pmin(
    1,
    pmax(
      0,
      a
    )
  )

  2 *
    earth_radius_km *
    asin(
      sqrt(
        a
      )
    )
}

# Calculate particle maximum horizontal displacement

if (file.exists(particle_file)) {

  # Preferred workflow:
  # calculate maximum horizontal displacement directly from the
  # complete particle trajectories.

  particle_trajectories <- readRDS(
    particle_file
  )

  # Identify each particle's initial release position.
  #
  # The trajectory record with the lowest particle age is taken as
  # the release position.

  release_positions <- particle_trajectories %>%

    group_by(
      particle_id
    ) %>%

    slice_min(
      order_by = age_hours,
      n = 1,
      with_ties = FALSE
    ) %>%

    ungroup() %>%

    transmute(
      particle_id,
      release_id,
      source_depth_m,
      release_lon = lon,
      release_lat = lat
    )

  # Calculate horizontal displacement from the release position at
  # every trajectory time step, then retain the largest displacement
  # reached by each particle.

  particle_displacement <- particle_trajectories %>%

    select(
      particle_id,
      lon,
      lat
    ) %>%

    left_join(
      release_positions,
      by = "particle_id"
    ) %>%

    mutate(
      horizontal_displacement_km = haversine_km(
        release_lon,
        release_lat,
        lon,
        lat
      )
    ) %>%

    group_by(
      particle_id,
      release_id,
      source_depth_m
    ) %>%

    summarise(
      maximum_horizontal_displacement_km = max(
        horizontal_displacement_km,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

} else if (file.exists(summary_file)) {

  # Fallback workflow:
  # use the previously calculated particle-level summary.

  particle_displacement <- read.csv(
    summary_file
  ) %>%

    as_tibble() %>%

    select(
      particle_id,
      release_id,
      source_depth_m,
      maximum_horizontal_displacement_km
    )

} else {

  stop(
    paste0(
      "Could not find either expected input file.\n\n",
      "Checked:\n",
      particle_file,
      "\n",
      summary_file,
      "\n\nCurrent working directory:\n",
      getwd()
    )
  )
}

# Add release month

particle_displacement <- particle_displacement %>%

  left_join(
    release_month_lookup,
    by = "release_id"
  ) %>%

  # The project focuses on deep-sea habitat at or below 200 m.

  filter(
    source_depth_m >= 200
  )

# Validate the release-month lookup before plotting.

if (any(is.na(particle_displacement$release_month))) {

  stop(
    "At least one release_id could not be matched to June, July or August."
  )
}

# Define analytical depth strata

particle_displacement <- particle_displacement %>%

  mutate(
    depth_stratum = case_when(
      source_depth_m >= 200 &
        source_depth_m < 700 ~
        "200-700 m",
      source_depth_m >= 700 &
        source_depth_m < 1500 ~
        "700-1500 m",
      source_depth_m >= 1500 &
        source_depth_m < 2500 ~
        "1500-2500 m",
      source_depth_m >= 2500 ~
        ">2500 m"
    ),
    depth_stratum = factor(
      depth_stratum,
      levels = c(
        "200-700 m",
        "700-1500 m",
        "1500-2500 m",
        ">2500 m"
      )
    )
  )

# Mean maximum displacement at each release depth

# The mean is used here as the summary shown by the black points
# and connecting line.

mean_displacement_by_depth <- particle_displacement %>%

  group_by(
    release_month,
    source_depth_m,
    depth_stratum
  ) %>%

  summarise(
    mean_maximum_horizontal_displacement_km = mean(
      maximum_horizontal_displacement_km,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  arrange(
    release_month,
    source_depth_m
  )

# Spearman depth-displacement relationship

# These values are saved for reporting in the proposal text but are
# intentionally not printed on the figure.

displacement_correlations <- particle_displacement %>%

  group_by(
    release_month
  ) %>%

  summarise(
    N_particles = n(),
    spearman_rho = cor(
      source_depth_m,
      maximum_horizontal_displacement_km,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

write.csv(
  displacement_correlations,
  file.path(
    figure_output_dir,
    "taupo_depth_displacement_spearman.csv"
  ),
  row.names = FALSE
)

# Background depth bands

# The limits extend slightly beyond the shallowest and deepest release
# points so the plotting region has visual breathing room.

depth_bands <- tibble(
  depth_stratum = factor(
    c(
      "200-700 m",
      "700-1500 m",
      "1500-2500 m",
      ">2500 m"
    ),
    levels = c(
      "200-700 m",
      "700-1500 m",
      "1500-2500 m",
      ">2500 m"
    )
  ),
  ymin = c(
    150,
    700,
    1500,
    2500
  ),
  ymax = c(
    700,
    1500,
    2500,
    4150
  )
)

# Soft sequential blue palette

# The background is deliberately pale so the strata are immediately
# identifiable without competing with the particle observations.

depth_band_colours <- c(
  "200-700 m"   = "#F7FAFC",
  "700-1500 m"  = "#E5EFF5",
  "1500-2500 m" = "#CDDFEA",
  ">2500 m"     = "#B3CDDC"
)

# Y-axis release-depth breaks

depth_breaks <- c(
  200,
  500,
  1000,
  1500,
  2000,
  2500,
  3000,
  3500,
  4000
)

# X-axis range

# Round the largest observed displacement upward to the nearest 100 km.
# This prevents clipping while keeping the scale consistent among months.

maximum_plot_displacement <- max(
  particle_displacement$maximum_horizontal_displacement_km,
  na.rm = TRUE
)

x_upper <- ceiling(
  maximum_plot_displacement / 100
) * 100

# Figure 2

figure_2 <- ggplot() +

  # Depth-stratum background

  geom_rect(
    data = depth_bands,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax,
      fill = depth_stratum
    ),
    inherit.aes = FALSE,
    colour = NA
  ) +

  # Individual particles
  #
  # Point size and opacity are chosen so individual observations remain visible
  # when the figure is reduced to proposal-page size.
  #
  # Transparency still allows overlapping particles to remain visible.

  geom_point(
    data = particle_displacement,
    aes(
      x = maximum_horizontal_displacement_km,
      y = source_depth_m
    ),
    colour = "#344A58",
    alpha = 0.62,
    size = 1.20
  ) +

  # Mean maximum displacement
  #
  # The line summarises the systematic change in dispersal with
  # release depth while retaining the complete particle distribution
  # underneath.

  geom_line(
    data = mean_displacement_by_depth,
    aes(
      x = mean_maximum_horizontal_displacement_km,
      y = source_depth_m,
      group = release_month
    ),
    colour = "black",
    linewidth = 0.90,
    lineend = "round",
    linejoin = "round"
  ) +

  # Mean points

  geom_point(
    data = mean_displacement_by_depth,
    aes(
      x = mean_maximum_horizontal_displacement_km,
      y = source_depth_m
    ),
    shape = 16,
    colour = "black",
    size = 2.90
  ) +

  # One panel per release month

  facet_wrap(
    ~ release_month,
    nrow = 1
  ) +

  # Depth-band legend

  scale_fill_manual(
    values = depth_band_colours,
    name = "Depth stratum",
    breaks = c(
      "200-700 m",
      "700-1500 m",
      "1500-2500 m",
      ">2500 m"
    ),
    labels = c(
      "200–700 m",
      "700–1500 m",
      "1500–2500 m",
      "≥ 2500 m"
    )
  ) +

  # Horizontal displacement axis

  scale_x_continuous(
    breaks = seq(
      0,
      x_upper,
      by = 500
    ),
    limits = c(
      0,
      x_upper
    ),
    expand = c(
      0,
      0
    )
  ) +

  # Reverse release depth so shallow water appears at the top and
  # deeper releases appear lower in the figure.

  scale_y_reverse(
    breaks = depth_breaks,
    limits = c(
      4150,
      150
    ),
    expand = c(
      0,
      0
    )
  ) +

  # Shared axis labels

  labs(
    x = "Maximum horizontal displacement (km)",
    y = "Release depth (m)"
  ) +

  # Publication-style theme

  theme_bw(
    base_size = 10.5
  ) +

  theme(
    # Thin full border around each panel.
    #
    # The separate axis lines are removed below so the left and bottom
    # borders are not drawn twice.
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    # No internal plotting grid.
    #
    # The depth bands already provide sufficient visual structure.
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Text
    text = element_text(
      colour = "black"
    ),
    # Plain axis titles keep them visually subordinate to the month
    # headings and data.
    axis.title = element_text(
      face = "plain",
      colour = "black",
      size = 11
    ),
    axis.text = element_text(
      colour = "black",
      size = 9.5
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.40
    ),
    axis.ticks.length = unit(
      0.12,
      "cm"
    ),
    # Remove default grey facet strips.
    strip.background = element_blank(),
    # Larger month labels without surrounding boxes.
    strip.text = element_text(
      face = "bold",
      colour = "black",
      size = 13.25,
      margin = margin(
        t = 1,
        b = 6
      )
    ),
    # Keep the three panels visually close enough to read as one figure
    # while preserving clear white gutters between months.
    panel.spacing.x = unit(
      0.30,
      "lines"
    ),
    # Compact horizontal legend.
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(
      face = "plain",
      size = 9.5
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.key.width = unit(
      0.72,
      "cm"
    ),
    legend.key.height = unit(
      0.30,
      "cm"
    ),
    legend.spacing.x = unit(
      0.16,
      "cm"
    ),
    legend.margin = margin(
      t = 2,
      r = 0,
      b = 0,
      l = 0
    ),
    plot.margin = margin(
      t = 8,
      r = 6,
      b = 6,
      l = 6
    )
  )

# Save 

figure_2_file <- file.path(
  figure_output_dir,
  "Figure_2_depth_dependent_dispersal.pdf"
)

ggsave(
  filename = figure_2_file,
  plot = figure_2,
  device = cairo_pdf,
  width = 8.7,
  height = 5.05,
  units = "in"
)

