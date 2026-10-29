# Supplementary Figure S1: Taupo Bank particle dispersal
#
# Exact reference design:
# - rows = June, July, August
# - columns = source depth strata
# - colour = proportion of released particles entering each 0.10 degree grid cell
# - GEBCO contours and study sites retained for geographic context
#
# This is a physical dispersal / visitation figure, not a habitat-encounter
# or connectivity calculation.

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(terra)
library(scales)

# Inputs

setup_candidates <- c(
  "00_setup_and_inputs_(7).R",
  "00_setup_and_inputs_.R",
  "00_setup_and_inputs.R"
)

setup_file <- setup_candidates[file.exists(setup_candidates)][1]

if (is.na(setup_file)) {
  stop(
    "Could not find the current 00 setup script. Expected one of: ",
    paste(setup_candidates, collapse = ", ")
  )
}

source(setup_file)

particle_candidates <- c(
  file.path("Outputs", "particle_trajectories.rds"),
  "particle_trajectories.rds"
)

particle_file <- particle_candidates[file.exists(particle_candidates)][1]

if (is.na(particle_file)) {
  stop(
    "Could not find particle_trajectories.rds. Run the current 03 script first."
  )
}

particle_trajectories <- readRDS(particle_file)

required_columns <- c(
  "particle_id",
  "release_id",
  "age_hours",
  "lon",
  "lat",
  "source_depth_stratum"
)

missing_columns <- setdiff(
  required_columns,
  names(particle_trajectories)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "particle_trajectories.rds is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Figure settings
#
# Keep these unchanged to reproduce the reference figure physically.

min_age_days <- 1
max_age_days <- 90
grid_size_deg <- 0.10
map_buffer_deg <- 0.30

month_levels <- c(
  "June",
  "July",
  "August"
)

depth_levels <- c(
  "200-700 m",
  "700-1500 m",
  "1500-2500 m",
  ">2500 m"
)

# Prepare trajectories

trajectory_data <- particle_trajectories %>%
  mutate(
    age_days = age_hours / 24,

    release_month = case_when(
      release_id == 1 ~ "June",
      release_id == 2 ~ "July",
      release_id == 3 ~ "August",
      TRUE ~ NA_character_
    ),

    release_month = factor(
      release_month,
      levels = month_levels
    ),

    source_depth_stratum = factor(
      source_depth_stratum,
      levels = depth_levels
    ),

    particle_key = paste(
      release_id,
      particle_id,
      sep = "_"
    )
  ) %>%
  filter(
    !is.na(release_month),
    !is.na(source_depth_stratum),
    is.finite(lon),
    is.finite(lat),
    is.finite(age_days),
    age_days >= min_age_days,
    age_days <= max_age_days
  )

if (nrow(trajectory_data) == 0) {
  stop("No trajectory positions remained after filtering.")
}

# Number of released particles represented in each panel

panel_particle_counts <- trajectory_data %>%
  distinct(
    release_month,
    source_depth_stratum,
    particle_key
  ) %>%
  count(
    release_month,
    source_depth_stratum,
    name = "n_particles_panel"
  )

# Spatial visitation grid
#
# Each particle contributes at most once to a grid cell within its panel.

particle_cell_visits <- trajectory_data %>%
  mutate(
    lon_bin =
      floor(lon / grid_size_deg) * grid_size_deg +
      grid_size_deg / 2,

    lat_bin =
      floor(lat / grid_size_deg) * grid_size_deg +
      grid_size_deg / 2
  ) %>%
  distinct(
    release_month,
    source_depth_stratum,
    particle_key,
    lon_bin,
    lat_bin
  )

visitation_grid <- particle_cell_visits %>%
  count(
    release_month,
    source_depth_stratum,
    lon_bin,
    lat_bin,
    name = "n_particles_cell"
  ) %>%
  left_join(
    panel_particle_counts,
    by = c(
      "release_month",
      "source_depth_stratum"
    )
  ) %>%
  mutate(
    visit_fraction =
      n_particles_cell / n_particles_panel
  )

# Retain all 12 month x depth panels

panel_grid <- tidyr::expand_grid(
  release_month = factor(
    month_levels,
    levels = month_levels
  ),
  source_depth_stratum = factor(
    depth_levels,
    levels = depth_levels
  )
)

# Map extent

x_limits <- range(
  trajectory_data$lon,
  na.rm = TRUE
) + c(-map_buffer_deg, map_buffer_deg)

y_limits <- range(
  trajectory_data$lat,
  na.rm = TRUE
) + c(-map_buffer_deg, map_buffer_deg)

x_limits[1] <- max(x_limits[1], terra::xmin(gebco))
x_limits[2] <- min(x_limits[2], terra::xmax(gebco))
y_limits[1] <- max(y_limits[1], terra::ymin(gebco))
y_limits[2] <- min(y_limits[2], terra::ymax(gebco))

# Bathymetry contours

gebco_crop <- terra::crop(
  gebco,
  terra::ext(
    x_limits[1],
    x_limits[2],
    y_limits[1],
    y_limits[2]
  )
)

gebco_plot <- terra::aggregate(
  gebco_crop,
  fact = 4,
  fun = mean,
  na.rm = TRUE
)

gebco_df <- as.data.frame(
  gebco_plot,
  xy = TRUE,
  na.rm = FALSE
)

names(gebco_df)[3] <- "elevation"

background_contours <- c(
  -4000,
  -3000,
  -2500,
  -1500,
  -700,
  -200
)

# Taupo source location and surrounding study sites

taupo_node <- candidate_nodes %>%
  filter(node_id == "T22") %>%
  slice(1)

if (nrow(taupo_node) == 0) {
  stop("T22 could not be found in candidate_nodes.")
}

context_nodes <- candidate_nodes %>%
  filter(
    lon >= x_limits[1],
    lon <= x_limits[2],
    lat >= y_limits[1],
    lat <= y_limits[2],
    node_id != "T22"
  )

# Shared colour scale

max_visit_fraction <- max(
  visitation_grid$visit_fraction,
  na.rm = TRUE
)

candidate_breaks <- c(
  0.01,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  1.00
)

fill_breaks <- candidate_breaks[
  candidate_breaks <= max_visit_fraction
]

if (length(fill_breaks) < 2) {
  fill_breaks <- pretty(
    c(0, max_visit_fraction),
    n = 4
  )

  fill_breaks <- fill_breaks[
    fill_breaks > 0 &
      fill_breaks <= max_visit_fraction
  ]
}

# Figure

figure_S1 <- ggplot() +
  geom_blank(
    data = panel_grid,
    aes(
      x = x_limits[1],
      y = y_limits[1]
    )
  ) +

  geom_tile(
    data = visitation_grid,
    aes(
      x = lon_bin,
      y = lat_bin,
      fill = visit_fraction
    ),
    width = grid_size_deg,
    height = grid_size_deg
  ) +

  geom_contour(
    data = gebco_df,
    aes(
      x = x,
      y = y,
      z = elevation
    ),
    breaks = background_contours,
    colour = "grey60",
    linewidth = 0.28,
    alpha = 0.65
  ) +

  geom_point(
    data = context_nodes,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    size = 1.6,
    stroke = 0.45,
    fill = "white",
    colour = "black"
  ) +

  geom_point(
    data = taupo_node,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    size = 3.2,
    stroke = 0.8,
    fill = "black",
    colour = "white"
  ) +

  geom_text(
    data = taupo_node,
    aes(
      x = lon,
      y = lat,
      label = "T22"
    ),
    nudge_x = 0.12,
    nudge_y = -0.10,
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold",
    colour = "black"
  ) +

  facet_grid(
    release_month ~ source_depth_stratum,
    drop = FALSE
  ) +

  scale_fill_viridis_c(
    option = "C",
    direction = 1,
    trans = "sqrt",
    limits = c(
      0,
      max_visit_fraction
    ),
    breaks = fill_breaks,
    labels = scales::label_percent(
      accuracy = 1
    ),
    name = "Particles entering\ngrid cell"
  ) +

  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE
  ) +

  labs(
    x = "Longitude",
    y = "Latitude"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.title = element_text(
      face = "bold"
    ),

    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.4
    ),

    strip.text = element_text(
      face = "bold",
      size = 10
    ),

    strip.text.y = element_text(
      angle = 0
    ),

    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),

    panel.spacing = grid::unit(
      0.45,
      "lines"
    ),

    legend.position = "right",

    legend.title = element_text(
      face = "bold"
    ),

    legend.key.height = grid::unit(
      2.5,
      "cm"
    )
  )

# Save

supplementary_dir <- file.path(
  "Outputs",
  "Proposal_figures",
  "Supplementary"
)

dir.create(
  supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    supplementary_dir,
    "Figure_S1_particle_dispersal.pdf"
  ),
  figure_S1,
  width = 13.5,
  height = 9.5,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    supplementary_dir,
    "Figure_S1_particle_dispersal.png"
  ),
  figure_S1,
  width = 13.5,
  height = 9.5,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

print(figure_S1)

message(
  "Saved Supplementary Figure S1 to: ",
  supplementary_dir
)
