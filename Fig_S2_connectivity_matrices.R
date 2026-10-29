# Supplementary Figure S2: all nine Taupo Bank connectivity matrices
#
# Exact reference design:
# - rows = PLD = 30, 60, 90 d
# - columns = PC = 14, 21, 28 d
# - each scenario contains T22, T21, T20 and L11 destination groups
# - four destination depth columns beneath each seamount
# - four Taupo source-depth rows
# - shared 0-10% greyscale connectivity legend

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# Inputs

analysis_file <- file.path(
  "Outputs",
  "connectivity_analysis.rds"
)

if (!file.exists(analysis_file)) {
  stop(
    "Could not find Outputs/connectivity_analysis.rds. ",
    "Run the current validated 04 network-analysis script first."
  )
}

connectivity_analysis <- readRDS(
  analysis_file
)

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

depth_levels <- c(
  "200-700 m",
  "700-1500 m",
  "1500-2500 m",
  ">2500 m"
)

destination_order <- c(
  "T22",
  "T21",
  "T20",
  "L11"
)

# Theme retained from the reference figure

theme_publication <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.text = element_text(
        face = "bold",
        colour = "black"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      legend.text = element_text(
        colour = "black"
      ),
      legend.key = element_blank(),
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
}

# Pooled source-normalised connectivity

pooled_connectivity <- connectivity_analysis$connectivity %>%
  filter(
    analysis_period == "pooled",
    period_type == "pooled",
    source_node == "T22"
  )

# Include every destination reached in at least one biological scenario.

all_active_destination_ids <- pooled_connectivity %>%
  filter(
    pld_days %in% c(30, 60, 90),
    precompetency_days %in% c(14, 21, 28),
    connectivity > 0
  ) %>%
  pull(destination_node) %>%
  as.character() %>%
  unique()

all_active_destinations <- destination_order[
  destination_order %in% all_active_destination_ids
]

extra_destinations <- setdiff(
  all_active_destination_ids,
  destination_order
)

all_active_destinations <- c(
  all_active_destinations,
  extra_destinations
)

# Complete scenario matrix

supp_edges <- pooled_connectivity %>%
  filter(
    pld_days %in% c(30, 60, 90),
    precompetency_days %in% c(14, 21, 28),
    destination_node %in% all_active_destinations
  ) %>%
  select(
    pld_days,
    precompetency_days,
    source_depth_stratum,
    destination_node,
    destination_depth_stratum,
    connectivity
  )

supp_matrix_data <- tidyr::expand_grid(
  pld_days = c(30, 60, 90),
  precompetency_days = c(14, 21, 28),
  source_depth_stratum = depth_levels,
  destination_node = all_active_destinations,
  destination_depth_stratum = depth_levels
) %>%
  left_join(
    supp_edges,
    by = c(
      "pld_days",
      "precompetency_days",
      "source_depth_stratum",
      "destination_node",
      "destination_depth_stratum"
    )
  ) %>%
  mutate(
    connectivity = coalesce(
      connectivity,
      0
    ),

    pld_label = factor(
      paste0(
        "PLD = ",
        pld_days,
        " d"
      ),
      levels = c(
        "PLD = 30 d",
        "PLD = 60 d",
        "PLD = 90 d"
      )
    ),

    pc_label = factor(
      paste0(
        "PC = ",
        precompetency_days,
        " d"
      ),
      levels = c(
        "PC = 14 d",
        "PC = 21 d",
        "PC = 28 d"
      )
    ),

    source_depth_stratum = factor(
      source_depth_stratum,
      levels = rev(depth_levels)
    ),

    destination_node = factor(
      destination_node,
      levels = all_active_destinations
    ),

    destination_depth_stratum = factor(
      destination_depth_stratum,
      levels = depth_levels
    )
  )

# Destination-depth positions beneath each seamount heading

seamount_positions <- tibble(
  destination_node = factor(
    all_active_destinations,
    levels = all_active_destinations
  ),
  site_index = seq_along(
    all_active_destinations
  )
) %>%
  mutate(
    base_x =
      (site_index - 1) * 4.8
  )

depth_positions <- tibble(
  destination_depth_stratum = factor(
    depth_levels,
    levels = depth_levels
  ),
  depth_index = seq_along(
    depth_levels
  ),
  depth_label = c(
    "200-700",
    "700-1500",
    "1500-2500",
    ">2500"
  )
)

supp_matrix_data <- supp_matrix_data %>%
  left_join(
    seamount_positions,
    by = "destination_node"
  ) %>%
  left_join(
    depth_positions,
    by = "destination_depth_stratum"
  ) %>%
  mutate(
    x_pos = base_x + depth_index,

    y_pos = case_when(
      source_depth_stratum == ">2500 m" ~ 1,
      source_depth_stratum == "1500-2500 m" ~ 2,
      source_depth_stratum == "700-1500 m" ~ 3,
      source_depth_stratum == "200-700 m" ~ 4
    ),

    connectivity_label = if_else(
      connectivity > 0,
      sprintf(
        "%.1f%%",
        100 * connectivity
      ),
      ""
    ),

    label_colour = if_else(
      connectivity >= 0.055,
      "white",
      "black"
    )
  )

supp_x_axis <- tidyr::expand_grid(
  destination_node = factor(
    all_active_destinations,
    levels = all_active_destinations
  ),
  destination_depth_stratum = factor(
    depth_levels,
    levels = depth_levels
  )
) %>%
  left_join(
    seamount_positions,
    by = "destination_node"
  ) %>%
  left_join(
    depth_positions,
    by = "destination_depth_stratum"
  ) %>%
  mutate(
    x_pos = base_x + depth_index
  ) %>%
  arrange(
    site_index,
    depth_index
  )

supp_group_labels <- tidyr::expand_grid(
  pld_label = factor(
    c(
      "PLD = 30 d",
      "PLD = 60 d",
      "PLD = 90 d"
    ),
    levels = c(
      "PLD = 30 d",
      "PLD = 60 d",
      "PLD = 90 d"
    )
  ),

  pc_label = factor(
    c(
      "PC = 14 d",
      "PC = 21 d",
      "PC = 28 d"
    ),
    levels = c(
      "PC = 14 d",
      "PC = 21 d",
      "PC = 28 d"
    )
  ),

  destination_node = factor(
    all_active_destinations,
    levels = all_active_destinations
  )
) %>%
  left_join(
    seamount_positions,
    by = "destination_node"
  ) %>%
  mutate(
    x_mid = base_x + 2.5,
    y_pos = 4.62
  )

separator_positions <-
  seamount_positions$base_x[-1] + 0.1

# Figure

figure_S2 <- ggplot(
  supp_matrix_data,
  aes(
    x = x_pos,
    y = y_pos,
    fill = connectivity
  )
) +
  geom_tile(
    width = 0.95,
    height = 0.95,
    colour = "grey78",
    linewidth = 0.32
  ) +

  geom_text(
    aes(
      label = connectivity_label,
      colour = label_colour
    ),
    size = 2.45,
    show.legend = FALSE
  ) +

  geom_vline(
    xintercept = separator_positions,
    colour = "grey25",
    linewidth = 0.6
  ) +

  geom_text(
    data = supp_group_labels,
    aes(
      x = x_mid,
      y = y_pos,
      label = destination_node
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.0
  ) +

  facet_grid(
    rows = vars(pld_label),
    cols = vars(pc_label),
    switch = "y"
  ) +

  scale_colour_identity() +

  scale_x_continuous(
    breaks = supp_x_axis$x_pos,
    labels = supp_x_axis$depth_label,
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +

  scale_y_continuous(
    breaks = c(
      1,
      2,
      3,
      4
    ),
    labels = c(
      ">2500 m",
      "1500-2500 m",
      "700-1500 m",
      "200-700 m"
    ),
    limits = c(
      0.5,
      4.72
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      0.10
    ),
    breaks = c(
      0,
      0.025,
      0.05,
      0.075,
      0.10
    ),
    labels = function(x) {
      paste0(
        round(
          100 * x,
          1
        ),
        "%"
      )
    }
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  labs(
    x = "Destination seamount and depth stratum",
    y = "Taupo Bank source depth stratum",
    fill = "Potential connectivity"
  ) +

  theme_publication(
    base_size = 10
  ) +

  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),

    panel.grid = element_blank(),

    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.45
    ),

    strip.text.x = element_text(
      face = "bold",
      size = 9
    ),

    strip.text.y.left = element_text(
      face = "bold",
      size = 9,
      angle = 90
    ),

    strip.placement = "outside",

    axis.text.x = element_text(
      angle = 48,
      hjust = 1,
      size = 6.7
    ),

    axis.text.y = element_text(
      size = 8
    ),

    panel.spacing = grid::unit(
      0.65,
      "lines"
    ),

    legend.position = "bottom",
    legend.direction = "horizontal",

    plot.margin = margin(
      12,
      16,
      24,
      14
    )
  ) +

  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(
        7.0,
        "cm"
      ),
      barheight = grid::unit(
        0.28,
        "cm"
      )
    )
  )

# Save at the same physical dimensions as the reference figure

ggsave(
  file.path(
    supplementary_dir,
    "Figure_S2_all_9_connectivity_matrices.pdf"
  ),
  figure_S2,
  width = 15.0,
  height = 10.5,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(
    supplementary_dir,
    "Figure_S2_all_9_connectivity_matrices.pdf"
  ),
  figure_S2,
  width = 15.0,
  height = 10.5,
  units = "in",
  dpi = 350,
  limitsize = FALSE
)

print(figure_S2)


