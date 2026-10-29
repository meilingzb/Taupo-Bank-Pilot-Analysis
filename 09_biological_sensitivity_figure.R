# 09 Figure 5: biological-scenario sensitivity
#
# Purpose:
# Plot how Taupo Bank pilot outcomes vary across the 3 x 3 combination of
# PLD (30, 60, 90 d) and pre-competency (14, 21, 28 d).

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
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

analysis <- readRDS(
  analysis_file
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

# Pooled biological scenarios

pooled_connectivity <- analysis$connectivity %>%
  filter(
    analysis_period == "pooled",
    period_type == "pooled"
  )

# Habitat-encounter probability
#
# N_released is repeated for destination rows in the connectivity table.
# Retain one denominator per source site-depth node before summing.

scenario_releases <- pooled_connectivity %>%
  distinct(
    pld_days,
    precompetency_days,
    source_node,
    source_depth_stratum,
    N_released
  ) %>%
  group_by(
    pld_days,
    precompetency_days
  ) %>%
  summarise(
    N_released =
      sum(
        N_released
      ),
    .groups = "drop"
  )

scenario_encounters <- pooled_connectivity %>%
  group_by(
    pld_days,
    precompetency_days
  ) %>%
  summarise(
    N_habitat_encounters =
      sum(
        N_settled
      ),
    .groups = "drop"
  )

# Cross-depth connectivity fraction
#
# F_cross is read directly from 04, where it is calculated from source-normalised
# connectivity totals rather than raw counts of successful encounters.

scenario_cross_depth <- analysis$cross_depth_fraction %>%
  filter(
    analysis_period == "pooled",
    period_type == "pooled"
  ) %>%
  select(
    pld_days,
    precompetency_days,
    total_connectivity,
    cross_depth_connectivity,
    same_seamount_cross_depth_connectivity,
    between_seamount_cross_depth_connectivity,
    F_cross
  )

# Combined scenario table

taupo_scenario_summary <- scenario_releases %>%
  left_join(
    scenario_encounters,
    by = c(
      "pld_days",
      "precompetency_days"
    )
  ) %>%
  left_join(
    scenario_cross_depth,
    by = c(
      "pld_days",
      "precompetency_days"
    )
  ) %>%
  mutate(
    habitat_encounter_probability =
      N_habitat_encounters /
      N_released
  ) %>%
  arrange(
    pld_days,
    precompetency_days
  )

write.csv(
  taupo_scenario_summary,
  file.path(
    output_dir,
    "taupo_biological_sensitivity.csv"
  ),
  row.names = FALSE
)

# Factor order
#
# With factor levels 14, 21, 28, ggplot places 14 at the bottom and 28 at the
# top, matching the reference Figure 5.

plot_data <- taupo_scenario_summary %>%
  mutate(
    pld_days =
      factor(
        pld_days,
        levels = c(
          30,
          60,
          90
        )
      ),
    precompetency_days =
      factor(
        precompetency_days,
        levels = c(
          14,
          21,
          28
        )
      )
  )

# Colour scales
#
# These blue anchors were recovered from the reference PDF so that the rebuilt
# figure retains the same visual palette. The left and right panels use their
# original independent ranges.

encounter_anchor_values <- c(
  0.000,
  0.010,
  0.038,
  0.042,
  0.052,
  0.061,
  0.062,
  0.076,
  0.081,
  0.094,
  0.100
)

encounter_anchor_colours <- c(
  "#EEF3F6",
  "#E5ECF0",
  "#C7D5E0",
  "#C2D2DE",
  "#B4C8D7",
  "#A3BCCE",
  "#A1BACD",
  "#88A8C0",
  "#7FA1BB",
  "#6990AE",
  "#5D87A8"
)

cross_anchor_values <- c(
  0.000,
  0.067,
  0.096,
  0.116,
  0.122,
  0.131,
  0.135,
  0.139,
  0.145,
  0.150
)

cross_anchor_colours <- c(
  "#EEF3F6",
  "#BFD0DC",
  "#9EB8CC",
  "#85A6BF",
  "#7EA0BA",
  "#7398B5",
  "#6F95B2",
  "#6A91AF",
  "#628BAB",
  "#5D87A8"
)

# Common physical styling

theme_sensitivity_heatmap <- theme_bw(
  base_size = 11
) +
  theme(
    panel.grid =
      element_blank(),
    panel.border =
      element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.45
      ),
    axis.title =
      element_text(
        colour = "black",
        face = "plain",
        size = 11
      ),
    axis.text =
      element_text(
        colour = "black",
        size = 9.5
      ),
    axis.ticks =
      element_blank(),
    plot.title =
      element_text(
        colour = "black",
        face = "bold",
        size = 11.5,
        hjust = 0.5,
        margin = margin(
          b = 6
        )
      ),
    legend.position =
      "bottom",
    legend.direction =
      "horizontal",
    legend.title =
      element_text(
        colour = "black",
        face = "plain",
        size = 9.5,
        hjust = 0.5
      ),
    legend.text =
      element_text(
        colour = "black",
        size = 8.5
      ),
    legend.margin =
      margin(
        t = 4,
        r = 0,
        b = 0,
        l = 0
      ),
    legend.box.margin =
      margin(
        0,
        0,
        0,
        0
      ),
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
  )

# Panel a: habitat-encounter probability

figure_5_a <- ggplot(
  plot_data,
  aes(
    x = pld_days,
    y = precompetency_days,
    fill =
      habitat_encounter_probability
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.90
  ) +

  geom_text(
    aes(
      label =
        sprintf(
          "%.1f%%",
          100 *
            habitat_encounter_probability
        )
    ),
    colour = "black",
    fontface = "bold",
    size = 4.20
  ) +

  scale_fill_gradientn(
    colours =
      encounter_anchor_colours,
    values =
      scales::rescale(
        encounter_anchor_values,
        from = c(
          0,
          0.10
        )
      ),
    limits = c(
      0,
      0.10
    ),
    breaks = c(
      0,
      0.025,
      0.050,
      0.075,
      0.100
    ),
    labels =
      label_percent(
        accuracy = 0.1
      ),
    oob =
      scales::squish,
    name =
      "Habitat-encounter probability"
  ) +

  coord_equal(
    expand = FALSE
  ) +

  labs(
    title =
      "a) Habitat-encounter probability",
    x =
      "Pelagic larval duration (days)",
    y =
      "Pre-competency (days)"
  ) +

  theme_sensitivity_heatmap +

  guides(
    fill =
      guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(
          4.20,
          "cm"
        ),
        barheight = unit(
          0.30,
          "cm"
        ),
        ticks = FALSE,
        frame.colour = NA
      )
  )

# Panel b: source-normalised cross-depth connectivity fraction

figure_5_b <- ggplot(
  plot_data,
  aes(
    x = pld_days,
    y = precompetency_days,
    fill = F_cross
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.90
  ) +

  geom_text(
    aes(
      label =
        sprintf(
          "%.1f%%",
          100 *
            F_cross
        )
    ),
    colour = "black",
    fontface = "bold",
    size = 4.20
  ) +

  scale_fill_gradientn(
    colours =
      cross_anchor_colours,
    values =
      scales::rescale(
        cross_anchor_values,
        from = c(
          0,
          0.15
        )
      ),
    limits = c(
      0,
      0.15
    ),
    breaks = c(
      0,
      0.05,
      0.10,
      0.15
    ),
    labels =
      label_percent(
        accuracy = 0.1
      ),
    oob =
      scales::squish,
    name =
      expression(
        F[cross]
      )
  ) +

  coord_equal(
    expand = FALSE
  ) +

  labs(
    title =
      "b) Cross-depth fraction of connectivity",
    x =
      "Pelagic larval duration (days)",
    y =
      "Pre-competency (days)"
  ) +

  theme_sensitivity_heatmap +

  guides(
    fill =
      guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(
          3.80,
          "cm"
        ),
        barheight = unit(
          0.30,
          "cm"
        ),
        ticks = FALSE,
        frame.colour = NA
      )
  )

# Assemble

figure_5 <- figure_5_a +
  figure_5_b +
  plot_layout(
    widths = c(
      1,
      1
    )
  )

# Save
#
# The dimensions match the existing proposal PDF (~626 x 331 pt).

figure_5_file <- file.path(
  output_dir,
  "Figure_5_biological_sensitivity.pdf"
)

ggsave(
  filename = figure_5_file,
  plot = figure_5,
  device = cairo_pdf,
  width = 8.694,
  height = 4.597,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_5_biological_sensitivity.png"
  ),
  plot = figure_5,
  width = 8.694,
  height = 4.597,
  units = "in",
  dpi = 600,
  bg = "white"
)

