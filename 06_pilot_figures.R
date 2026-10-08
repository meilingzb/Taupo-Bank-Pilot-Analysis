# Taupo Bank Analysis pt. 7 :)
# Final static figures

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(terra)
library(scales)
library(patchwork)

source("00_setup_and_inputs_.R")

dir.create(
  "figures",
  showWarnings = FALSE
)

analysis <- readRDS(
  "connectivity_analysis.rds"
)

connectivity <- analysis$connectivity
network_roles <- analysis$network_roles
depth_levels <- analysis$depth_levels

pooled_connectivity <- connectivity %>%
  filter(
    analysis_period == "pooled"
  )


# Plotting data

plot_nodes <- candidate_nodes %>%
  mutate(
    plot_lon = lon,
    plot_lat = lat
  )

depth_map_edges <- pooled_connectivity %>%
  filter(
    source_node == "T22",
    destination_node != source_node,
    connectivity > 0
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
  ) %>%
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


# Colours

depth_colours <- c(
  "200-700 m" = "#0072B2",
  "700-1500 m" = "#009E73",
  "1500-2500 m" = "#E69F00",
  ">2500 m" = "#7A3E9D"
)

curve_values <- c(
  "200-700 m" = -0.20,
  "700-1500 m" = -0.07,
  "1500-2500 m" = 0.07,
  ">2500 m" = 0.20
)

pooled_network_edges <- pooled_connectivity %>%
  filter(
    source_node == "T22",
    connectivity > 0
  )

max_pooled_connectivity <- max(
  pooled_network_edges$connectivity,
  na.rm = TRUE
)


# Common map extent

connected_ids <- unique(
  c(
    "T22",
    depth_map_edges$destination_node
  )
)

extent_nodes <- plot_nodes %>%
  filter(
    node_id %in% connected_ids
  )

x_limits <- range(
  extent_nodes$plot_lon,
  na.rm = TRUE
) + c(-0.6, 0.6)

y_limits <- range(
  extent_nodes$plot_lat,
  na.rm = TRUE
) + c(-0.8, 0.8)

background_nodes <- plot_nodes %>%
  filter(
    plot_lon >= x_limits[1],
    plot_lon <= x_limits[2],
    plot_lat >= y_limits[1],
    plot_lat <= y_limits[2]
  )


# GEBCO contours

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
  fact = 3,
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
  -3500,
  -2500,
  -1500,
  -700,
  -200
)



# Figure A
# Geographic connectivity by source depth
# External arrows plus depth-resolved retention around Taupo Bank

depth_legend_dummy <- tibble(
  destination_depth_stratum = factor(
    depth_levels,
    levels = depth_levels
  ),
  x = NA_real_,
  y = NA_real_,
  xend = NA_real_,
  yend = NA_real_
)

taupo_point <- plot_nodes %>%
  filter(
    node_id == "T22"
  )

taupo_lon <- taupo_point$plot_lon[1]
taupo_lat <- taupo_point$plot_lat[1]

map_loop_radius <- c(
  "200-700 m" = 0.08,
  "700-1500 m" = 0.11,
  "1500-2500 m" = 0.14,
  ">2500 m" = 0.17
)

make_map_retention_loops <- function(edges) {
  
  if (nrow(edges) == 0) {
    return(
      tibble(
        loop_id = character(),
        source_depth_stratum = factor(
          character(),
          levels = depth_levels
        ),
        destination_depth_stratum = factor(
          character(),
          levels = depth_levels
        ),
        connectivity = numeric(),
        x = numeric(),
        y = numeric()
      )
    )
  }
  
  bind_rows(
    lapply(
      seq_len(
        nrow(edges)
      ),
      function(i) {
        
        edge <- edges[i, ]
        
        depth_current <- as.character(
          edge$destination_depth_stratum
        )
        
        theta <- seq(
          pi / 4,
          2 * pi + pi / 4 - 0.25,
          length.out = 70
        )
        
        radius_current <- map_loop_radius[
          depth_current
        ]
        
        tibble(
          loop_id = paste(
            edge$pld_days,
            edge$source_depth_stratum,
            edge$destination_depth_stratum,
            sep = "_"
          ),
          source_depth_stratum = edge$source_depth_stratum,
          destination_depth_stratum = edge$destination_depth_stratum,
          connectivity = edge$connectivity,
          x = taupo_lon +
            radius_current * cos(theta),
          y = taupo_lat +
            radius_current * sin(theta)
        )
      }
    )
  )
}


plot_geographic_connectivity <- function(pld_current) {
  
  edges_current <- depth_map_edges %>%
    filter(
      pld_days == pld_current
    )
  
  retention_current <- pooled_network_edges %>%
    filter(
      pld_days == pld_current,
      destination_node == "T22"
    ) %>%
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
  
  retention_loops <- make_map_retention_loops(
    retention_current
  )
  
  connected_current <- unique(
    c(
      "T22",
      edges_current$destination_node
    )
  )
  
  connected_nodes <- plot_nodes %>%
    filter(
      node_id %in% connected_current
    )
  
  p <- ggplot() +
    
    # Keeps all four destination depths in every legend
    
    geom_segment(
      data = depth_legend_dummy,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = destination_depth_stratum
      ),
      linewidth = 1.5,
      alpha = 0,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.12,
          "cm"
        ),
        type = "closed"
      ),
      show.legend = TRUE
    ) +
    
    geom_contour(
      data = gebco_df,
      aes(
        x = x,
        y = y,
        z = elevation
      ),
      breaks = background_contours,
      colour = "grey82",
      linewidth = 0.3
    ) +
    
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
      colour = "grey70"
    )
  
  for (depth_current in depth_levels) {
    
    edge_subset <- edges_current %>%
      filter(
        destination_depth_stratum == depth_current
      )
    
    if (nrow(edge_subset) > 0) {
      
      p <- p +
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
          curvature = curve_values[[depth_current]],
          alpha = 0.92,
          lineend = "round",
          arrow = grid::arrow(
            length = grid::unit(
              0.12,
              "cm"
            ),
            type = "closed"
          )
        )
    }
  }
  
  if (nrow(retention_loops) > 0) {
    
    p <- p +
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
  
  p +
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
      segment.colour = "grey40",
      seed = 1,
      max.overlaps = Inf
    ) +
    
    geom_text(
      data = taupo_point,
      aes(
        x = plot_lon,
        y = plot_lat,
        label = node_id
      ),
      nudge_x = 0.12,
      nudge_y = -0.12,
      size = 3
    ) +
    
    facet_wrap(
      ~source_depth_stratum,
      nrow = 1,
      drop = FALSE
    ) +
    
    scale_colour_manual(
      values = depth_colours,
      breaks = depth_levels,
      limits = depth_levels,
      name = "Destination depth",
      drop = FALSE
    ) +
    
    scale_linewidth_continuous(
      range = c(
        0.35,
        2.8
      ),
      limits = c(
        0,
        max_pooled_connectivity
      ),
      labels = label_percent(
        accuracy = 0.1
      ),
      name = "Connectivity"
    ) +
    
    guides(
      linewidth = guide_legend(
        order = 1,
        title.position = "left"
      ),
      colour = guide_legend(
        order = 2,
        title.position = "left",
        override.aes = list(
          linewidth = 1.8,
          alpha = 1
        )
      )
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
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      ),
      panel.spacing = grid::unit(
        0.65,
        "lines"
      )
    )
}

figure_A_30 <- plot_geographic_connectivity(30)
figure_A_60 <- plot_geographic_connectivity(60)
figure_A_90 <- plot_geographic_connectivity(90)


# Add PLD labels to each row
figure_A_30_clean <- figure_A_30 +
  labs(
    title = "30-day PLD",
    x = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.margin = margin(2, 2, 2, 2)
  )

figure_A_60_clean <- figure_A_60 +
  labs(
    title = "60-day PLD",
    x = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.margin = margin(2, 2, 2, 2)
  )

figure_A_90_clean <- figure_A_90 +
  labs(
    title = "90-day PLD"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.margin = margin(2, 2, 2, 2)
  )


# Combine into one figure
figure_A_all <- (
  figure_A_30_clean /
    figure_A_60_clean /
    figure_A_90_clean
) +
  plot_layout(
    guides = "collect",
    axes = "collect",
    axis_titles = "collect",
    heights = c(1, 1, 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(2, 2, 2, 2),
    legend.box.margin = margin(0, 0, 0, 0)
  )

figure_A_all

ggsave(
  "figures/Figure_2_Taupo_connectivity_all_PLDs.png",
  figure_A_all,
  width = 14,
  height = 10.5,
  units = "in",
  dpi = 400,
  bg = "white"
)


# Figure B
# Pooled multilayer network including retention

layer_positions <- tibble(
  depth_stratum = factor(
    depth_levels,
    levels = depth_levels
  ),
  y = c(
    4,
    3,
    2,
    1
  )
)

make_network_self_loops <- function(edges) {
  
  if (nrow(edges) == 0) {
    return(
      tibble(
        loop_id = character(),
        connectivity = numeric(),
        x = numeric(),
        y = numeric()
      )
    )
  }
  
  bind_rows(
    lapply(
      seq_len(
        nrow(edges)
      ),
      function(i) {
        
        edge <- edges[i, ]
        
        theta <- seq(
          pi / 4,
          2 * pi + pi / 4 - 0.30,
          length.out = 70
        )
        
        tibble(
          loop_id = paste(
            edge$pld_days,
            edge$source_depth_stratum,
            sep = "_"
          ),
          connectivity = edge$connectivity,
          x = 0.88 +
            0.14 * cos(theta),
          y = edge$source_y +
            0.16 * sin(theta)
        )
      }
    )
  )
}


plot_multilayer_network <- function(pld_current) {
  
  edges_current <- pooled_network_edges %>%
    filter(
      pld_days == pld_current
    )
  
  external_current <- edges_current %>%
    filter(
      destination_node != "T22"
    )
  
  retention_current <- edges_current %>%
    filter(
      destination_node == "T22"
    ) %>%
    left_join(
      layer_positions %>%
        transmute(
          source_depth_stratum = depth_stratum,
          source_y = y
        ),
      by = "source_depth_stratum"
    ) %>%
    left_join(
      layer_positions %>%
        transmute(
          destination_depth_stratum = depth_stratum,
          destination_y = y
        ),
      by = "destination_depth_stratum"
    )
  
  receiver_ids <- unique(
    external_current$destination_node
  )
  
  tasmantid_ids <- plot_nodes %>%
    filter(
      node_id %in% receiver_ids,
      grepl("^T", node_id)
    ) %>%
    arrange(
      desc(lat)
    ) %>%
    pull(
      node_id
    )
  
  lord_howe_ids <- plot_nodes %>%
    filter(
      node_id %in% receiver_ids,
      grepl("^L", node_id)
    ) %>%
    arrange(
      desc(lat)
    ) %>%
    pull(
      node_id
    )
  
  tasmantid_positions <- tibble(
    node_id = tasmantid_ids,
    x = seq_along(
      tasmantid_ids
    ) + 1
  )
  
  lord_howe_start <- if (
    nrow(tasmantid_positions) > 0
  ) {
    max(
      tasmantid_positions$x
    ) + 2
  } else {
    3
  }
  
  lord_howe_positions <- tibble(
    node_id = lord_howe_ids,
    x = seq_along(
      lord_howe_ids
    ) + lord_howe_start - 1
  )
  
  seamount_positions <- bind_rows(
    tibble(
      node_id = "T22",
      x = 1
    ),
    tasmantid_positions,
    lord_howe_positions
  )
  
  network_edges <- external_current %>%
    left_join(
      seamount_positions %>%
        transmute(
          source_node = node_id,
          source_x = x
        ),
      by = "source_node"
    ) %>%
    left_join(
      seamount_positions %>%
        transmute(
          destination_node = node_id,
          destination_x = x
        ),
      by = "destination_node"
    ) %>%
    left_join(
      layer_positions %>%
        transmute(
          source_depth_stratum = depth_stratum,
          source_y = y
        ),
      by = "source_depth_stratum"
    ) %>%
    left_join(
      layer_positions %>%
        transmute(
          destination_depth_stratum = depth_stratum,
          destination_y = y
        ),
      by = "destination_depth_stratum"
    )
  
  within_depth_edges <- network_edges %>%
    filter(
      source_depth_stratum ==
        destination_depth_stratum
    )
  
  between_depth_edges <- network_edges %>%
    filter(
      source_depth_stratum !=
        destination_depth_stratum
    )
  
  retention_same <- retention_current %>%
    filter(
      source_depth_stratum ==
        destination_depth_stratum
    )
  
  retention_up <- retention_current %>%
    filter(
      source_depth_stratum !=
        destination_depth_stratum,
      destination_y > source_y
    )
  
  retention_down <- retention_current %>%
    filter(
      source_depth_stratum !=
        destination_depth_stratum,
      destination_y < source_y
    )
  
  retention_loops <- make_network_self_loops(
    retention_same
  )
  
  source_nodes <- layer_positions %>%
    transmute(
      node_id = "T22",
      x = 1,
      depth_stratum,
      y
    )
  
  destination_nodes <- network_edges %>%
    distinct(
      node_id = destination_node,
      x = destination_x,
      depth_stratum = destination_depth_stratum,
      y = destination_y
    )
  
  chain_labels <- tibble(
    label = character(),
    x = numeric(),
    y = numeric()
  )
  
  if (nrow(tasmantid_positions) > 0) {
    
    chain_labels <- bind_rows(
      chain_labels,
      tibble(
        label = "Tasmantid",
        x = mean(
          tasmantid_positions$x
        ),
        y = 4.5
      )
    )
  }
  
  if (nrow(lord_howe_positions) > 0) {
    
    chain_labels <- bind_rows(
      chain_labels,
      tibble(
        label = "Lord Howe",
        x = mean(
          lord_howe_positions$x
        ),
        y = 4.5
      )
    )
  }
  
  ggplot() +
    
    geom_hline(
      yintercept = c(
        1,
        2,
        3,
        4
      ),
      colour = "grey90",
      linewidth = 0.3
    ) +
    
    geom_segment(
      data = between_depth_edges,
      aes(
        x = source_x,
        y = source_y,
        xend = destination_x,
        yend = destination_y,
        linewidth = connectivity
      ),
      colour = "black",
      alpha = 0.82,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.09,
          "cm"
        ),
        type = "closed"
      )
    ) +
    
    geom_curve(
      data = within_depth_edges,
      aes(
        x = source_x,
        y = source_y,
        xend = destination_x,
        yend = destination_y,
        linewidth = connectivity
      ),
      curvature = -0.12,
      colour = "black",
      alpha = 0.82,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.09,
          "cm"
        ),
        type = "closed"
      )
    ) +
    
    geom_curve(
      data = retention_up,
      aes(
        x = 1,
        y = source_y,
        xend = 1,
        yend = destination_y,
        linewidth = connectivity
      ),
      curvature = 0.35,
      colour = "black",
      alpha = 0.90,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.09,
          "cm"
        ),
        type = "closed"
      )
    ) +
    
    geom_curve(
      data = retention_down,
      aes(
        x = 1,
        y = source_y,
        xend = 1,
        yend = destination_y,
        linewidth = connectivity
      ),
      curvature = -0.35,
      colour = "black",
      alpha = 0.90,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.09,
          "cm"
        ),
        type = "closed"
      )
    ) +
    
    geom_path(
      data = retention_loops,
      aes(
        x = x,
        y = y,
        group = loop_id,
        linewidth = connectivity
      ),
      colour = "black",
      alpha = 0.90,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.08,
          "cm"
        ),
        type = "closed"
      )
    ) +
    
    geom_point(
      data = destination_nodes,
      aes(
        x = x,
        y = y
      ),
      shape = 21,
      size = 3.1,
      stroke = 0.8,
      fill = "white",
      colour = "black"
    ) +
    
    geom_point(
      data = source_nodes,
      aes(
        x = x,
        y = y
      ),
      shape = 21,
      size = 3.4,
      stroke = 0.8,
      fill = "black",
      colour = "black"
    ) +
    
    geom_text(
      data = chain_labels,
      aes(
        x = x,
        y = y,
        label = label
      ),
      fontface = "bold",
      size = 3.5
    ) +
    
    scale_x_continuous(
      breaks = seamount_positions$x,
      labels = seamount_positions$node_id,
      expand = expansion(
        add = c(
          0.55,
          0.30
        )
      )
    ) +
    
    scale_y_continuous(
      breaks = c(
        4,
        3,
        2,
        1
      ),
      labels = depth_levels,
      limits = c(
        0.7,
        4.7
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    scale_linewidth_continuous(
      range = c(
        0.35,
        2.8
      ),
      limits = c(
        0,
        max_pooled_connectivity
      ),
      labels = label_percent(
        accuracy = 0.1
      ),
      name = "Connectivity"
    ) +
    
    guides(
      linewidth = guide_legend(
        title.position = "left"
      )
    ) +
    
    labs(
      title = paste0(
        pld_current,
        " d PLD"
      ),
      x = "Seamount",
      y = "Depth stratum"
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold"
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      )
    )
}


multilayer_30 <- plot_multilayer_network(30)
multilayer_60 <- plot_multilayer_network(60)
multilayer_90 <- plot_multilayer_network(90)

figure_B <- (
  multilayer_30 |
    multilayer_60 |
    multilayer_90
) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )


# Figure C
# Monthly multilayer networks
# Rows = PLD, columns = release month
# Includes external connectivity and depth-resolved retention

monthly_network_edges <- connectivity %>%
  filter(
    period_type == "release",
    source_node == "T22",
    connectivity > 0
  ) %>%
  mutate(
    month = case_when(
      analysis_period == "release_1" ~ "June",
      analysis_period == "release_2" ~ "July",
      analysis_period == "release_3" ~ "August"
    ),
    month = factor(
      month,
      levels = c(
        "June",
        "July",
        "August"
      )
    ),
    pld_label = factor(
      paste0(
        pld_days,
        " d PLD"
      ),
      levels = c(
        "30 d PLD",
        "60 d PLD",
        "90 d PLD"
      )
    ),
    source_depth_stratum = factor(
      source_depth_stratum,
      levels = depth_levels
    ),
    destination_depth_stratum = factor(
      destination_depth_stratum,
      levels = depth_levels
    )
  )

monthly_external <- monthly_network_edges %>%
  filter(
    destination_node != "T22"
  )

monthly_retention <- monthly_network_edges %>%
  filter(
    destination_node == "T22"
  )

monthly_receiver_ids <- unique(
  monthly_external$destination_node
)

monthly_tasmantid_ids <- plot_nodes %>%
  filter(
    node_id %in% monthly_receiver_ids,
    grepl("^T", node_id)
  ) %>%
  arrange(
    desc(lat)
  ) %>%
  pull(
    node_id
  )

monthly_lord_howe_ids <- plot_nodes %>%
  filter(
    node_id %in% monthly_receiver_ids,
    grepl("^L", node_id)
  ) %>%
  arrange(
    desc(lat)
  ) %>%
  pull(
    node_id
  )

monthly_tasmantid_positions <- tibble(
  node_id = monthly_tasmantid_ids,
  x = seq_along(
    monthly_tasmantid_ids
  ) + 1
)

monthly_lord_howe_start <- if (
  nrow(monthly_tasmantid_positions) > 0
) {
  max(
    monthly_tasmantid_positions$x
  ) + 2
} else {
  3
}

monthly_lord_howe_positions <- tibble(
  node_id = monthly_lord_howe_ids,
  x = seq_along(
    monthly_lord_howe_ids
  ) + monthly_lord_howe_start - 1
)

monthly_seamount_positions <- bind_rows(
  tibble(
    node_id = "T22",
    x = 1
  ),
  monthly_tasmantid_positions,
  monthly_lord_howe_positions
)

monthly_external_plot <- monthly_external %>%
  left_join(
    monthly_seamount_positions %>%
      transmute(
        source_node = node_id,
        source_x = x
      ),
    by = "source_node"
  ) %>%
  left_join(
    monthly_seamount_positions %>%
      transmute(
        destination_node = node_id,
        destination_x = x
      ),
    by = "destination_node"
  ) %>%
  left_join(
    layer_positions %>%
      transmute(
        source_depth_stratum = depth_stratum,
        source_y = y
      ),
    by = "source_depth_stratum"
  ) %>%
  left_join(
    layer_positions %>%
      transmute(
        destination_depth_stratum = depth_stratum,
        destination_y = y
      ),
    by = "destination_depth_stratum"
  )

monthly_within_depth <- monthly_external_plot %>%
  filter(
    source_depth_stratum ==
      destination_depth_stratum
  ) %>%
  mutate(
    destination_x_plot = destination_x - 0.10
  )

monthly_between_depth <- monthly_external_plot %>%
  filter(
    source_depth_stratum !=
      destination_depth_stratum
  ) %>%
  mutate(
    destination_x_plot = destination_x - 0.10
  )

monthly_retention_plot <- monthly_retention %>%
  left_join(
    layer_positions %>%
      transmute(
        source_depth_stratum = depth_stratum,
        source_y = y
      ),
    by = "source_depth_stratum"
  ) %>%
  left_join(
    layer_positions %>%
      transmute(
        destination_depth_stratum = depth_stratum,
        destination_y = y
      ),
    by = "destination_depth_stratum"
  )

monthly_retention_same <- monthly_retention_plot %>%
  filter(
    source_depth_stratum ==
      destination_depth_stratum
  )

monthly_retention_up <- monthly_retention_plot %>%
  filter(
    source_depth_stratum !=
      destination_depth_stratum,
    destination_y > source_y
  )

monthly_retention_down <- monthly_retention_plot %>%
  filter(
    source_depth_stratum !=
      destination_depth_stratum,
    destination_y < source_y
  )

make_monthly_self_loops <- function(edges) {
  
  if (nrow(edges) == 0) {
    return(
      tibble(
        loop_id = character(),
        month = factor(
          character(),
          levels = c(
            "June",
            "July",
            "August"
          )
        ),
        pld_label = factor(
          character(),
          levels = c(
            "30 d PLD",
            "60 d PLD",
            "90 d PLD"
          )
        ),
        connectivity = numeric(),
        x = numeric(),
        y = numeric()
      )
    )
  }
  
  bind_rows(
    lapply(
      seq_len(
        nrow(edges)
      ),
      function(i) {
        
        edge <- edges[i, ]
        
        theta <- seq(
          pi / 4,
          2 * pi + pi / 4 - 0.30,
          length.out = 60
        )
        
        tibble(
          loop_id = paste(
            edge$analysis_period,
            edge$pld_days,
            edge$source_depth_stratum,
            sep = "_"
          ),
          month = edge$month,
          pld_label = edge$pld_label,
          connectivity = edge$connectivity,
          x = 0.88 +
            0.14 * cos(theta),
          y = edge$source_y +
            0.16 * sin(theta)
        )
      }
    )
  )
}

monthly_self_loops <- make_monthly_self_loops(
  monthly_retention_same
)

panel_grid <- expand_grid(
  pld_label = factor(
    c(
      "30 d PLD",
      "60 d PLD",
      "90 d PLD"
    ),
    levels = c(
      "30 d PLD",
      "60 d PLD",
      "90 d PLD"
    )
  ),
  month = factor(
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

monthly_receiver_depth_nodes <- monthly_external_plot %>%
  distinct(
    node_id = destination_node,
    x = destination_x,
    depth_stratum = destination_depth_stratum,
    y = destination_y
  )

monthly_destination_nodes <- tidyr::crossing(
  panel_grid,
  monthly_receiver_depth_nodes
)

monthly_source_nodes <- tidyr::crossing(
  panel_grid,
  layer_positions %>%
    transmute(
      node_id = "T22",
      x = 1,
      depth_stratum,
      y
    )
)

max_monthly_connectivity <- max(
  monthly_network_edges$connectivity,
  na.rm = TRUE
)

figure_C <- ggplot() +
  
  geom_hline(
    yintercept = c(
      1,
      2,
      3,
      4
    ),
    colour = "grey90",
    linewidth = 0.25
  ) +
  
  geom_segment(
    data = monthly_between_depth,
    aes(
      x = source_x,
      y = source_y,
      xend = destination_x_plot,
      yend = destination_y,
      linewidth = connectivity
    ),
    colour = "black",
    alpha = 0.80,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(
        0.13,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_curve(
    data = monthly_within_depth,
    aes(
      x = source_x,
      y = source_y,
      xend = destination_x_plot,
      yend = destination_y,
      linewidth = connectivity
    ),
    curvature = -0.12,
    colour = "black",
    alpha = 0.80,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(
        0.13,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_curve(
    data = monthly_retention_up,
    aes(
      x = 1,
      y = source_y,
      xend = 0.91,
      yend = destination_y,
      linewidth = connectivity
    ),
    curvature = 0.35,
    colour = "black",
    alpha = 0.88,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(
        0.13,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_curve(
    data = monthly_retention_down,
    aes(
      x = 1,
      y = source_y,
      xend = 1.09,
      yend = destination_y,
      linewidth = connectivity
    ),
    curvature = -0.35,
    colour = "black",
    alpha = 0.88,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(
        0.13,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_path(
    data = monthly_self_loops,
    aes(
      x = x,
      y = y,
      group = loop_id,
      linewidth = connectivity
    ),
    colour = "black",
    alpha = 0.88,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(
        0.11,
        "cm"
      ),
      type = "closed"
    )
  ) +
  
  geom_point(
    data = monthly_destination_nodes,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 2.1,
    stroke = 0.55,
    fill = "white",
    colour = "black"
  ) +
  
  geom_point(
    data = monthly_source_nodes,
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 2.4,
    stroke = 0.55,
    fill = "black",
    colour = "black"
  ) +
  
  facet_grid(
    pld_label ~ month
  ) +
  
  scale_x_continuous(
    breaks = monthly_seamount_positions$x,
    labels = monthly_seamount_positions$node_id,
    expand = expansion(
      add = c(
        0.55,
        0.30
      )
    )
  ) +
  
  scale_y_continuous(
    breaks = c(
      4,
      3,
      2,
      1
    ),
    labels = depth_levels,
    limits = c(
      0.65,
      4.65
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  scale_linewidth_continuous(
    range = c(
      0.25,
      2.0
    ),
    limits = c(
      0,
      max_monthly_connectivity
    ),
    labels = label_percent(
      accuracy = 0.1
    ),
    name = "Connectivity"
  ) +
  
  guides(
    linewidth = guide_legend(
      title.position = "left"
    )
  ) +
  
  labs(
    x = "Seamount",
    y = "Depth stratum"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.5
    ),
    axis.text.y = element_text(
      size = 8
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    panel.spacing = grid::unit(
      0.35,
      "lines"
    )
  )

# Display figures

figure_A_30
figure_A_60
figure_A_90
figure_B
figure_C


# Save figures

ggsave(
  "figures/Figure_A_30d_geographic_connectivity.png",
  figure_A_30,
  width = 12,
  height = 4.8,
  dpi = 400
)

ggsave(
  "figures/Figure_A_60d_geographic_connectivity.png",
  figure_A_60,
  width = 12,
  height = 4.8,
  dpi = 400
)

ggsave(
  "figures/Figure_A_90d_geographic_connectivity.png",
  figure_A_90,
  width = 12,
  height = 4.8,
  dpi = 400
)

ggsave(
  "figures/Figure_B_multilayer_network.png",
  figure_B,
  width = 14,
  height = 5.5,
  dpi = 400
)

ggsave(
  "figures/Figure_C_monthly_multilayer_networks.png",
  figure_C,
  width = 15,
  height = 10,
  dpi = 400
)

