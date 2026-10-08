# Taupo Bank Analysis pt. 6 :) 
# Pretty figures !!

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(scales)

source("00_setup_and_inputs_.R")

dir.create("figures", showWarnings = FALSE)
dir.create("html", showWarnings = FALSE)

analysis <- readRDS("connectivity_analysis.rds")
settlement_events <- readRDS("settlement_events.rds")
network_nodes <- readRDS("network_nodes.rds")
particle_trajectories <- readRDS("particle_trajectories.rds")

connectivity <- analysis$connectivity
network_roles <- analysis$network_roles
depth_connectivity <- analysis$depth_connectivity
within_between_depth <- analysis$within_between_depth
chain_connectivity <- analysis$chain_connectivity
collapsed_connectivity <- analysis$collapsed_connectivity
depth_levels <- analysis$depth_levels

pld_values <- c(30, 60, 90)

scenario_lookup <- settlement_events %>%
  distinct(pld_days, precompetency_days) %>%
  arrange(pld_days)


# Colours

pld_colours <- c(
  "30" = "#20A486",
  "60" = "#7B2CBF",
  "90" = "#E76F91"
)

depth_layer_colours <- c(
  "200-700 m" = "#BFD7EA",
  "700-1500 m" = "#8FB9D4",
  "1500-2500 m" = "#5F9BC6",
  ">2500 m" = "#2F6FA3"
)

transition_colours <- c(
  "within-depth" = "#355C7D",
  "between-depth" = "#E0A458"
)

connectivity_colours <- c(
  "#F7F4F9",
  "#E7E1EF",
  "#C994C7",
  "#88419D",
  "#4A1486"
)

chain_names <- unique(candidate_nodes$chain)

chain_palette <- c(
  "#2F6B5F",
  "#C47C2C",
  "#6A4C93",
  "#4D908E"
)

chain_colours <- setNames(
  chain_palette[seq_along(chain_names)],
  chain_names
)


# Static figure theme

theme_project <- theme_classic(
  base_size = 12
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      face = "bold"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )


# Pooled plotting tables

pooled_connectivity <- connectivity %>%
  filter(
    analysis_period == "pooled"
  )

pooled_network_roles <- network_roles %>%
  filter(
    analysis_period == "pooled"
  )

pooled_depth_connectivity <- depth_connectivity %>%
  filter(
    analysis_period == "pooled"
  )

pooled_within_between_depth <- within_between_depth %>%
  filter(
    analysis_period == "pooled"
  )

pooled_chain_connectivity <- chain_connectivity %>%
  filter(
    analysis_period == "pooled"
  )

pooled_collapsed_connectivity <- collapsed_connectivity %>%
  filter(
    analysis_period == "pooled"
  )


# Figure 1
# Within-depth versus between-depth settlement

figure_within_between <- pooled_within_between_depth %>%
  mutate(
    settlement_percent =
      100 * settlement_fraction,
    
    pld_label =
      factor(
        paste0(
          pld_days,
          " d"
        ),
        levels = c(
          "30 d",
          "60 d",
          "90 d"
        )
      )
  ) %>%
  ggplot(
    aes(
      pld_label,
      settlement_percent,
      fill = depth_transition
    )
  ) +
  geom_col(
    position =
      position_dodge(
        width = 0.75
      ),
    width = 0.65
  ) +
  scale_fill_manual(
    values = transition_colours,
    labels = c(
      "within-depth" = "Within depth",
      "between-depth" = "Between depth"
    )
  ) +
  labs(
    title =
      "Within- and between-depth settlement",
    
    subtitle =
      "Pooled across release periods",
    
    x =
      "Pelagic larval duration",
    
    y =
      "Settlement probability (%)",
    
    fill =
      NULL
  ) +
  theme_project

ggsave(
  "figures/01_within_between_depth.png",
  figure_within_between,
  width = 7,
  height = 4.5,
  dpi = 400
)


# Figure 2
# Connectivity among depth strata

figure_depth_heatmap <- pooled_depth_connectivity %>%
  mutate(
    source_depth_stratum =
      factor(
        source_depth_stratum,
        levels =
          rev(
            depth_levels
          )
      ),
    
    destination_depth_stratum =
      factor(
        destination_depth_stratum,
        levels =
          depth_levels
      ),
    
    pld_label =
      factor(
        paste0(
          pld_days,
          " d PLD"
        ),
        levels = c(
          "30 d PLD",
          "60 d PLD",
          "90 d PLD"
        )
      )
  ) %>%
  ggplot(
    aes(
      destination_depth_stratum,
      source_depth_stratum,
      fill = depth_connectivity
    )
  ) +
  geom_tile(
    colour = "white",
    linewidth = 1
  ) +
  geom_text(
    aes(
      label =
        percent(
          depth_connectivity,
          accuracy = 0.1
        )
    ),
    size = 3
  ) +
  scale_fill_gradientn(
    colours =
      connectivity_colours,
    
    labels =
      label_percent(
        accuracy = 0.1
      )
  ) +
  facet_wrap(
    ~pld_label,
    nrow = 1
  ) +
  coord_equal() +
  labs(
    title =
      "Connectivity among depth strata",
    
    subtitle =
      "Rows are source depths; columns are settlement depths",
    
    x =
      "Destination depth",
    
    y =
      "Source depth",
    
    fill =
      "Connectivity"
  ) +
  theme_project +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

ggsave(
  "figures/02_depth_transition_heatmap.png",
  figure_depth_heatmap,
  width = 11,
  height = 4.5,
  dpi = 400
)


# Figure 3
# Within- and between-chain connectivity

figure_chain <- pooled_chain_connectivity %>%
  mutate(
    connectivity_percent =
      100 * chain_connectivity,
    
    pld_label =
      factor(
        paste0(
          pld_days,
          " d"
        ),
        levels = c(
          "30 d",
          "60 d",
          "90 d"
        )
      )
  ) %>%
  ggplot(
    aes(
      pld_label,
      connectivity_percent,
      fill = destination_chain
    )
  ) +
  geom_col(
    position =
      position_dodge(
        width = 0.75
      ),
    width = 0.65
  ) +
  facet_wrap(
    ~source_chain
  ) +
  scale_fill_manual(
    values =
      chain_colours
  ) +
  labs(
    title =
      "Connectivity within and between seamount chains",
    
    subtitle =
      "Pooled across release periods",
    
    x =
      "Pelagic larval duration",
    
    y =
      "Connectivity probability (%)",
    
    fill =
      "Destination chain"
  ) +
  theme_project

ggsave(
  "figures/03_chain_connectivity.png",
  figure_chain,
  width = 7.5,
  height = 4.5,
  dpi = 400
)


# Figure 4
# External export and retention by depth

figure_roles <- pooled_network_roles %>%
  select(
    node_id,
    pld_days,
    depth_stratum,
    s_out,
    retention
  ) %>%
  pivot_longer(
    c(
      s_out,
      retention
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value_percent =
      100 * value,
    
    depth_stratum =
      factor(
        depth_stratum,
        levels =
          depth_levels
      ),
    
    pld_label =
      factor(
        paste0(
          pld_days,
          " d"
        ),
        levels = c(
          "30 d",
          "60 d",
          "90 d"
        )
      ),
    
    metric =
      recode(
        metric,
        s_out =
          "External out-strength",
        retention =
          "Retention"
      )
  ) %>%
  ggplot(
    aes(
      depth_stratum,
      value_percent,
      group =
        interaction(
          node_id,
          pld_label
        ),
      colour =
        pld_label
    )
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~metric,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "30 d" =
        pld_colours[["30"]],
      
      "60 d" =
        pld_colours[["60"]],
      
      "90 d" =
        pld_colours[["90"]]
    )
  ) +
  labs(
    title =
      "Connectivity roles across source depths",
    
    subtitle =
      "Pooled Taupo Bank pilot",
    
    x =
      "Source depth stratum",
    
    y =
      "Probability / strength (%)",
    
    colour =
      "PLD"
  ) +
  theme_project +
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      )
  )

ggsave(
  "figures/04_export_retention_by_depth.png",
  figure_roles,
  width = 9,
  height = 4.5,
  dpi = 400
)


# Figure 5
# Geographic depth-collapsed connectivity

map_edges <- pooled_collapsed_connectivity %>%
  filter(
    connectivity > 0,
    source_node != destination_node
  ) %>%
  left_join(
    candidate_nodes %>%
      transmute(
        source_node =
          node_id,
        source_lon =
          lon,
        source_lat =
          lat
      ),
    by =
      "source_node"
  ) %>%
  left_join(
    candidate_nodes %>%
      transmute(
        destination_node =
          node_id,
        destination_lon =
          lon,
        destination_lat =
          lat
      ),
    by =
      "destination_node"
  ) %>%
  mutate(
    pld_label =
      factor(
        paste0(
          pld_days,
          " d PLD"
        ),
        levels = c(
          "30 d PLD",
          "60 d PLD",
          "90 d PLD"
        )
      )
  )

figure_connectivity_map <- ggplot() +
  geom_segment(
    data =
      map_edges,
    
    aes(
      x =
        source_lon,
      
      y =
        source_lat,
      
      xend =
        destination_lon,
      
      yend =
        destination_lat,
      
      linewidth =
        connectivity
    ),
    
    colour =
      "grey25",
    
    alpha =
      0.75,
    
    arrow =
      grid::arrow(
        length =
          grid::unit(
            0.12,
            "cm"
          ),
        type =
          "closed"
      )
  ) +
  geom_point(
    data =
      candidate_nodes,
    
    aes(
      lon,
      lat
    ),
    
    shape =
      21,
    
    size =
      2.5,
    
    stroke =
      0.6,
    
    fill =
      "white",
    
    colour =
      "black"
  ) +
  geom_point(
    data =
      candidate_nodes %>%
      filter(
        node_id == "T22"
      ),
    
    aes(
      lon,
      lat
    ),
    
    shape =
      21,
    
    size =
      3.3,
    
    fill =
      "black",
    
    colour =
      "black"
  ) +
  geom_text(
    data =
      candidate_nodes,
    
    aes(
      lon,
      lat,
      label = node_id
    ),
    
    size =
      2.5,
    
    nudge_x =
      0.08,
    
    check_overlap =
      TRUE
  ) +
  facet_wrap(
    ~pld_label,
    nrow = 1
  ) +
  scale_linewidth_continuous(
    range = c(
      0.3,
      2.2
    ),
    labels =
      label_percent(
        accuracy = 0.1
      )
  ) +
  coord_quickmap() +
  labs(
    title =
      "Geographic connectivity from Taupo Bank",
    
    subtitle =
      "Depth-collapsed connectivity; line width represents probability",
    
    x =
      "Longitude",
    
    y =
      "Latitude",
    
    linewidth =
      "Connectivity"
  ) +
  theme_project

ggsave(
  "figures/05_geographic_connectivity.png",
  figure_connectivity_map,
  width = 11,
  height = 6,
  dpi = 400
)


# Figure 6
# Observed depth-resolved connectivity

connected_destinations <- pooled_connectivity %>%
  filter(
    connectivity > 0
  ) %>%
  distinct(
    destination_node,
    destination_depth_stratum
  ) %>%
  mutate(
    destination_label =
      paste(
        destination_node,
        destination_depth_stratum,
        sep = " | "
      ),
    
    depth_order =
      match(
        destination_depth_stratum,
        depth_levels
      )
  ) %>%
  arrange(
    destination_node,
    depth_order
  ) %>%
  pull(
    destination_label
  )

figure_connectivity_heatmap <- pooled_connectivity %>%
  filter(
    connectivity > 0
  ) %>%
  mutate(
    destination_label =
      paste(
        destination_node,
        destination_depth_stratum,
        sep = " | "
      ),
    
    destination_label =
      factor(
        destination_label,
        levels =
          connected_destinations
      ),
    
    source_depth_stratum =
      factor(
        source_depth_stratum,
        levels =
          rev(
            depth_levels
          )
      ),
    
    pld_label =
      factor(
        paste0(
          pld_days,
          " d PLD"
        ),
        levels = c(
          "30 d PLD",
          "60 d PLD",
          "90 d PLD"
        )
      )
  ) %>%
  ggplot(
    aes(
      destination_label,
      source_depth_stratum,
      fill = connectivity
    )
  ) +
  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +
  scale_fill_gradientn(
    colours =
      connectivity_colours,
    
    labels =
      label_percent(
        accuracy = 0.1
      )
  ) +
  facet_wrap(
    ~pld_label,
    ncol = 1
  ) +
  labs(
    title =
      "Observed depth-resolved connectivity",
    
    subtitle =
      "Non-zero pooled connections from the Taupo Bank pilot",
    
    x =
      "Destination seamount and depth",
    
    y =
      "Source depth",
    
    fill =
      "Connectivity"
  ) +
  theme_project +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        size = 8
      )
  )

ggsave(
  "figures/06_depth_resolved_connectivity_heatmap.png",
  figure_connectivity_heatmap,
  width = 12,
  height = 8,
  dpi = 400
)


# 3-D helpers

make_bathymetry_grid <- function(
    gebco,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    aggregate_factor = 8) {
  
  bathy <- terra::crop(
    gebco,
    terra::ext(
      lon_min,
      lon_max,
      lat_min,
      lat_max
    )
  )
  
  bathy <- terra::aggregate(
    bathy,
    fact =
      aggregate_factor,
    fun =
      mean,
    na.rm =
      TRUE
  )
  
  list(
    x =
      terra::xFromCol(
        bathy,
        seq_len(
          ncol(
            bathy
          )
        )
      ),
    
    y =
      terra::yFromRow(
        bathy,
        seq_len(
          nrow(
            bathy
          )
        )
      ),
    
    z =
      as.matrix(
        bathy,
        wide = TRUE
      )
  )
}


make_plot_extent <- function(
    lon,
    lat,
    buffer_deg = 0.25) {
  
  c(
    lon_min =
      max(
        min(
          lon,
          na.rm = TRUE
        ) -
          buffer_deg,
        terra::xmin(
          gebco
        )
      ),
    
    lon_max =
      min(
        max(
          lon,
          na.rm = TRUE
        ) +
          buffer_deg,
        terra::xmax(
          gebco
        )
      ),
    
    lat_min =
      max(
        min(
          lat,
          na.rm = TRUE
        ) -
          buffer_deg,
        terra::ymin(
          gebco
        )
      ),
    
    lat_max =
      min(
        max(
          lat,
          na.rm = TRUE
        ) +
          buffer_deg,
        terra::ymax(
          gebco
        )
      )
  )
}


add_bathymetry <- function(
    figure,
    bathymetry_grid) {
  
  figure %>%
    plotly::add_surface(
      x =
        bathymetry_grid$x,
      
      y =
        bathymetry_grid$y,
      
      z =
        bathymetry_grid$z,
      
      colorscale =
        "Viridis",
      
      showscale =
        TRUE,
      
      colorbar =
        list(
          title =
            "GEBCO<br>elevation (m)",
          
          len =
            0.55,
          
          y =
            0.25
        ),
      
      hovertemplate =
        paste0(
          "Longitude: %{x:.3f}",
          "<br>Latitude: %{y:.3f}",
          "<br>Elevation: %{z:.0f} m",
          "<extra></extra>"
        ),
      
      showlegend =
        FALSE
    )
}


add_depth_planes <- function(
    figure,
    bathymetry_grid) {
  
  depth_boundaries <- c(
    200,
    700,
    1500,
    2500
  )
  
  plane_colours <- unname(
    depth_layer_colours
  )
  
  x_plane <- range(
    bathymetry_grid$x
  )
  
  y_plane <- range(
    bathymetry_grid$y
  )
  
  for (
    i in seq_along(
      depth_boundaries
    )
  ) {
    
    figure <- figure %>%
      plotly::add_surface(
        x =
          x_plane,
        
        y =
          y_plane,
        
        z =
          matrix(
            -depth_boundaries[i],
            nrow = 2,
            ncol = 2
          ),
        
        colorscale =
          list(
            c(
              0,
              plane_colours[i]
            ),
            c(
              1,
              plane_colours[i]
            )
          ),
        
        showscale =
          FALSE,
        
        opacity =
          0.15,
        
        hoverinfo =
          "skip",
        
        showlegend =
          FALSE
      )
  }
  
  label_x <-
    min(
      bathymetry_grid$x
    ) +
    0.025 *
    diff(
      range(
        bathymetry_grid$x
      )
    )
  
  label_y <-
    max(
      bathymetry_grid$y
    ) -
    0.025 *
    diff(
      range(
        bathymetry_grid$y
      )
    )
  
  figure %>%
    plotly::add_trace(
      x =
        rep(
          label_x,
          4
        ),
      
      y =
        rep(
          label_y,
          4
        ),
      
      z =
        -depth_boundaries,
      
      type =
        "scatter3d",
      
      mode =
        "text",
      
      text =
        paste0(
          depth_boundaries,
          " m"
        ),
      
      textfont =
        list(
          color =
            "#17324D",
          
          size =
            11
        ),
      
      hoverinfo =
        "skip",
      
      showlegend =
        FALSE
    )
}


# GEBCO elevation at node coordinates
#
# terra::extract is explicit because tidyr also has extract().

add_node_elevation <- function(
    nodes) {
  
  node_points <- terra::vect(
    nodes,
    geom = c(
      "lon",
      "lat"
    ),
    crs =
      "EPSG:4326"
  )
  
  node_elevation <- terra::extract(
    gebco,
    node_points
  )
  
  nodes %>%
    mutate(
      surface_elevation_m =
        node_elevation[[2]]
    )
}


# Separate particle lines with an NA row.

make_trajectory_segments <- function(
    trajectory_data) {
  
  trajectory_data <- trajectory_data %>%
    arrange(
      release_id,
      particle_id,
      age_hours
    )
  
  groups <- split(
    trajectory_data,
    interaction(
      trajectory_data$release_id,
      trajectory_data$particle_id,
      drop = TRUE
    )
  )
  
  bind_rows(
    lapply(
      groups,
      function(x) {
        
        separator <-
          x[
            1,
            ,
            drop = FALSE
          ]
        
        separator[,] <- NA
        
        bind_rows(
          x,
          separator
        )
      }
    )
  )
}


# Separate network edges with NA rows.

make_edge_segments <- function(
    edge_data) {
  
  bind_rows(
    lapply(
      seq_len(
        nrow(
          edge_data
        )
      ),
      function(i) {
        
        tibble(
          lon =
            c(
              edge_data$source_lon[i],
              edge_data$destination_lon[i],
              NA_real_
            ),
          
          lat =
            c(
              edge_data$source_lat[i],
              edge_data$destination_lat[i],
              NA_real_
            ),
          
          z =
            c(
              edge_data$source_z[i],
              edge_data$destination_z[i],
              NA_real_
            ),
          
          hover_text =
            rep(
              edge_data$edge_hover[i],
              3
            )
        )
      }
    )
  )
}


# Interactive trajectory HTML

external_settlements <- settlement_events %>%
  filter(
    destination_node !=
      source_node
  ) %>%
  left_join(
    candidate_nodes %>%
      transmute(
        destination_node =
          node_id,
        
        destination_name =
          name
      ),
    by =
      "destination_node"
  )


trajectory_limits <- external_settlements %>%
  group_by(
    release_id,
    particle_id
  ) %>%
  summarise(
    maximum_settlement_age_days =
      max(
        settlement_age_days
      ),
    
    .groups =
      "drop"
  )


trajectory_extent_data <- particle_trajectories %>%
  inner_join(
    trajectory_limits,
    by = c(
      "release_id",
      "particle_id"
    )
  ) %>%
  filter(
    age_hours / 24 <=
      maximum_settlement_age_days
  )


trajectory_extent <- make_plot_extent(
  trajectory_extent_data$lon,
  trajectory_extent_data$lat,
  buffer_deg = 0.25
)


trajectory_bathymetry <- make_bathymetry_grid(
  gebco =
    gebco,
  
  lon_min =
    trajectory_extent[["lon_min"]],
  
  lon_max =
    trajectory_extent[["lon_max"]],
  
  lat_min =
    trajectory_extent[["lat_min"]],
  
  lat_max =
    trajectory_extent[["lat_max"]],
  
  aggregate_factor =
    8
)


trajectory_3d <- plotly::plot_ly()

trajectory_trace_group <- character()


# Bathymetry stays visible.

trajectory_3d <- add_bathymetry(
  trajectory_3d,
  trajectory_bathymetry
)

trajectory_trace_group <- c(
  trajectory_trace_group,
  "base"
)


# Four depth surfaces + one text trace stay visible.

trajectory_3d <- add_depth_planes(
  trajectory_3d,
  trajectory_bathymetry
)

trajectory_trace_group <- c(
  trajectory_trace_group,
  rep(
    "base",
    5
  )
)


# Seamount labels

trajectory_nodes <- candidate_nodes %>%
  filter(
    lon >=
      trajectory_extent[["lon_min"]],
    
    lon <=
      trajectory_extent[["lon_max"]],
    
    lat >=
      trajectory_extent[["lat_min"]],
    
    lat <=
      trajectory_extent[["lat_max"]]
  ) %>%
  add_node_elevation()


trajectory_3d <- trajectory_3d %>%
  plotly::add_trace(
    data =
      trajectory_nodes,
    
    x =
      ~lon,
    
    y =
      ~lat,
    
    z =
      ~surface_elevation_m,
    
    type =
      "scatter3d",
    
    mode =
      "text",
    
    text =
      ~node_id,
    
    textfont =
      list(
        color =
          "#222222",
        
        size =
          10
      ),
    
    hoverinfo =
      "text",
    
    showlegend =
      FALSE
  )

trajectory_trace_group <- c(
  trajectory_trace_group,
  "base"
)


# PLD-specific trajectories

for (
  pld_current in pld_values
) {
  
  competency_current <- scenario_lookup %>%
    filter(
      pld_days ==
        pld_current
    ) %>%
    pull(
      precompetency_days
    )
  
  settlement_current <- external_settlements %>%
    filter(
      pld_days ==
        pld_current
    )
  
  trajectory_current <- particle_trajectories %>%
    inner_join(
      settlement_current %>%
        select(
          release_id,
          particle_id,
          precompetency_days,
          destination_node,
          destination_name,
          destination_depth_stratum,
          settlement_age_days,
          settlement_depth_m
        ),
      by = c(
        "release_id",
        "particle_id"
      )
    ) %>%
    filter(
      age_hours / 24 <=
        settlement_age_days
    ) %>%
    mutate(
      hover_text =
        paste0(
          "<b>Particle ",
          particle_id,
          "</b>",
          
          "<br>PLD: ",
          pld_current,
          " d",
          
          "<br>Competency from: ",
          precompetency_days,
          " d",
          
          "<br>Particle age: ",
          round(
            age_hours / 24,
            2
          ),
          " d",
          
          "<br>Particle depth: ",
          round(
            depth_m
          ),
          " m",
          
          "<br>Release depth: ",
          round(
            source_depth_m
          ),
          " m",
          
          "<br>Release stratum: ",
          source_depth_stratum,
          
          "<br>Receiver: ",
          destination_node,
          " — ",
          destination_name,
          
          "<br>Destination stratum: ",
          destination_depth_stratum
        )
    )
  
  trajectory_segments <- make_trajectory_segments(
    trajectory_current
  )
  
  
  # Trajectory lines
  
  trajectory_3d <- trajectory_3d %>%
    plotly::add_trace(
      data =
        trajectory_segments,
      
      x =
        ~lon,
      
      y =
        ~lat,
      
      z =
        ~-depth_m,
      
      type =
        "scatter3d",
      
      mode =
        "lines",
      
      name =
        paste0(
          "PLD ",
          pld_current,
          " d / competency ",
          competency_current,
          " d"
        ),
      
      legendgroup =
        paste0(
          "PLD_",
          pld_current
        ),
      
      text =
        ~hover_text,
      
      hoverinfo =
        "text",
      
      line =
        list(
          color =
            pld_colours[[
              as.character(
                pld_current
              )
            ]],
          
          width =
            3
        ),
      
      opacity =
        0.75,
      
      visible =
        pld_current ==
        pld_values[1],
      
      showlegend =
        TRUE
    )
  
  trajectory_trace_group <- c(
    trajectory_trace_group,
    as.character(
      pld_current
    )
  )
  
  
  # Release circles
  
  release_points <- trajectory_current %>%
    group_by(
      release_id,
      particle_id
    ) %>%
    slice_min(
      age_hours,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    mutate(
      release_hover =
        paste0(
          "<b>Release</b>",
          
          "<br>Particle: ",
          particle_id,
          
          "<br>PLD: ",
          pld_current,
          " d",
          
          "<br>Release depth: ",
          round(
            source_depth_m
          ),
          " m",
          
          "<br>Release stratum: ",
          source_depth_stratum
        )
    )
  
  trajectory_3d <- trajectory_3d %>%
    plotly::add_trace(
      data =
        release_points,
      
      x =
        ~lon,
      
      y =
        ~lat,
      
      z =
        ~-depth_m,
      
      type =
        "scatter3d",
      
      mode =
        "markers",
      
      legendgroup =
        paste0(
          "PLD_",
          pld_current
        ),
      
      text =
        ~release_hover,
      
      hoverinfo =
        "text",
      
      marker =
        list(
          color =
            pld_colours[[
              as.character(
                pld_current
              )
            ]],
          
          size =
            3.5,
          
          symbol =
            "circle",
          
          line =
            list(
              color =
                "white",
              
              width =
                1
            )
        ),
      
      visible =
        pld_current ==
        pld_values[1],
      
      showlegend =
        FALSE
    )
  
  trajectory_trace_group <- c(
    trajectory_trace_group,
    as.character(
      pld_current
    )
  )
  
  
  # First valid external settlement diamonds
  
  settlement_current <- settlement_current %>%
    mutate(
      settlement_hover =
        paste0(
          "<b>First valid settlement</b>",
          
          "<br>Particle: ",
          particle_id,
          
          "<br>PLD: ",
          pld_days,
          " d",
          
          "<br>Competency from: ",
          precompetency_days,
          " d",
          
          "<br>Settlement age: ",
          round(
            settlement_age_days,
            2
          ),
          " d",
          
          "<br>Particle depth: ",
          round(
            particle_depth_m
          ),
          " m",
          
          "<br>Receiver: ",
          destination_node,
          " — ",
          destination_name,
          
          "<br>Receiver benthic depth: ",
          round(
            settlement_depth_m
          ),
          " m",
          
          "<br>Destination stratum: ",
          destination_depth_stratum
        )
    )
  
  trajectory_3d <- trajectory_3d %>%
    plotly::add_trace(
      data =
        settlement_current,
      
      x =
        ~settlement_lon,
      
      y =
        ~settlement_lat,
      
      z =
        ~-particle_depth_m,
      
      type =
        "scatter3d",
      
      mode =
        "markers",
      
      legendgroup =
        paste0(
          "PLD_",
          pld_current
        ),
      
      text =
        ~settlement_hover,
      
      hoverinfo =
        "text",
      
      marker =
        list(
          color =
            pld_colours[[
              as.character(
                pld_current
              )
            ]],
          
          size =
            5,
          
          symbol =
            "diamond",
          
          line =
            list(
              color =
                "black",
              
              width =
                1
            )
        ),
      
      visible =
        pld_current ==
        pld_values[1],
      
      showlegend =
        FALSE
    )
  
  trajectory_trace_group <- c(
    trajectory_trace_group,
    as.character(
      pld_current
    )
  )
}


# PLD buttons for trajectories

trajectory_buttons <- lapply(
  pld_values,
  function(
    pld_current) {
    
    competency_current <- scenario_lookup %>%
      filter(
        pld_days ==
          pld_current
      ) %>%
      pull(
        precompetency_days
      )
    
    list(
      method =
        "update",
      
      label =
        paste0(
          pld_current,
          " d"
        ),
      
      args =
        list(
          list(
            visible =
              trajectory_trace_group ==
              "base" |
              trajectory_trace_group ==
              as.character(
                pld_current
              )
          ),
          
          list(
            title =
              paste0(
                "Taupo Bank pilot: 3-D external settlement trajectories",
                
                "<br><sup>PLD ",
                pld_current,
                " d; competency from ",
                competency_current,
                " d; colour = PLD; circles = release; diamonds = first valid settlement</sup>"
              )
          )
        )
    )
  }
)


initial_competency <- scenario_lookup %>%
  filter(
    pld_days ==
      pld_values[1]
  ) %>%
  pull(
    precompetency_days
  )


trajectory_3d <- trajectory_3d %>%
  plotly::layout(
    title =
      paste0(
        "Taupo Bank pilot: 3-D external settlement trajectories",
        
        "<br><sup>PLD ",
        pld_values[1],
        " d; competency from ",
        initial_competency,
        " d; colour = PLD; circles = release; diamonds = first valid settlement</sup>"
      ),
    
    updatemenus =
      list(
        list(
          type =
            "buttons",
          
          direction =
            "right",
          
          x =
            0.02,
          
          y =
            1.10,
          
          xanchor =
            "left",
          
          yanchor =
            "top",
          
          buttons =
            trajectory_buttons,
          
          bgcolor =
            "rgba(255,255,255,0.85)",
          
          bordercolor =
            "#666666"
        )
      ),
    
    legend =
      list(
        title =
          list(
            text =
              "Biological scenario"
          ),
        
        groupclick =
          "togglegroup",
        
        yanchor =
          "top",
        
        y =
          0.55
      ),
    
    scene =
      list(
        xaxis =
          list(
            title =
              "Longitude"
          ),
        
        yaxis =
          list(
            title =
              "Latitude"
          ),
        
        zaxis =
          list(
            title =
              "Depth / elevation (m)"
          ),
        
        aspectmode =
          "manual",
        
        aspectratio =
          list(
            x =
              0.7,
            
            y =
              2,
            
            z =
              0.9
          ),
        
        camera =
          list(
            eye =
              list(
                x =
                  1.4,
                
                y =
                  -1.5,
                
                z =
                  1
              )
          )
      ),
    
    margin =
      list(
        b =
          40,
        
        l =
          60,
        
        t =
          80,
        
        r =
          10
      )
  ) %>%
  plotly::config(
    displaylogo =
      FALSE,
    
    responsive =
      TRUE
  )


htmlwidgets::saveWidget(
  trajectory_3d,
  "html/trajectory_3D.html",
  selfcontained = TRUE
)


# Interactive 3-D multilayer connectivity
#
# Representative z positions are used only to display the four
# analytical strata.
#
# Edge colour = PLD
# Edge width = connectivity strength
# Diamond = destination endpoint

depth_node_z <- c(
  "200-700 m" = -450,
  "700-1500 m" = -1100,
  "1500-2500 m" = -2000,
  ">2500 m" = -3000
)


network_extent_nodes <- candidate_nodes %>%
  semi_join(
    network_nodes %>%
      distinct(
        node_id
      ),
    by =
      "node_id"
  )


connectivity_extent <- make_plot_extent(
  network_extent_nodes$lon,
  network_extent_nodes$lat,
  buffer_deg = 0.35
)


connectivity_bathymetry <- make_bathymetry_grid(
  gebco =
    gebco,
  
  lon_min =
    connectivity_extent[["lon_min"]],
  
  lon_max =
    connectivity_extent[["lon_max"]],
  
  lat_min =
    connectivity_extent[["lat_min"]],
  
  lat_max =
    connectivity_extent[["lat_max"]],
  
  aggregate_factor =
    12
)


node_positions <- network_nodes %>%
  left_join(
    candidate_nodes %>%
      select(
        node_id,
        lon,
        lat
      ),
    by =
      "node_id"
  ) %>%
  mutate(
    z_layer =
      unname(
        depth_node_z[
          depth_stratum
        ]
      ),
    
    node_hover =
      paste0(
        "<b>",
        name,
        " (",
        node_id,
        ")</b>",
        
        "<br>Chain: ",
        chain,
        
        "<br>Depth stratum: ",
        depth_stratum
      )
  )


connectivity_edges <- pooled_connectivity %>%
  filter(
    connectivity > 0,
    
    !(
      source_node ==
        destination_node &
        source_depth_stratum ==
        destination_depth_stratum
    )
  ) %>%
  left_join(
    node_positions %>%
      transmute(
        source_node =
          node_id,
        
        source_depth_stratum =
          depth_stratum,
        
        source_lon =
          lon,
        
        source_lat =
          lat,
        
        source_z =
          z_layer
      ),
    
    by = c(
      "source_node",
      "source_depth_stratum"
    )
  ) %>%
  left_join(
    node_positions %>%
      transmute(
        destination_node =
          node_id,
        
        destination_depth_stratum =
          depth_stratum,
        
        destination_lon =
          lon,
        
        destination_lat =
          lat,
        
        destination_z =
          z_layer
      ),
    
    by = c(
      "destination_node",
      "destination_depth_stratum"
    )
  ) %>%
  filter(
    is.finite(
      source_lon
    ),
    
    is.finite(
      source_lat
    ),
    
    is.finite(
      source_z
    ),
    
    is.finite(
      destination_lon
    ),
    
    is.finite(
      destination_lat
    ),
    
    is.finite(
      destination_z
    )
  ) %>%
  group_by(
    pld_days
  ) %>%
  mutate(
    strength_bin =
      ntile(
        connectivity,
        3
      ),
    
    strength_class =
      case_when(
        strength_bin ==
          1 ~
          "Weaker",
        
        strength_bin ==
          2 ~
          "Intermediate",
        
        TRUE ~
          "Stronger"
      ),
    
    edge_hover =
      paste0(
        "<b>",
        source_node,
        " | ",
        source_depth_stratum,
        "</b> → <b>",
        destination_node,
        " | ",
        destination_depth_stratum,
        "</b>",
        
        "<br>PLD: ",
        pld_days,
        " d",
        
        "<br>Connectivity: ",
        percent(
          connectivity,
          accuracy = 0.01
        ),
        
        "<br>Settled: ",
        N_settled,
        " / ",
        N_released
      )
  ) %>%
  ungroup()


edge_widths <- c(
  "Weaker" = 1.5,
  "Intermediate" = 3.5,
  "Stronger" = 6
)


connectivity_3d <- plotly::plot_ly()

connectivity_trace_group <- character()


# Bathymetry

connectivity_3d <- add_bathymetry(
  connectivity_3d,
  connectivity_bathymetry
)

connectivity_trace_group <- c(
  connectivity_trace_group,
  "base"
)


# Depth planes and labels

connectivity_3d <- add_depth_planes(
  connectivity_3d,
  connectivity_bathymetry
)

connectivity_trace_group <- c(
  connectivity_trace_group,
  rep(
    "base",
    5
  )
)


# Multilayer nodes

for (
  depth_current in depth_levels
) {
  
  node_data <- node_positions %>%
    filter(
      depth_stratum ==
        depth_current
    )
  
  connectivity_3d <- connectivity_3d %>%
    plotly::add_trace(
      data =
        node_data,
      
      x =
        ~lon,
      
      y =
        ~lat,
      
      z =
        ~z_layer,
      
      type =
        "scatter3d",
      
      mode =
        "markers",
      
      name =
        depth_current,
      
      text =
        ~node_hover,
      
      hoverinfo =
        "text",
      
      marker =
        list(
          size =
            5,
          
          color =
            depth_layer_colours[[
              depth_current
            ]],
          
          line =
            list(
              color =
                "#17324D",
              
              width =
                0.8
            )
        ),
      
      showlegend =
        TRUE
    )
  
  connectivity_trace_group <- c(
    connectivity_trace_group,
    "base"
  )
}


# PLD-specific connectivity edges

for (
  pld_current in pld_values
) {
  
  pld_edges <- connectivity_edges %>%
    filter(
      pld_days ==
        pld_current
    )
  
  for (
    strength_current in names(
      edge_widths
    )
  ) {
    
    edge_data <- pld_edges %>%
      filter(
        strength_class ==
          strength_current
      )
    
    if (
      nrow(
        edge_data
      ) ==
      0
    ) {
      next
    }
    
    edge_segments <- make_edge_segments(
      edge_data
    )
    
    connectivity_3d <- connectivity_3d %>%
      plotly::add_trace(
        data =
          edge_segments,
        
        x =
          ~lon,
        
        y =
          ~lat,
        
        z =
          ~z,
        
        type =
          "scatter3d",
        
        mode =
          "lines",
        
        text =
          ~hover_text,
        
        hoverinfo =
          "text",
        
        line =
          list(
            color =
              pld_colours[[
                as.character(
                  pld_current
                )
              ]],
            
            width =
              edge_widths[[
                strength_current
              ]]
          ),
        
        opacity =
          0.85,
        
        visible =
          pld_current ==
          pld_values[1],
        
        showlegend =
          FALSE
      )
    
    connectivity_trace_group <- c(
      connectivity_trace_group,
      as.character(
        pld_current
      )
    )
  }
  
  
  # Destination diamonds
  
  connectivity_3d <- connectivity_3d %>%
    plotly::add_trace(
      data =
        pld_edges,
      
      x =
        ~destination_lon,
      
      y =
        ~destination_lat,
      
      z =
        ~destination_z,
      
      type =
        "scatter3d",
      
      mode =
        "markers",
      
      text =
        ~edge_hover,
      
      hoverinfo =
        "text",
      
      marker =
        list(
          size =
            4,
          
          color =
            pld_colours[[
              as.character(
                pld_current
              )
            ]],
          
          symbol =
            "diamond",
          
          line =
            list(
              color =
                "white",
              
              width =
                0.7
            )
        ),
      
      visible =
        pld_current ==
        pld_values[1],
      
      showlegend =
        FALSE
    )
  
  connectivity_trace_group <- c(
    connectivity_trace_group,
    as.character(
      pld_current
    )
  )
}


# Connectivity PLD buttons

connectivity_buttons <- lapply(
  pld_values,
  function(
    pld_current) {
    
    list(
      method =
        "update",
      
      label =
        paste0(
          pld_current,
          " d"
        ),
      
      args =
        list(
          list(
            visible =
              connectivity_trace_group ==
              "base" |
              connectivity_trace_group ==
              as.character(
                pld_current
              )
          ),
          
          list(
            title =
              paste0(
                "Depth-resolved seamount connectivity — ",
                pld_current,
                " d PLD",
                
                "<br><sup>colour = PLD; line width = connectivity strength; diamonds = destination; blue planes = analytical depth boundaries</sup>"
              )
          )
        )
    )
  }
)


connectivity_3d <- connectivity_3d %>%
  plotly::layout(
    title =
      paste0(
        "Depth-resolved seamount connectivity — ",
        pld_values[1],
        " d PLD",
        
        "<br><sup>colour = PLD; line width = connectivity strength; diamonds = destination; blue planes = analytical depth boundaries</sup>"
      ),
    
    updatemenus =
      list(
        list(
          type =
            "buttons",
          
          direction =
            "right",
          
          x =
            0.02,
          
          y =
            1.10,
          
          xanchor =
            "left",
          
          yanchor =
            "top",
          
          buttons =
            connectivity_buttons,
          
          bgcolor =
            "rgba(255,255,255,0.85)",
          
          bordercolor =
            "#666666"
        )
      ),
    
    legend =
      list(
        title =
          list(
            text =
              "Depth stratum"
          ),
        
        yanchor =
          "top",
        
        y =
          0.55
      ),
    
    scene =
      list(
        xaxis =
          list(
            title =
              "Longitude"
          ),
        
        yaxis =
          list(
            title =
              "Latitude"
          ),
        
        zaxis =
          list(
            title =
              "Depth / elevation (m)"
          ),
        
        aspectmode =
          "manual",
        
        aspectratio =
          list(
            x =
              0.7,
            
            y =
              2,
            
            z =
              0.9
          ),
        
        camera =
          list(
            eye =
              list(
                x =
                  1.4,
                
                y =
                  -1.5,
                
                z =
                  1
              )
          )
      ),
    
    margin =
      list(
        b =
          40,
        
        l =
          60,
        
        t =
          80,
        
        r =
          10
      )
  ) %>%
  plotly::config(
    displaylogo =
      FALSE,
    
    responsive =
      TRUE
  )


htmlwidgets::saveWidget(
  connectivity_3d,
  "html/connectivity_3D.html",
  selfcontained = TRUE
)


# Display static figures

print(
  figure_within_between
)

print(
  figure_depth_heatmap
)

print(
  figure_chain
)

print(
  figure_roles
)

print(
  figure_connectivity_map
)

print(
  figure_connectivity_heatmap
)