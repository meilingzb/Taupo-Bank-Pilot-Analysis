# Study site figure 

library(terra)
library(sf)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

# Load current project inputs

source("00_setup_and_inputs_.R")


# Prepare included nodes

Lord_howe_nodes <- candidate_nodes %>%
  filter(
    chain == "Lord Howe",
    status == "include"
  ) %>%
  mutate(
    chain_plot = "Lord Howe"
  )

Tasmantid_nodes <- candidate_nodes %>%
  filter(
    chain == "Tasmantid",
    status == "include"
  ) %>%
  mutate(
    chain_plot = "Tasmantid"
  )

all_nodes <- bind_rows(
  Lord_howe_nodes,
  Tasmantid_nodes
)


# Load larger GEBCO file for study-site figure

map_gebco_path <- "gebco_2026_n-10.0_s-45.0_w145.0_e165.0.nc"

if (!file.exists(map_gebco_path)) {
  
  message(
    "Select the larger GEBCO file from the old project."
  )
  
  map_gebco_path <- file.choose()
}

map_gebco <- rast(
  map_gebco_path,
  subds = "elevation"
)


# Crop GEBCO to figure extent

map_extent <- ext(
  147,
  164,
  -38,
  -14
)

gebco_crop <- crop(
  map_gebco,
  map_extent
)


# Reduce resolution for plotting

gebco_small <- aggregate(
  gebco_crop,
  fact = 4,
  fun = mean,
  na.rm = TRUE
)


# Create land polygons

land_r <- ifel(
  gebco_small >= 0,
  1,
  NA
)

land_poly <- as.polygons(
  land_r,
  dissolve = TRUE,
  na.rm = TRUE
)

land_sf <- st_as_sf(
  land_poly
)


# Create bathymetric contours

contours_minor <- as.contour(
  gebco_small,
  levels = seq(
    -6000,
    0,
    by = 1000
  )
)

contours_major <- as.contour(
  gebco_small,
  levels = seq(
    -6000,
    0,
    by = 2000
  )
)

contours_minor_sf <- st_as_sf(
  contours_minor
)

contours_major_sf <- st_as_sf(
  contours_major
)


# Australian EEZ

eez_file <- "australian_eez.rds"

if (file.exists(eez_file)) {
  
  australian_eez <- readRDS(
    eez_file
  )
  
} else {
  
  australian_eez <- st_read(
    paste0(
      "https://geo.vliz.be/geoserver/MarineRegions/wfs?",
      "service=WFS&",
      "version=1.0.0&",
      "request=GetFeature&",
      "typeName=eez&",
      "cql_filter=mrgid=8323&",
      "outputFormat=application/json"
    ),
    quiet = TRUE
  ) %>%
    st_make_valid() %>%
    st_transform(
      4326
    )
  
  saveRDS(
    australian_eez,
    eez_file
  )
}


# Crop EEZ to figure extent

eez_bbox <- st_bbox(
  c(
    xmin = 147,
    xmax = 164,
    ymin = -38,
    ymax = -14
  ),
  crs = st_crs(4326)
)

australian_eez_crop <- st_crop(
  australian_eez,
  eez_bbox
)


# Prepare site lists

tasm_table_data <- Tasmantid_nodes %>%
  arrange(
    node_id
  ) %>%
  transmute(
    label = paste(
      node_id,
      name
    )
  )

lord_table_data <- Lord_howe_nodes %>%
  arrange(
    node_id
  ) %>%
  transmute(
    label = paste(
      node_id,
      name
    )
  )

tasm_table_data$y <- rev(
  seq_len(
    nrow(
      tasm_table_data
    )
  )
)

lord_table_data$y <- rev(
  seq_len(
    nrow(
      lord_table_data
    )
  )
)


# Create map

map_plot <- ggplot() +
  
  # Land
  geom_sf(
    data = land_sf,
    fill = "grey92",
    colour = "grey35",
    linewidth = 0.4
  ) +
  
  # Minor bathymetric contours
  geom_sf(
    data = contours_minor_sf,
    colour = "grey75",
    linewidth = 0.20,
    alpha = 0.45
  ) +
  
  # Major bathymetric contours
  geom_sf(
    data = contours_major_sf,
    colour = "grey50",
    linewidth = 0.40,
    alpha = 0.60
  ) +
  
  # Australian EEZ
  geom_sf(
    data = australian_eez_crop,
    aes(
      linetype = "Australian EEZ"
    ),
    fill = NA,
    colour = "grey40",
    linewidth = 0.55,
    inherit.aes = FALSE,
    show.legend = "line"
  ) +
  
  # Tasmantid sites
  geom_point(
    data = Tasmantid_nodes,
    aes(
      x = lon,
      y = lat,
      shape = chain_plot
    ),
    size = 2.8,
    colour = "black",
    fill = "black",
    stroke = 0.5
  ) +
  
  # Lord Howe sites
  geom_point(
    data = Lord_howe_nodes,
    aes(
      x = lon,
      y = lat,
      shape = chain_plot
    ),
    size = 3.6,
    colour = "black",
    fill = "black",
    stroke = 0.5
  ) +
  
  # Tasmantid labels with leader lines
  geom_text_repel(
    data = Tasmantid_nodes,
    aes(
      x = lon,
      y = lat,
      label = node_id
    ),
    nudge_x = -0.75,
    direction = "y",
    hjust = 1,
    size = 3.1,
    fontface = "bold",
    colour = "black",
    box.padding = 0.12,
    point.padding = 0.20,
    segment.color = "grey60",
    segment.size = 0.25,
    min.segment.length = 0,
    force = 2.5,
    max.iter = 10000,
    max.overlaps = Inf,
    seed = 1
  ) +
  
  # Lord Howe labels with leader lines
  geom_text_repel(
    data = Lord_howe_nodes,
    aes(
      x = lon,
      y = lat,
      label = node_id
    ),
    nudge_x = 0.65,
    direction = "y",
    hjust = 0,
    size = 3.1,
    fontface = "bold",
    colour = "black",
    box.padding = 0.12,
    point.padding = 0.18,
    segment.color = "grey60",
    segment.size = 0.25,
    min.segment.length = 0,
    force = 2,
    max.iter = 10000,
    max.overlaps = Inf,
    seed = 2
  ) +
  
  # Chain symbols
  scale_shape_manual(
    values = c(
      "Lord Howe" = 17,
      "Tasmantid" = 16
    ),
    name = "Chain"
  ) +
  
  # EEZ line style
  scale_linetype_manual(
    values = c(
      "Australian EEZ" = "dashed"
    ),
    name = "EEZ"
  ) +
  
  # Keep shape and line legends separate
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        linetype = 0
      )
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        shape = NA
      )
    )
  ) +
  
  coord_sf(
    xlim = c(
      147,
      164
    ),
    ylim = c(
      -38,
      -14
    ),
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(
      148,
      152,
      156,
      160,
      164
    ),
    labels = c(
      "148°E",
      "152°E",
      "156°E",
      "160°E",
      "164°E"
    )
  ) +
  
  scale_y_continuous(
    breaks = c(
      -15,
      -20,
      -25,
      -30,
      -35
    ),
    labels = c(
      "15°S",
      "20°S",
      "25°S",
      "30°S",
      "35°S"
    )
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.55
    ),
    
    axis.title = element_text(
      size = 16
    ),
    
    axis.text = element_text(
      size = 14,
      colour = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 13
    ),
    
    legend.key.width = grid::unit(
      1,
      "cm"
    ),
    
    plot.margin = margin(
      10,
      10,
      10,
      10
    )
  )


# Create Tasmantid site list

tasm_table_plot <- ggplot(
  tasm_table_data,
  aes(
    x = 0,
    y = y,
    label = label
  )
) +
  
  geom_text(
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = max(
      tasm_table_data$y
    ) + 1.3,
    label = "Tasmantid sites",
    hjust = 0,
    fontface = "bold",
    size = 5
  ) +
  
  xlim(
    0,
    1
  ) +
  
  ylim(
    0.5,
    max(
      tasm_table_data$y
    ) + 2
  ) +
  
  theme_void() +
  
  theme(
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )


# Create Lord Howe site list

lord_table_plot <- ggplot(
  lord_table_data,
  aes(
    x = 0,
    y = y,
    label = label
  )
) +
  
  geom_text(
    hjust = 0,
    size = 4
  ) +
  
  annotate(
    "text",
    x = 0,
    y = max(
      lord_table_data$y
    ) + 1.3,
    label = "Lord Howe sites",
    hjust = 0,
    fontface = "bold",
    size = 5
  ) +
  
  xlim(
    0,
    1
  ) +
  
  ylim(
    0.5,
    max(
      lord_table_data$y
    ) + 2
  ) +
  
  theme_void() +
  
  theme(
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )


# Combine site lists

table_plot <- tasm_table_plot /
  lord_table_plot +
  
  plot_layout(
    heights = c(
      2.1,
      1.3
    )
  )


# Combine map and site lists

final_plot <- (
  map_plot |
    table_plot
) +
  
  plot_layout(
    widths = c(
      1.3,
      1
    ),
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )


# View figure

final_plot


# Save PDF

ggsave(
  "study_site_figure_final_EEZ.pdf",
  plot = final_plot,
  width = 14,
  height = 10,
  units = "in",
  bg = "white"
)


# Save PNG

ggsave(
  "study_site_figure_final_EEZ.png",
  plot = final_plot,
  width = 14,
  height = 10,
  units = "in",
  dpi = 400,
  bg = "white"
)
