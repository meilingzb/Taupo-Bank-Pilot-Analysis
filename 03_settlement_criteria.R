# Taupo Bank Analysis pt. 4 :)
# Depth-resolved potential habitat encounters
#
# Potential settlement events are identified for competent particles using
# the same bathymetry and grid resolution represented in the LLC4320v2
# circulation model.
#
# The GEBCO-derived NETVOLC polygons from 02 are retained as the horizontal
# seamount footprints. They define which part of the model seabed belongs to
# each seamount, while LLC4320v2 supplies the benthic depth used to evaluate
# the particle trajectories.
#
# A particle position qualifies as a potential settlement event when:
#
# 1. the particle is biologically competent and has not exceeded its PLD;
# 2. its horizontal distance from a mapped benthic patch is within one
#    area-equivalent local LLC4320v2 model-cell width;
# 3. the particle is above the benthic surface; and
# 4. its height above the benthic surface is within one LLC4320v2 vertical
#    model layer at that benthic depth.
#
# If more than one benthic surface qualifies at the same timestep, the nearest
# surface in 3-D Euclidean distance is retained. Only the first qualifying
# encounter is then kept for each particle and PLD scenario.
#
# This was updated from the previous Euclidean distance method as biologically 
# we could not justify any reasonable encounter distance thresholds, so this
# method basically considers "settlement" as being within the model-resolved 
# grid cell of the benthic surface on the seamount.
# 
# This is to do with the fact that we don't know any sort of encounter distances
# required for deep-sea coral (e.g. if they can sense when they're near a seamount), 
# but we know of course contact with the benthos needs to occur for settlement.
# This does have limitations, such as the fact the grid cells change size with 
# depth/horizontal distance, however, it is more biologically reasonable and 
# represents benthic encounterment rather than larval sensing distance/thresholds. 

library(sf)
library(dplyr)
library(tibble)
library(ncdf4)


# Load project inputs

candidate_nodes <- readRDS(
  file.path(
    "Outputs",
    "candidate_nodes.rds"
  )
)

final_footprints <- readRDS(
  file.path(
    "Outputs",
    "final_seamount_footprints.rds"
  )
)


# Raw model files

particle_file <- file.path(
  "Data",
  "tasman_seamount_2023_merged.nc"
)

horizontal_grid_file <- file.path(
  "Data",
  "LLC4320v2_tasman_horizontal_grid.nc"
)

vertical_grid_file <- file.path(
  "Data",
  "LLC4320v2_vertical_grid.nc"
)


# Current pilot source
#
# Source seamount is defined by Dustin's release simulation rather than
# inferred from the later NETVOLC footprint.

source_node_current <- "T22"


# Depth strata
#
# Depth is expressed as positive metres below sea level.
# The upper 200 m is excluded from the depth-resolved network.

depth_breaks_m <- c(
  200,
  700,
  1500,
  2500,
  Inf
)

depth_labels <- c(
  "200-700 m",
  "700-1500 m",
  "1500-2500 m",
  ">2500 m"
)


assign_depth_stratum <- function(depth_m) {
  
  cut(
    depth_m,
    breaks = depth_breaks_m,
    labels = depth_labels,
    right = FALSE,
    include.lowest = TRUE
  )
}


# Biological scenarios
#
# The full biological design crosses three pelagic larval durations (PLDs)
# with three pre-competency periods. This gives nine biological scenarios:
#
# PLD = 30, 60, 90 days
# pre-competency = 14, 21, 28 days
#
# Pre-competency is treated independently of PLD so its effect on potential
# settlement can be examined without changing the total larval duration.

biological_scenarios <- dplyr::cross_join(
  tibble(
    pld_days = c(
      30,
      60,
      90
    )
  ),
  tibble(
    precompetency_days = c(
      14,
      21,
      28
    )
  )
) %>%
  arrange(
    pld_days,
    precompetency_days
  )


# LLC4320v2 horizontal grid
#
# Depth gives the model bathymetry and RAC gives the horizontal tracer-cell
# area. sqrt(RAC) is used as an area-equivalent local model-cell width.

orient_model_matrix <- function(
    x,
    n_j,
    n_i) {
  
  if (
    length(dim(x)) == 2 &&
    all(dim(x) == c(n_j, n_i))
  ) {
    return(x)
  }
  
  if (
    length(dim(x)) == 2 &&
    all(dim(x) == c(n_i, n_j))
  ) {
    return(t(x))
  }
  
  stop(
    "Unexpected LLC4320v2 grid dimensions."
  )
}


horizontal_nc <- ncdf4::nc_open(
  horizontal_grid_file
)

n_j <- horizontal_nc$dim$j$len
n_i <- horizontal_nc$dim$i$len

XC <- orient_model_matrix(
  ncdf4::ncvar_get(
    horizontal_nc,
    "XC"
  ),
  n_j,
  n_i
)

YC <- orient_model_matrix(
  ncdf4::ncvar_get(
    horizontal_nc,
    "YC"
  ),
  n_j,
  n_i
)

model_depth_m <- orient_model_matrix(
  ncdf4::ncvar_get(
    horizontal_nc,
    "Depth"
  ),
  n_j,
  n_i
)

cell_area_m2 <- orient_model_matrix(
  ncdf4::ncvar_get(
    horizontal_nc,
    "RAC"
  ),
  n_j,
  n_i
)

ncdf4::nc_close(
  horizontal_nc
)

# Dustin's Tasman Sea subset is rectilinear in longitude and latitude, so
# the 2-D centre coordinates can be represented by one longitude vector and
# one latitude vector when constructing local model-cell polygons.

llc_longitude <- as.numeric(
  XC[1, ]
)

llc_latitude <- as.numeric(
  YC[, 1]
)

rm(
  XC,
  YC
)


# Convert model-cell centres to cell edges.
#
# Interior edges lie halfway between neighbouring centres. The two outer
# edges are extended by half of the nearest centre spacing.

centres_to_edges <- function(centres) {
  
  interior_edges <- (
    centres[-1] +
      centres[-length(centres)]
  ) / 2
  
  first_edge <-
    centres[1] -
    (interior_edges[1] - centres[1])
  
  last_edge <-
    centres[length(centres)] +
    (centres[length(centres)] -
       interior_edges[length(interior_edges)])
  
  c(
    first_edge,
    interior_edges,
    last_edge
  )
}


llc_longitude_edges <- centres_to_edges(
  llc_longitude
)

llc_latitude_edges <- centres_to_edges(
  llc_latitude
)


# LLC4320v2 vertical grid
#
# DRF gives model-layer thickness and depth_F gives layer-interface depths.
# The layer thickness at each benthic depth is used as the vertical encounter
# tolerance.

vertical_nc <- ncdf4::nc_open(
  vertical_grid_file
)

llc_layer_width_m <- abs(
  as.vector(
    ncdf4::ncvar_get(
      vertical_nc,
      "DRF"
    )
  )
)

llc_depth_interfaces_m <- abs(
  as.vector(
    ncdf4::ncvar_get(
      vertical_nc,
      "depth_F"
    )
  )
)

ncdf4::nc_close(
  vertical_nc
)

vertical_grid <- tibble(
  layer = seq_along(
    llc_layer_width_m
  ),
  top_depth_m =
    llc_depth_interfaces_m[
      -length(llc_depth_interfaces_m)
    ],
  bottom_depth_m =
    llc_depth_interfaces_m[-1],
  layer_width_m =
    llc_layer_width_m
)


get_vertical_resolution_at_depth <- function(
    depth_m,
    vertical_grid) {
  
  layer_index <- findInterval(
    depth_m,
    vertical_grid$top_depth_m
  )
  
  layer_width_m <- rep(
    NA_real_,
    length(depth_m)
  )
  
  valid <-
    layer_index >= 1 &
    layer_index <= nrow(vertical_grid)
  
  layer_width_m[valid] <-
    vertical_grid$layer_width_m[
      layer_index[valid]
    ]
  
  layer_width_m
}


# Local projected coordinate system

get_local_utm_epsg <- function(lon, lat) {
  
  zone <- floor(
    (lon + 180) / 6
  ) + 1
  
  if (lat < 0) {
    32700 + zone
  } else {
    32600 + zone
  }
}


# Build the benthic surface for one seamount.
#
# The GEBCO-derived NETVOLC polygon defines the horizontal footprint of the
# seamount. LLC4320v2 cells intersecting that footprint provide the benthic
# depth used by the particle model and the local horizontal model resolution.

build_benthic_surface <- function(
    footprint,
    node_lon,
    node_lat,
    minimum_depth_m = 200) {
  
  footprint_ll <- st_transform(
    footprint,
    4326
  )
  
  footprint_geometry <- st_union(
    st_geometry(
      footprint_ll
    )
  )
  
  footprint_bbox <- st_bbox(
    footprint_geometry
  )
  
  lon_left <- pmin(
    llc_longitude_edges[
      -length(llc_longitude_edges)
    ],
    llc_longitude_edges[-1]
  )
  
  lon_right <- pmax(
    llc_longitude_edges[
      -length(llc_longitude_edges)
    ],
    llc_longitude_edges[-1]
  )
  
  lat_bottom <- pmin(
    llc_latitude_edges[
      -length(llc_latitude_edges)
    ],
    llc_latitude_edges[-1]
  )
  
  lat_top <- pmax(
    llc_latitude_edges[
      -length(llc_latitude_edges)
    ],
    llc_latitude_edges[-1]
  )
  
  i_index <- which(
    lon_left <= footprint_bbox["xmax"] &
      lon_right >= footprint_bbox["xmin"]
  )
  
  j_index <- which(
    lat_bottom <= footprint_bbox["ymax"] &
      lat_top >= footprint_bbox["ymin"]
  )
  
  cell_index <- expand.grid(
    global_j = j_index,
    global_i = i_index,
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    as_tibble()
  
  cell_index$depth_m <- model_depth_m[
    cbind(
      cell_index$global_j,
      cell_index$global_i
    )
  ]
  
  cell_index$cell_area_m2 <- cell_area_m2[
    cbind(
      cell_index$global_j,
      cell_index$global_i
    )
  ]
  
  cell_index <- cell_index %>%
    filter(
      is.finite(depth_m),
      depth_m >= minimum_depth_m,
      is.finite(cell_area_m2),
      cell_area_m2 > 0
    )
  
  cell_polygons <- lapply(
    seq_len(nrow(cell_index)),
    function(i) {
      
      grid_i <- cell_index$global_i[i]
      grid_j <- cell_index$global_j[i]
      
      xmin <- min(
        llc_longitude_edges[
          c(grid_i, grid_i + 1)
        ]
      )
      
      xmax <- max(
        llc_longitude_edges[
          c(grid_i, grid_i + 1)
        ]
      )
      
      ymin <- min(
        llc_latitude_edges[
          c(grid_j, grid_j + 1)
        ]
      )
      
      ymax <- max(
        llc_latitude_edges[
          c(grid_j, grid_j + 1)
        ]
      )
      
      st_polygon(
        list(
          matrix(
            c(
              xmin, ymin,
              xmax, ymin,
              xmax, ymax,
              xmin, ymax,
              xmin, ymin
            ),
            ncol = 2,
            byrow = TRUE
          )
        )
      )
    }
  )
  
  surface_cells <- st_sf(
    cell_index,
    geometry = st_sfc(
      cell_polygons,
      crs = 4326
    )
  )
  
  surface_cells <- surface_cells[
    lengths(
      st_intersects(
        surface_cells,
        footprint_geometry
      )
    ) > 0,
  ]
  
  # Clip cells at the edge of the seamount so horizontal distance is measured
  # to the part of the model cell that actually lies within the NETVOLC footprint.
  
  surface_cells <- suppressWarnings(
    st_intersection(
      surface_cells,
      footprint_geometry
    )
  )
  
  surface_cells <- surface_cells[
    !st_is_empty(surface_cells),
  ]
  
  # One local horizontal model-cell width.
  #
  # RAC is the horizontal area of each tracer cell. Its square root provides
  # an area-equivalent characteristic cell width in metres.
  
  surface_cells$horizontal_tolerance_m <- sqrt(
    surface_cells$cell_area_m2
  )
  
  surface_cells$vertical_tolerance_m <-
    get_vertical_resolution_at_depth(
      surface_cells$depth_m,
      vertical_grid
    )
  
  surface_cells$depth_stratum <-
    as.character(
      assign_depth_stratum(
        surface_cells$depth_m
      )
    )
  
  surface_cells <- surface_cells %>%
    filter(
      is.finite(horizontal_tolerance_m),
      is.finite(vertical_tolerance_m),
      !is.na(depth_stratum)
    )
  
  surface_cells$surface_cell_id <-
    seq_len(
      nrow(surface_cells)
    )
  
  # Work in metres for distance calculations.
  
  utm_epsg <- get_local_utm_epsg(
    node_lon,
    node_lat
  )
  
  st_transform(
    surface_cells,
    utm_epsg
  )
}


# Build benthic surfaces for all seamounts.

build_all_benthic_surfaces <- function(
    final_footprints,
    candidate_nodes) {
  
  surface_list <- vector(
    "list",
    nrow(final_footprints)
  )
  
  names(surface_list) <-
    final_footprints$node_id
  
  for (i in seq_len(nrow(final_footprints))) {
    
    node_id_current <-
      final_footprints$node_id[i]
    
    node <- candidate_nodes %>%
      filter(
        node_id == node_id_current
      )
    
    surface_list[[node_id_current]] <-
      build_benthic_surface(
        footprint =
          final_footprints[i, ],
        
        node_lon =
          node$lon,
        
        node_lat =
          node$lat
      )
  }
  
  surface_list
}


# Habitat-encounter search zones
#
# These buffers are only a computational pre-screen. Each footprint is
# expanded by the largest local model-cell width on that seamount. The final
# encounter test still uses the tolerance attached to each individual cell.

build_settlement_search_zones <- function(
    final_footprints,
    candidate_nodes,
    surface_list) {
  
  zone_list <- vector(
    "list",
    nrow(final_footprints)
  )
  
  for (i in seq_len(nrow(final_footprints))) {
    
    node_id_current <-
      final_footprints$node_id[i]
    
    node <- candidate_nodes %>%
      filter(
        node_id == node_id_current
      )
    
    maximum_tolerance_m <- max(
      surface_list[[node_id_current]]$horizontal_tolerance_m,
      na.rm = TRUE
    )
    
    utm_epsg <- get_local_utm_epsg(
      node$lon,
      node$lat
    )
    
    footprint_m <- final_footprints[i, ] %>%
      st_transform(
        utm_epsg
      )
    
    zone_m <- st_buffer(
      st_union(
        st_geometry(
          footprint_m
        )
      ),
      dist = maximum_tolerance_m
    )
    
    zone_list[[i]] <- st_sf(
      node_id =
        node_id_current,
      
      geometry =
        st_transform(
          zone_m,
          4326
        )
    )
  }
  
  do.call(
    rbind,
    zone_list
  )
}


# Network node inventory
#
# Each seamount-depth combination containing benthic habitat
# becomes one node in the multilayer network.

build_network_node_table <- function(
    surface_list,
    candidate_nodes) {
  
  node_tables <- lapply(
    names(surface_list),
    function(node_id_current) {
      
      surface_list[[node_id_current]] %>%
        st_drop_geometry() %>%
        distinct(
          depth_stratum
        ) %>%
        transmute(
          node_id =
            node_id_current,
          
          depth_stratum =
            as.character(
              depth_stratum
            )
        )
    }
  )
  
  
  bind_rows(
    node_tables
  ) %>%
    left_join(
      candidate_nodes %>%
        select(
          node_id,
          name,
          chain
        ),
      by = "node_id"
    ) %>%
    arrange(
      node_id,
      match(
        depth_stratum,
        depth_labels
      )
    )
}


# Read particle trajectories

read_particle_trajectories <- function(
    particle_file) {
  
  particle_nc <- ncdf4::nc_open(
    particle_file
  )
  
  on.exit(
    ncdf4::nc_close(
      particle_nc
    )
  )
  
  
  n_obs <-
    particle_nc$dim$obs$len
  
  n_particles <-
    particle_nc$dim$particle$len
  
  
  lon <- ncdf4::ncvar_get(
    particle_nc,
    "lon"
  )
  
  lat <- ncdf4::ncvar_get(
    particle_nc,
    "lat"
  )
  
  depth <- ncdf4::ncvar_get(
    particle_nc,
    "depth"
  )
  
  
  particle_id <- as.vector(
    ncdf4::ncvar_get(
      particle_nc,
      "merged_id"
    )
  )
  
  release_id <- as.vector(
    ncdf4::ncvar_get(
      particle_nc,
      "release_id"
    )
  )
  
  above_summit <- as.vector(
    ncdf4::ncvar_get(
      particle_nc,
      "above_summit"
    )
  )
  
  
  age_raw <- as.vector(
    ncdf4::ncvar_get(
      particle_nc,
      "age_hours"
    )
  )
  
  
  # age_hours may be shared between particles or stored
  # separately for every observation.
  
  if (length(age_raw) == n_obs) {
    
    age_matrix <- matrix(
      rep(
        age_raw,
        n_particles
      ),
      nrow = n_obs,
      ncol = n_particles
    )
    
  } else {
    
    age_matrix <- matrix(
      age_raw,
      nrow = n_obs,
      ncol = n_particles
    )
  }
  
  
  trajectories <- tibble(
    
    particle_id =
      rep(
        particle_id,
        each = n_obs
      ),
    
    observation_index =
      rep(
        seq_len(n_obs),
        times = n_particles
      ),
    
    age_hours =
      as.vector(
        age_matrix
      ),
    
    lon =
      as.vector(
        lon
      ),
    
    lat =
      as.vector(
        lat
      ),
    
    depth_m =
      abs(
        as.vector(
          depth
        )
      ),
    
    above_summit =
      rep(
        above_summit,
        each = n_obs
      ),
    
    release_id =
      rep(
        release_id,
        each = n_obs
      )
  )
  
  
  trajectories %>%
    filter(
      above_summit == 0,
      is.finite(age_hours),
      is.finite(lon),
      is.finite(lat),
      is.finite(depth_m)
    )
}


# Assign source node and source depth.
#
# Source node comes from Dustin's release simulation.
# Only particles released at 200 m depth or deeper are retained.
# Release positions are not filtered using the NETVOLC footprint.

assign_particle_sources <- function(
    trajectory_data,
    source_node,
    minimum_depth_m = 200) {
  
  release_data <- trajectory_data %>%
    group_by(
      particle_id
    ) %>%
    slice_min(
      order_by = age_hours,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    filter(
      depth_m >= minimum_depth_m
    ) %>%
    mutate(
      source_node =
        source_node,
      
      source_depth_m =
        depth_m,
      
      source_depth_stratum =
        as.character(
          assign_depth_stratum(
            source_depth_m
          )
        )
    )
  
  
  # Remove complete trajectories for releases shallower than 200 m.
  
  trajectory_data %>%
    filter(
      particle_id %in%
        release_data$particle_id
    ) %>%
    left_join(
      release_data %>%
        select(
          particle_id,
          source_node,
          source_depth_m,
          source_depth_stratum
        ),
      by = "particle_id"
    )
}


# Apply the PLD and pre-competency scenarios.

prepare_biological_trajectories <- function(
    trajectory_data,
    biological_scenarios) {
  
  trajectory_data <- trajectory_data %>%
    mutate(
      particle_age_days =
        age_hours / 24
    )
  
  
  trajectory_data <- dplyr::cross_join(
    trajectory_data,
    biological_scenarios
  )
  
  
  trajectory_data %>%
    mutate(
      settlement_eligible =
        particle_age_days >=
        precompetency_days &
        particle_age_days <=
        pld_days
    ) %>%
    filter(
      particle_age_days >= 0,
      particle_age_days <= pld_days,
      !is.na(
        source_depth_stratum
      )
    )
}


# Number of particles released from each source node and depth.
#
# These counts form the denominator of the connectivity probability.

calculate_release_counts <- function(
    biological_trajectories) {
  
  biological_trajectories %>%
    distinct(
      particle_id,
      release_id,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    count(
      release_id,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum,
      name = "N_released"
    )
}

# Potential settlement events
#
# Each unique competent particle position is tested against the model benthic
# surface once. PLD and competency are then used below to retain the first
# qualifying encounter under each biological scenario.

calculate_position_encounters <- function(
    biological_trajectories,
    surface_list,
    search_zones) {
  
  positions <- biological_trajectories %>%
    filter(
      settlement_eligible
    ) %>%
    distinct(
      particle_id,
      release_id,
      observation_index,
      lon,
      lat,
      depth_m
    ) %>%
    mutate(
      position_id =
        row_number()
    )
  
  points_ll <- st_as_sf(
    positions,
    coords = c(
      "lon",
      "lat"
    ),
    crs = 4326,
    remove = FALSE
  )
  
  # Use the buffered footprints to identify which seamount surfaces need to
  # be checked for each particle position.
  
  zone_hits <- st_intersects(
    points_ll,
    search_zones
  )
  
  candidate_pairs <- tibble(
    position_id =
      rep(
        positions$position_id,
        lengths(zone_hits)
      ),
    
    zone_index =
      unlist(
        zone_hits,
        use.names = FALSE
      )
  )
  
  if (nrow(candidate_pairs) == 0) {
    return(
      tibble()
    )
  }
  
  candidate_pairs <- candidate_pairs %>%
    mutate(
      destination_node =
        search_zones$node_id[
          zone_index
        ]
    )
  
  encounter_list <- vector(
    "list",
    length(surface_list)
  )
  
  node_ids <- names(
    surface_list
  )
  
  for (n in seq_along(node_ids)) {
    
    node_id_current <-
      node_ids[n]
    
    node_positions <- candidate_pairs %>%
      filter(
        destination_node ==
          node_id_current
      ) %>%
      pull(
        position_id
      )
    
    if (length(node_positions) == 0) {
      next
    }
    
    surface_cells <-
      surface_list[[
        node_id_current
      ]]
    
    position_rows <- match(
      node_positions,
      positions$position_id
    )
    
    points_m <- points_ll[
      position_rows,
    ] %>%
      st_transform(
        st_crs(
          surface_cells
        )
      )
    
    maximum_tolerance_m <- max(
      surface_cells$horizontal_tolerance_m,
      na.rm = TRUE
    )
    
    nearby_cells <- st_is_within_distance(
      points_m,
      surface_cells,
      dist = maximum_tolerance_m
    )
    
    pair_table <- tibble(
      point_index =
        rep(
          seq_along(nearby_cells),
          lengths(nearby_cells)
        ),
      
      cell_index =
        unlist(
          nearby_cells,
          use.names = FALSE
        )
    )
    
    if (nrow(pair_table) == 0) {
      next
    }
    
    horizontal_distance_m <- as.numeric(
      st_distance(
        points_m[
          pair_table$point_index,
        ],
        surface_cells[
          pair_table$cell_index,
        ],
        by_element = TRUE
      )
    )
    
    particle_depth_m <-
      positions$depth_m[
        position_rows[
          pair_table$point_index
        ]
      ]
    
    surface_depth_m <-
      surface_cells$depth_m[
        pair_table$cell_index
      ]
    
    # Positive values mean the particle is above the seabed. Negative values
    # would place the particle below the mapped benthic surface and are not
    # considered valid habitat encounters.
    
    height_above_benthos_m <-
      surface_depth_m -
      particle_depth_m
    
    vertical_distance_m <- abs(
      height_above_benthos_m
    )
    
    horizontal_tolerance_m <-
      surface_cells$horizontal_tolerance_m[
        pair_table$cell_index
      ]
    
    vertical_tolerance_m <-
      surface_cells$vertical_tolerance_m[
        pair_table$cell_index
      ]
    
    distance_3d_m <- sqrt(
      horizontal_distance_m^2 +
        vertical_distance_m^2
    )
    
    pair_table <- pair_table %>%
      mutate(
        position_id =
          node_positions[
            point_index
          ],
        
        destination_node =
          node_id_current,
        
        settlement_depth_m =
          surface_depth_m,
        
        destination_depth_stratum =
          as.character(
            surface_cells$depth_stratum[
              cell_index
            ]
          ),
        
        horizontal_distance_m =
          horizontal_distance_m,
        
        horizontal_tolerance_m =
          horizontal_tolerance_m,
        
        vertical_distance_m =
          vertical_distance_m,
        
        height_above_benthos_m =
          height_above_benthos_m,
        
        vertical_tolerance_m =
          vertical_tolerance_m,
        
        distance_3d_m =
          distance_3d_m,
        
        # Keep the original distance_m field name for compatibility with
        # existing downstream code. It is the full 3-D distance.
        
        distance_m =
          distance_3d_m,
        
        surface_cell_id =
          surface_cells$surface_cell_id[
            cell_index
          ]
      ) %>%
      filter(
        horizontal_distance_m <=
          horizontal_tolerance_m,
        
        height_above_benthos_m >= 0,
        
        height_above_benthos_m <=
          vertical_tolerance_m
      )
    
    if (nrow(pair_table) == 0) {
      next
    }
    
    # If several benthic cells on the same seamount qualify, retain the
    # nearest surface in 3-D.
    
    encounter_list[[n]] <- pair_table %>%
      arrange(
        position_id,
        distance_3d_m,
        surface_cell_id
      ) %>%
      group_by(
        position_id
      ) %>%
      slice_head(
        n = 1
      ) %>%
      ungroup()
  }
  
  encounters <- bind_rows(
    encounter_list
  )
  
  if (nrow(encounters) == 0) {
    return(
      tibble()
    )
  }
  
  # If more than one seamount qualifies at the same timestep, retain the
  # nearest qualifying benthic surface in 3-D.
  
  encounters %>%
    arrange(
      position_id,
      distance_3d_m,
      destination_node
    ) %>%
    group_by(
      position_id
    ) %>%
    slice_head(
      n = 1
    ) %>%
    ungroup() %>%
    left_join(
      positions %>%
        select(
          position_id,
          particle_id,
          release_id,
          observation_index
        ),
      by = "position_id"
    )
}


# Empty potential-settlement-event table

empty_settlement_events <- function() {
  
  tibble(
    particle_id = numeric(),
    release_id = numeric(),
    pld_days = numeric(),
    precompetency_days = numeric(),
    source_node = character(),
    source_depth_m = numeric(),
    source_depth_stratum = character(),
    settlement_age_days = numeric(),
    settlement_lon = numeric(),
    settlement_lat = numeric(),
    particle_depth_m = numeric(),
    destination_node = character(),
    settlement_depth_m = numeric(),
    destination_depth_stratum = character(),
    horizontal_distance_m = numeric(),
    horizontal_tolerance_m = numeric(),
    vertical_distance_m = numeric(),
    height_above_benthos_m = numeric(),
    vertical_tolerance_m = numeric(),
    distance_3d_m = numeric(),
    distance_m = numeric(),
    surface_cell_id = integer()
  )
}


# First qualifying potential-settlement-event rule
#
# Geometric encounters are calculated once for each particle position. PLD
# and competency are then used to identify the first qualifying encounter
# under each biological scenario. No persistence criterion is applied.

calculate_settlement_events <- function(
    biological_trajectories,
    surface_list,
    search_zones) {
  
  position_encounters <- calculate_position_encounters(
    biological_trajectories =
      biological_trajectories,
    
    surface_list =
      surface_list,
    
    search_zones =
      search_zones
  )
  
  if (nrow(position_encounters) == 0) {
    return(
      empty_settlement_events()
    )
  }
  
  biological_trajectories %>%
    filter(
      settlement_eligible
    ) %>%
    left_join(
      position_encounters %>%
        select(
          particle_id,
          release_id,
          observation_index,
          destination_node,
          settlement_depth_m,
          destination_depth_stratum,
          horizontal_distance_m,
          horizontal_tolerance_m,
          vertical_distance_m,
          height_above_benthos_m,
          vertical_tolerance_m,
          distance_3d_m,
          distance_m,
          surface_cell_id
        ),
      by = c(
        "particle_id",
        "release_id",
        "observation_index"
      )
    ) %>%
    filter(
      !is.na(
        destination_node
      )
    ) %>%
    arrange(
      particle_id,
      release_id,
      pld_days,
      precompetency_days,
      particle_age_days,
      observation_index
    ) %>%
    group_by(
      particle_id,
      release_id,
      pld_days,
      precompetency_days
    ) %>%
    slice_head(
      n = 1
    ) %>%
    ungroup() %>%
    transmute(
      particle_id,
      release_id,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_m,
      source_depth_stratum,
      
      settlement_age_days =
        particle_age_days,
      
      settlement_lon =
        lon,
      
      settlement_lat =
        lat,
      
      particle_depth_m =
        depth_m,
      
      destination_node,
      settlement_depth_m,
      destination_depth_stratum,
      horizontal_distance_m,
      horizontal_tolerance_m,
      vertical_distance_m,
      height_above_benthos_m,
      vertical_tolerance_m,
      distance_3d_m,
      distance_m,
      surface_cell_id
    )
}


# Build benthic surfaces

benthic_surfaces <- build_all_benthic_surfaces(
  final_footprints =
    final_footprints,
  
  candidate_nodes =
    candidate_nodes
)


# Build horizontal search zones

settlement_search_zones <-
  build_settlement_search_zones(
    final_footprints =
      final_footprints,
    
    candidate_nodes =
      candidate_nodes,
    
    surface_list =
      benthic_surfaces
  )


# Create seamount-depth nodes

network_nodes <- build_network_node_table(
  surface_list =
    benthic_surfaces,
  
  candidate_nodes =
    candidate_nodes
)


# Read particle trajectories

particle_trajectories <- read_particle_trajectories(
  particle_file
)


# Add the known source node and remove releases shallower than 200 m

particle_trajectories <- assign_particle_sources(
  trajectory_data =
    particle_trajectories,
  
  source_node =
    source_node_current
)


# Apply PLD and competency scenarios

biological_trajectories <-
  prepare_biological_trajectories(
    trajectory_data =
      particle_trajectories,
    
    biological_scenarios =
      biological_scenarios
  )


# Count released particles

release_counts <- calculate_release_counts(
  biological_trajectories
)


# Find first qualifying potential settlement events

settlement_events <- calculate_settlement_events(
  biological_trajectories =
    biological_trajectories,
  
  surface_list =
    benthic_surfaces,
  
  search_zones =
    settlement_search_zones
)


# Save outputs used by the connectivity analysis

saveRDS(
  network_nodes,
  file.path(
    "Outputs",
    "network_nodes.rds"
  )
)

saveRDS(
  release_counts,
  file.path(
    "Outputs",
    "release_counts.rds"
  )
)

saveRDS(
  settlement_events,
  file.path(
    "Outputs",
    "settlement_events.rds"
  )
)

saveRDS(
  biological_scenarios,
  file.path(
    "Outputs",
    "biological_scenarios.rds"
  )
)

# Save the model benthic surfaces and vertical grid used by the encounter
# calculation so the spatial criteria can be inspected and reproduced later.

saveRDS(
  benthic_surfaces,
  file.path(
    "Outputs",
    "benthic_surfaces.rds"
  )
)

saveRDS(
  vertical_grid,
  file.path(
    "Outputs",
    "llc_vertical_grid.rds"
  )
)

saveRDS(
  particle_trajectories,
  file.path(
    "Outputs",
    "particle_trajectories.rds"
  )
)
