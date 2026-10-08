# Taupo Bank Analysis pt. 4 :)
# Depth-resolved settlement 
#
# Settlement is defined as the first encounter between a competent
# particle and a seamount benthic surface within the settlement distance.
#
# For each competent particle position, the minimum 3-D Euclidean
# distance to the benthic surface is calculated:
# D = sqrt(horizontal_distance^2 + vertical_distance^2)
#
# Horizontal distance is measured from the particle to the GEBCO
# benthic surface patch, while vertical distance is the difference
# between particle depth and seafloor depth at that location.
#
# If more than one seamount satisfies the settlement criterion at the
# same time, the closest benthic surface is used. The seafloor depth of
# this closest surface defines the destination depth stratum.
#
# The settlement distance is to be based on the coarsest vertical resolution
# of the LLC4320 ocean model, so settlement is not resolved at a finer
# scale than the circulation model itself.
#
# However I don't yet have the LLC4320 vertical-grid info.
# I used a settlement distance of 500 m for now as a placeholder since Dustin 
# used this value in his analysis, and mentioned the resolution at the bottom 
# layer was about 500 m

library(sf)
library(terra)
library(dplyr)
library(tibble)
library(ncdf4)


# Load project inputs

source("00_setup_and_inputs_.R")


# Final seamount footprints

final_footprints <- readRDS(
  "final_seamount_footprints.rds"
)


# Current pilot source
#
# Source seamount is defined by Dustin's release simulation,
# rather than inferred from the NETVOLC footprints.

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
# Pre-competency is set to one-third of each PLD.

biological_scenarios <- tibble(
  pld_days = c(
    30,
    60,
    90
  )
) %>%
  mutate(
    precompetency_days =
      pld_days / 3
  )


# Settlement distance
#
# Replace this with the maximum LLC4320 vertical layer thickness
# once the vertical-grid information is available.

settlement_radius_m <- 500


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
# Each GEBCO raster cell inside the final footprint is treated
# as a small benthic surface patch.

build_benthic_surface <- function(
    footprint,
    gebco,
    node_lon,
    node_lat,
    minimum_depth_m = 200) {
  
  footprint_ll <- st_transform(
    footprint,
    4326
  )
  
  
  # Crop GEBCO around the footprint.
  
  local_gebco <- crop(
    gebco,
    vect(footprint_ll),
    snap = "out"
  )
  
  
  elevation_name <- names(
    local_gebco
  )[1]
  
  
  # Convert GEBCO cells to polygons so horizontal distance
  # is measured to the cell surface rather than its centre.
  
  surface_cells <- as.polygons(
    local_gebco,
    values = TRUE,
    na.rm = TRUE,
    aggregate = FALSE
  ) %>%
    st_as_sf()
  
  
  # Clip cells to the final footprint.
  
  surface_cells <- suppressWarnings(
    st_intersection(
      surface_cells,
      st_union(
        st_geometry(
          footprint_ll
        )
      )
    )
  )
  
  surface_cells <- surface_cells[
    !st_is_empty(surface_cells),
  ]
  
  
  # GEBCO elevation is negative below sea level.
  
  surface_cells$z_m <-
    surface_cells[[elevation_name]]
  
  surface_cells$depth_m <-
    -surface_cells$z_m
  
  
  # Only benthic habitat at 200 m depth or deeper is retained.
  
  surface_cells <- surface_cells %>%
    filter(
      depth_m >= minimum_depth_m
    )
  
  
  surface_cells$depth_stratum <-
    as.character(
      assign_depth_stratum(
        surface_cells$depth_m
      )
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
    candidate_nodes,
    gebco) {
  
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
        
        gebco =
          gebco,
        
        node_lon =
          node$lon,
        
        node_lat =
          node$lat
      )
  }
  
  
  surface_list
}


# Settlement search zones
#
# These buffers are only used as a horizontal pre-screen.
# Actual settlement is determined using the 3-D Euclidean distance.

build_settlement_search_zones <- function(
    final_footprints,
    candidate_nodes,
    settlement_radius_m) {
  
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
      dist = settlement_radius_m
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

# Settlement encounters
#
# Each unique eligible particle position is tested against the
# benthic surfaces once. PLD and competency are applied afterwards.

calculate_position_encounters <- function(
    biological_trajectories,
    surface_list,
    search_zones,
    settlement_radius_m) {
  
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
  
  
  # Convert particle positions to spatial points once.
  
  points_ll <- st_as_sf(
    positions,
    coords = c(
      "lon",
      "lat"
    ),
    crs = 4326,
    remove = FALSE
  )
  
  
  # Horizontal pre-screen using the buffered seamount footprints.
  
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
    
    
    # Transform all candidate particle positions for this seamount
    # into the same projected CRS as its benthic surface.
    
    points_m <- points_ll[
      position_rows,
    ] %>%
      st_transform(
        st_crs(
          surface_cells
        )
      )
    
    
    # Find GEBCO patches within the horizontal settlement radius.
    
    nearby_cells <- st_is_within_distance(
      points_m,
      surface_cells,
      dist = settlement_radius_m
    )
    
    
    # Expand the spatial matches into particle-surface pairs.
    
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
    
    
    # Calculate horizontal distance for all candidate pairs together.
    
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
    
    
    surface_z_m <-
      surface_cells$z_m[
        pair_table$cell_index
      ]
    
    
    vertical_distance_m <- abs(
      -particle_depth_m -
        surface_z_m
    )
    
    
    # Full 3-D Euclidean distance.
    
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
        
        distance_m =
          distance_3d_m,
        
        settlement_depth_m =
          surface_cells$depth_m[
            cell_index
          ],
        
        destination_depth_stratum =
          as.character(
            surface_cells$depth_stratum[
              cell_index
            ]
          ),
        
        surface_cell_id =
          surface_cells$surface_cell_id[
            cell_index
          ]
      ) %>%
      filter(
        distance_m <=
          settlement_radius_m
      )
    
    
    if (nrow(pair_table) == 0) {
      next
    }
    
    
    # Retain the closest benthic patch for each particle position
    # on this seamount.
    
    encounter_list[[n]] <- pair_table %>%
      arrange(
        position_id,
        distance_m,
        surface_cell_id
      ) %>%
      group_by(
        position_id
      ) %>%
      slice_head(
        n = 1
      ) %>%
      ungroup() %>%
      select(
        position_id,
        destination_node,
        settlement_depth_m,
        destination_depth_stratum,
        distance_m,
        surface_cell_id
      )
  }
  
  
  encounters <- bind_rows(
    encounter_list
  )
  
  
  if (nrow(encounters) == 0) {
    return(
      tibble()
    )
  }
  
  
  # If more than one seamount qualifies at the same timestep,
  # use the nearest benthic surface.
  
  encounters %>%
    arrange(
      position_id,
      distance_m,
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


# Empty settlement table

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
    distance_m = numeric(),
    settlement_radius_m = numeric(),
    surface_cell_id = integer()
  )
}


# First-settlement rule
#
# Geometric encounters are calculated once for each particle position.
# PLD and competency are then used to identify the first valid settlement
# under each biological scenario.

calculate_settlement_events <- function(
    biological_trajectories,
    surface_list,
    search_zones,
    settlement_radius_m) {
  
  position_encounters <- calculate_position_encounters(
    biological_trajectories =
      biological_trajectories,
    
    surface_list =
      surface_list,
    
    search_zones =
      search_zones,
    
    settlement_radius_m =
      settlement_radius_m
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
      pld_days,
      particle_age_days,
      observation_index
    ) %>%
    group_by(
      particle_id,
      pld_days
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
      distance_m,
      
      settlement_radius_m =
        settlement_radius_m,
      
      surface_cell_id
    )
}


# Build benthic surfaces

benthic_surfaces <- build_all_benthic_surfaces(
  final_footprints =
    final_footprints,
  
  candidate_nodes =
    candidate_nodes,
  
  gebco =
    gebco
)


# Build horizontal search zones

settlement_search_zones <-
  build_settlement_search_zones(
    final_footprints =
      final_footprints,
    
    candidate_nodes =
      candidate_nodes,
    
    settlement_radius_m =
      settlement_radius_m
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


# Find first settlement events

settlement_events <- calculate_settlement_events(
  biological_trajectories =
    biological_trajectories,
  
  surface_list =
    benthic_surfaces,
  
  search_zones =
    settlement_search_zones,
  
  settlement_radius_m =
    settlement_radius_m
)


# Save outputs used by the network analysis

saveRDS(
  network_nodes,
  "network_nodes.rds"
)

saveRDS(
  release_counts,
  "release_counts.rds"
)

saveRDS(
  settlement_events,
  "settlement_events.rds"
)

saveRDS(
  biological_scenarios,
  "biological_scenarios.rds"
)

saveRDS(
  settlement_radius_m,
  "settlement_radius_m.rds"
)


# Quick look at the raw settlement results

settlement_events %>%
  count(
    release_id,
    pld_days,
    source_depth_stratum,
    destination_node,
    destination_depth_stratum,
    name = "n_settled"
  ) %>%
  print(n = Inf)

release_counts %>%
  arrange(
    release_id,
    pld_days,
    source_depth_stratum
  )

# Save trajectories 

saveRDS(
  particle_trajectories,
  "particle_trajectories.rds"
)
