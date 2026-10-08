# Taupo Bank Analysis pt. 2 :)
# Data exploration
#
# Purpose:
# Explore the cleaned project inputs before beginning
# seamount-footprint construction.

library(dplyr)
library(sf)
library(terra)
library(ncdf4)

# Load cleaned project inputs

source("00_setup_and_inputs_.R")

# Node inventory

# Number of nodes in each chain

candidate_nodes %>%
  count(
    chain,
    name = "n_nodes"
  )


# Number of nodes in each feature class

candidate_nodes %>%
  count(
    chain,
    feature_class,
    name = "n_nodes"
  ) %>%
  arrange(
    chain,
    feature_class
  )


# Full node inventory

candidate_nodes %>%
  select(
    node_id,
    name,
    chain,
    feature_class,
    lon,
    lat
  ) %>%
  arrange(
    chain,
    node_id
  ) %>%
  print(n = Inf)


# Plot candidate nodes over GEBCO bathymetry

plot(
  gebco,
  main = "Candidate seamount nodes"
)

plot(
  st_geometry(candidate_nodes_sf),
  add = TRUE,
  pch = 19,
  cex = 0.7
)


# Pilot particle data

# Open Dustin's particle NetCDF

particle_nc <- ncdf4::nc_open(
  particle_file
)


# Number of observations and particles

n_obs <- particle_nc$dim$obs$len
n_particles <- particle_nc$dim$particle$len


# Read particle-level metadata
# Convert NetCDF arrays to ordinary R vectors

merged_id <- as.vector(
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

age_hours <- as.vector(
  ncdf4::ncvar_get(
    particle_nc,
    "age_hours"
  )
)


# Read particle release depth
# Observation 1 is the release position.

release_depth <- as.vector(
  ncdf4::ncvar_get(
    particle_nc,
    "depth",
    start = c(1, 1),
    count = c(1, n_particles)
  )
)


# Read particle release coordinates

release_lon <- as.vector(
  ncdf4::ncvar_get(
    particle_nc,
    "lon",
    start = c(1, 1),
    count = c(1, n_particles)
  )
)

release_lat <- as.vector(
  ncdf4::ncvar_get(
    particle_nc,
    "lat",
    start = c(1, 1),
    count = c(1, n_particles)
  )
)


# Read particle positions at 90 days
# The final observation is observation n_obs.

endpoint_lon <- as.vector(
  ncdf4::ncvar_get(
    particle_nc,
    "lon",
    start = c(n_obs, 1),
    count = c(1, n_particles)
  )
)

endpoint_lat <- as.vector(
  ncdf4::ncvar_get(
    particle_nc,
    "lat",
    start = c(n_obs, 1),
    count = c(1, n_particles)
  )
)


# Close NetCDF after extracting required data

ncdf4::nc_close(
  particle_nc
)

rm(particle_nc)


# Create particle metadata table

particle_metadata <- tibble(
  merged_id = merged_id,
  release_id = release_id,
  
  release_month = case_when(
    release_id == 1 ~ "June",
    release_id == 2 ~ "July",
    release_id == 3 ~ "August",
    TRUE ~ "Unknown"
  ),
  
  above_summit = above_summit,
  
  release_depth_m = round(
    release_depth
  ),
  
  release_lon = release_lon,
  release_lat = release_lat,
  
  endpoint_lon = endpoint_lon,
  endpoint_lat = endpoint_lat
)


# Basic pilot summary

tibble(
  total_particles = n_particles,
  observations_per_particle = n_obs,
  
  maximum_age_days =
    max(
      age_hours,
      na.rm = TRUE
    ) / 24,
  
  retained_particles =
    sum(
      above_summit == 0,
      na.rm = TRUE
    ),
  
  excluded_above_summit =
    sum(
      above_summit == 1,
      na.rm = TRUE
    )
)


# Particle numbers by release month

particle_metadata %>%
  count(
    release_month,
    name = "n_particles"
  )


# Retained particles by release month

particle_metadata %>%
  filter(
    above_summit == 0
  ) %>%
  count(
    release_month,
    name = "retained_particles"
  )


# Release-depth distribution for retained particles

particle_metadata %>%
  filter(
    above_summit == 0
  ) %>%
  count(
    release_depth_m,
    name = "n_particles"
  ) %>%
  arrange(
    release_depth_m
  )


# Number of retained particles within the analysed
# depth range of 200 m and deeper

particle_metadata %>%
  filter(
    above_summit == 0,
    release_depth_m >= 200
  ) %>%
  summarise(
    network_particles = n()
  )


# Spatial extent of particle releases and endpoints

particle_metadata %>%
  summarise(
    release_lon_min = min(
      release_lon,
      na.rm = TRUE
    ),
    
    release_lon_max = max(
      release_lon,
      na.rm = TRUE
    ),
    
    release_lat_min = min(
      release_lat,
      na.rm = TRUE
    ),
    
    release_lat_max = max(
      release_lat,
      na.rm = TRUE
    ),
    
    endpoint_lon_min = min(
      endpoint_lon,
      na.rm = TRUE
    ),
    
    endpoint_lon_max = max(
      endpoint_lon,
      na.rm = TRUE
    ),
    
    endpoint_lat_min = min(
      endpoint_lat,
      na.rm = TRUE
    ),
    
    endpoint_lat_max = max(
      endpoint_lat,
      na.rm = TRUE
    )
  )


# Plot 90-day particle endpoints

valid_endpoints <- particle_metadata %>%
  filter(
    is.finite(endpoint_lon),
    is.finite(endpoint_lat)
  )


plot(
  valid_endpoints$endpoint_lon,
  valid_endpoints$endpoint_lat,
  pch = 16,
  cex = 0.3,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Taupo Bank pilot: 90-day particle endpoints"
)


# Add release locations

points(
  particle_metadata$release_lon,
  particle_metadata$release_lat,
  pch = 17,
  cex = 0.7
)


# Add candidate seamount nodes for spatial context

points(
  candidate_nodes$lon,
  candidate_nodes$lat,
  pch = 1,
  cex = 0.7
)


# Add legend

legend(
  "topright",
  legend = c(
    "90-day endpoints",
    "Release locations",
    "Candidate nodes"
  ),
  pch = c(
    16,
    17,
    1
  ),
  pt.cex = c(
    0.5,
    0.8,
    0.8
  ),
  bty = "n"
)


# Endpoint displacement

# Identify particles with valid release and endpoint coordinates

valid_distance <- (
  is.finite(particle_metadata$release_lon) &
    is.finite(particle_metadata$release_lat) &
    is.finite(particle_metadata$endpoint_lon) &
    is.finite(particle_metadata$endpoint_lat)
)


# Create release-point geometry

release_points <- st_as_sf(
  particle_metadata[
    valid_distance,
  ],
  coords = c(
    "release_lon",
    "release_lat"
  ),
  crs = 4326,
  remove = FALSE
)


# Create endpoint geometry

endpoint_points <- st_as_sf(
  particle_metadata[
    valid_distance,
  ],
  coords = c(
    "endpoint_lon",
    "endpoint_lat"
  ),
  crs = 4326,
  remove = FALSE
)


# Calculate straight-line geodesic distance from
# release position to the 90-day endpoint

endpoint_distance_m <- st_distance(
  st_geometry(release_points),
  st_geometry(endpoint_points),
  by_element = TRUE
)


# Add endpoint displacement to particle metadata

particle_metadata$endpoint_distance_km <- NA_real_

particle_metadata$endpoint_distance_km[
  valid_distance
] <- as.numeric(
  endpoint_distance_m
) / 1000


# Summarise 90-day endpoint displacement
# by individual release depth

depth_displacement_summary <- particle_metadata %>%
  filter(
    above_summit == 0
  ) %>%
  group_by(
    release_depth_m
  ) %>%
  summarise(
    n_particles = n(),
    
    mean_distance_km = mean(
      endpoint_distance_km,
      na.rm = TRUE
    ),
    
    median_distance_km = median(
      endpoint_distance_km,
      na.rm = TRUE
    ),
    
    max_distance_km = max(
      endpoint_distance_km,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    release_depth_m
  )

depth_displacement_summary


# Preliminary pattern:
# 90-day endpoint displacement generally decreases with release depth.
# Shallower releases show substantially greater horizontal displacement
# than deeper releases.
#
# Endpoint displacement represents straight-line distance between the
# release position and final position. It does not represent total
# distance travelled, settlement success, or connectivity.


# Assign particles to analysis depth strata
#
# Releases shallower than 200 m are excluded from the
# depth-resolved network analysis.

particle_metadata <- particle_metadata %>%
  mutate(
    source_stratum = case_when(
      
      release_depth_m >= 200 &
        release_depth_m < 700 ~
        "200-700 m",
      
      release_depth_m >= 700 &
        release_depth_m < 1500 ~
        "700-1500 m",
      
      release_depth_m >= 1500 &
        release_depth_m < 2500 ~
        "1500-2500 m",
      
      release_depth_m >= 2500 ~
        ">2500 m",
      
      TRUE ~ NA_character_
    )
  )


# Put depth strata into shallow-to-deep order

particle_metadata <- particle_metadata %>%
  mutate(
    source_stratum = factor(
      source_stratum,
      levels = c(
        "200-700 m",
        "700-1500 m",
        "1500-2500 m",
        ">2500 m"
      )
    )
  )


# Summarise endpoint displacement by analysis depth stratum

stratum_displacement_summary <- particle_metadata %>%
  filter(
    above_summit == 0,
    !is.na(source_stratum)
  ) %>%
  group_by(
    source_stratum
  ) %>%
  summarise(
    n_particles = n(),
    
    mean_distance_km = mean(
      endpoint_distance_km,
      na.rm = TRUE
    ),
    
    median_distance_km = median(
      endpoint_distance_km,
      na.rm = TRUE
    ),
    
    max_distance_km = max(
      endpoint_distance_km,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

stratum_displacement_summary


# Distribution of endpoint displacement
# by release-depth stratum

boxplot(
  endpoint_distance_km ~ source_stratum,
  
  data = particle_metadata %>%
    filter(
      above_summit == 0,
      !is.na(source_stratum)
    ),
  
  xlab = "Release-depth stratum",
  ylab = "90-day endpoint displacement (km)",
  
  main = "Particle displacement by release depth"
)

