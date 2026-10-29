# 02 Seamount basal-footprint delineation with NETVOLC
#
# Purpose:
# Delineate a horizontal basal footprint for each study feature from GEBCO
# bathymetry using an adaptation of the NETVOLC method of Euillades et al. (2013).
#
# Method overview
# 1. Fit a local quadratic surface to GEBCO bathymetry to estimate slope,
#    aspect and profile convexity.
# 2. Identify concave, outward-facing cells that can represent a basal break.
# 3. Construct a spatial graph and use the iterative NETVOLC minimum-cost path
#    search to obtain a continuous boundary around each feature.
# 4. Repeat delineation with 0.75 degree and 1.0 degree search windows and use
#    numerical overlap plus visual bathymetric checks to select final footprints.
# 5. Treat the Britannia complex separately because the three predefined nodes
#    share one NETVOLC outer boundary; the shared footprint is partitioned with
#    nearest-node Voronoi regions.
#
# The final sf polygons are used only as horizontal seamount footprints. Benthic
# depth for habitat-encounter calculations is subsequently taken from LLC4320v2.
#
# Method origin:
# NETVOLC was developed by Euillades et al. (2013) to delineate
# the basal boundaries of volcanic edifices from digital elevation models.
# Article:
# https://www.sciencedirect.com/science/article/abs/pii/S0098300413000757?via%3Dihub
#
# Provides a more rigorous way to estimate the basal perimeter, rather
# than using a fixed depth contour or simple geometric buffer.
# Identifies terrain consistent with a concave break in slope
# around the volcanic edifice, then uses a network-based minimum-cost
# path to identify a continuous basal boundary.
#
# Adaptations:
# - GEBCO bathymetry is used instead of a terrestrial DEM.
# - Applied to submarine seamounts and raised features instead of land volcanoes
# - Search-window convergence is checked using 0.75 degree
#   and 1.0 degree GEBCO crops.
# - Final footprints are visually inspected against GEBCO to check
#   that the boundary follows the intended basal break rather than
#   an internal slope break, neighbouring relief, or the crop edge.
#
# The boundary is then converted into a 2-D sf polygon
# representing the horizontal basal footprint of the seamount.

library(sf)
library(terra)
library(deldir)
library(igraph)
library(dplyr)
library(ggplot2)

# Load project inputs

candidate_nodes <- readRDS(
  "Outputs/candidate_nodes.rds"
)

candidate_nodes_sf <- readRDS(
  "Outputs/candidate_nodes_sf.rds"
)

gebco_file <- file.path(
  "Data",
  "GEBCO_23_Sep_2026_4a55cbf62670",
  "gebco_2026_n-16.0_s-38.0_w154.0_e161.0.nc"
)

if (!file.exists(gebco_file)) {
  stop("GEBCO file not found: ", gebco_file)
}

gebco <- terra::rast(
  gebco_file,
  subds = "elevation"
)
# Convert longitude/latitude coordinates from degrees
# into a local coordinate system measured in metres.
#
# UTM zone:
# zone = floor((longitude + 180) / 6) + 1

get_utm_epsg <- function(lon, lat) {

  zone <- floor((lon + 180) / 6) + 1

  if (lat < 0) {
    return(32700 + zone)
  } else {
    return(32600 + zone)
  }
}

# Local quadratic terrain fit
# z = ax^2 + by^2 + cxy + dx + ey + f
#
# A 7 x 7 cell neighbourhood is used around each GEBCO cell.
# This reduces whole-metre GEBCO quantisation artefacts and gives
# a smoother estimate of local terrain shape.
# 3 x 3 and 5 x 5 windows were initially tested but were too noisy.

make_quadratic_fit <- function(window_size = 7, dx, dy) {

  # The window must have one centre cell.
  if (window_size %% 2 == 0) {
    stop("window_size must be odd.")
  }

  half_window <- (window_size - 1) / 2

  # Position of each cell relative to the centre, in metres.
  x <- seq(-half_window, half_window) * dx
  y <- seq(half_window, -half_window) * dy

  coords <- expand.grid(
    x = x,
    y = y
  )

  # Matrix form of the quadratic terrain equation.
  # Each row represents one cell in the 7 x 7 neighbourhood.
  X <- cbind(
    x2 = coords$x^2,
    y2 = coords$y^2,
    xy = coords$x * coords$y,
    x = coords$x,
    y = coords$y,
    intercept = 1
  )

  # Pre-calculate the weights used to estimate the
  # quadratic coefficients from the 49 elevation values.
  fit_matrix <- qr.solve(
    X,
    diag(nrow(X))
  )

  return(fit_matrix)
}

# Iterative NETVOLC boundary search
#
# Find the minimum-cost path from start to finish.
# Unusually costly edges in that path are removed, then the
# minimum-cost path is calculated again.
#
# Euillades et al. (2013) found 10 iterations were sufficient.
# Edges with costs greater than mean + 2 SD are removed
# before the next iteration.

run_boundary_search <- function(
    real_edges,
    start_edges,
    finish_edges,
    n_iterations = 10,
    sd_factor = 2) {

  cancelled_edges <- character(0)

  solutions <- vector(
    "list",
    n_iterations
  )

  solution_summary <- data.frame()

  for (iteration in seq_len(n_iterations)) {

    # Give each edge the same ID regardless of direction.
    edge_keys <- paste(
      pmin(real_edges$from, real_edges$to),
      pmax(real_edges$from, real_edges$to),
      sep = "_"
    )

    # Record which edges had already been removed for this iteration.
    cancelled_before <- cancelled_edges

    # Remove edges cancelled in previous iterations.
    active_edges <- real_edges[
      !edge_keys %in% cancelled_edges,
    ]

    # Add the temporary START and FINISH connections.
    all_edges <- rbind(
      active_edges,
      start_edges,
      finish_edges
    )

    graph_current <- graph_from_data_frame(
      all_edges,
      directed = FALSE
    )

    # Find the minimum-cost path across the unwrapped seamount.
    path_current <- shortest_paths(
      graph_current,
      from = "START",
      to = "FINISH",
      weights = E(graph_current)$edge_cost,
      output = "both"
    )

    if (length(path_current$vpath[[1]]) == 0) {
      stop("No complete boundary path could be found.")
    }

    # Extract the ordered candidate cells along the path.
    path_nodes <- V(graph_current)$name[
      as.integer(path_current$vpath[[1]])
    ]

    # Remove the artificial START and FINISH nodes.
    real_nodes <- path_nodes[
      !path_nodes %in% c("START", "FINISH")
    ]

    # Get the edges and costs used by the current path.
    path_edge_ids <- path_current$epath[[1]]

    path_edge_costs <- E(graph_current)[
      path_edge_ids
    ]$edge_cost

    path_edge_ends <- ends(
      graph_current,
      path_edge_ids,
      names = TRUE
    )

    # Exclude the artificial start and finish connections.
    is_real_edge <- !(
      path_edge_ends[, 1] %in% c("START", "FINISH") |
        path_edge_ends[, 2] %in% c("START", "FINISH")
    )

    real_path_costs <- path_edge_costs[is_real_edge]

    # NETVOLC Eq. 6
    #
    # Normalise the total path cost by the number of
    # candidate cells making up the boundary.

    total_cost <- sum(path_edge_costs)

    normalised_cost <- (
      total_cost / length(real_nodes)
    )

    # Remove unusually costly edges before the next iteration.
    #
    # NETVOLC uses:
    # threshold = mean edge cost + 2 SD

    cancellation_threshold <- (
      mean(real_path_costs) +
        sd_factor * sd(real_path_costs)
    )

    cancel_current <- (
      is_real_edge &
        path_edge_costs > cancellation_threshold
    )

    # Give cancelled edges the same ID regardless of direction.
    new_cancelled_edges <- paste(
      pmin(
        path_edge_ends[cancel_current, 1],
        path_edge_ends[cancel_current, 2]
      ),
      pmax(
        path_edge_ends[cancel_current, 1],
        path_edge_ends[cancel_current, 2]
      ),
      sep = "_"
    )

    # Save the current candidate boundary.

    solutions[[iteration]] <- list(
      real_nodes = real_nodes,
      total_cost = total_cost,
      normalised_cost = normalised_cost,
      cancelled_before = cancelled_before
    )

    # Store the main values so the iterations can be compared later.
    solution_summary <- rbind(
      solution_summary,
      data.frame(
        iteration = iteration,
        n_nodes = length(real_nodes),
        total_cost = total_cost,
        normalised_cost = normalised_cost
      )
    )

    # Carry removed edges into the next iteration.
    cancelled_edges <- unique(
      c(
        cancelled_edges,
        new_cancelled_edges
      )
    )
  }

  # Choose the candidate boundary with the lowest normalised cost.

  best_iteration <- which.min(
    solution_summary$normalised_cost
  )

  best_solution <- solutions[[best_iteration]]

  # Return the selected boundary and the iteration summary.
  return(
    list(
      best_iteration = best_iteration,
      best_solution = best_solution,
      summary = solution_summary
    )
  )
}

# Delineate a single seamount

delineate_seamount <- function(seamount, gebco, crop_buffer_deg) {

  # Process one seamount at a time.
  if (nrow(seamount) != 1) {
    stop("seamount must contain exactly one row.")
  }

  # Extract the longitude and latitude of the seamount centre.
  coordinates <- st_coordinates(seamount)

  lon <- coordinates[1, "X"]
  lat <- coordinates[1, "Y"]

  # Find the local UTM coordinate system.
  utm_epsg <- get_utm_epsg(lon, lat)

  # Convert the seamount centre from degrees to metres.
  seamount_utm <- st_transform(
    seamount,
    crs = utm_epsg
  )

  # Crop GEBCO around the seamount centre.
  # crop_buffer_deg = distance from the centre to each edge.
  # 1.0 gives a 2 x 2 degree window; 0.75 gives 1.5 x 1.5 degrees.

  gebco_crop <- crop(
    gebco,
    ext(
      lon - crop_buffer_deg,
      lon + crop_buffer_deg,
      lat - crop_buffer_deg,
      lat + crop_buffer_deg
    )
  )

  if (ncell(gebco_crop) == 0) {
    stop("GEBCO crop contains no cells.")
  }

  # GEBCO cell size
  #
  # Convert the raster resolution from degrees to metres.
  # These distances are needed for the local terrain fit.

  cell_res <- res(gebco_crop)

  centre_point <- st_sfc(
    st_point(c(lon, lat)),
    crs = 4326
  )

  east_point <- st_sfc(
    st_point(c(lon + cell_res[1], lat)),
    crs = 4326
  )

  north_point <- st_sfc(
    st_point(c(lon, lat + cell_res[2])),
    crs = 4326
  )

  dx <- as.numeric(
    st_distance(centre_point, east_point)
  )

  dy <- as.numeric(
    st_distance(centre_point, north_point)
  )

  # dx = east-west cell width (m)
  # dy = north-south cell height (m)

  # Build the 7 x 7 quadratic terrain fit.

  fit_matrix <- make_quadratic_fit(
    window_size = 7,
    dx = dx,
    dy = dy
  )

  # Convert each coefficient row into a 7 x 7 kernel.
  # z = ax^2 + by^2 + cxy + dx + ey + f

  kernel_a <- matrix(fit_matrix[1, ], 7, 7, byrow = TRUE)
  kernel_b <- matrix(fit_matrix[2, ], 7, 7, byrow = TRUE)
  kernel_c <- matrix(fit_matrix[3, ], 7, 7, byrow = TRUE)
  kernel_d <- matrix(fit_matrix[4, ], 7, 7, byrow = TRUE)
  kernel_e <- matrix(fit_matrix[5, ], 7, 7, byrow = TRUE)

  # Apply the quadratic fit across the GEBCO crop.
  # This estimates coefficients a-e for each raster cell.

  coef_a <- focal(
    gebco_crop,
    w = kernel_a,
    fun = "sum",
    na.rm = FALSE
  )

  coef_b <- focal(
    gebco_crop,
    w = kernel_b,
    fun = "sum",
    na.rm = FALSE
  )

  coef_c <- focal(
    gebco_crop,
    w = kernel_c,
    fun = "sum",
    na.rm = FALSE
  )

  coef_d <- focal(
    gebco_crop,
    w = kernel_d,
    fun = "sum",
    na.rm = FALSE
  )

  coef_e <- focal(
    gebco_crop,
    w = kernel_e,
    fun = "sum",
    na.rm = FALSE
  )

  # Use the fitted quadratic coefficients to calculate
  # slope, aspect and profile convexity for each GEBCO cell.

  gradient_squared <- coef_d^2 + coef_e^2

  # Slope in degrees.
  slope <- atan(
    sqrt(gradient_squared)
  ) * 180 / pi

  # Downhill aspect in degrees.
  # 0 = north, 90 = east, 180 = south, 270 = west.
  aspect <- (
    atan2(
      -coef_d,
      -coef_e
    ) * 180 / pi + 360
  ) %% 360

  # Aspect is undefined on perfectly flat terrain.
  aspect[gradient_squared == 0] <- NA

  # Profile convexity.
  # Negative values = concave terrain.
  # Positive values = convex terrain.
  profile_convexity <- (
    -200 * (
      coef_a * coef_d^2 +
        coef_b * coef_e^2 +
        coef_c * coef_d * coef_e
    )
  ) / (
    gradient_squared *
      (1 + gradient_squared)^(3 / 2)
  )

  # Convexity is undefined on perfectly flat terrain.
  profile_convexity[gradient_squared == 0] <- NA

  # Candidate basal-break cells
  #
  # A candidate must:
  # 1. be concave (profile_convexity < 0), and
  # 2. face broadly away from the seamount centre.

  # Get the centre coordinates of every GEBCO cell.
  cell_numbers <- seq_len(ncell(aspect))
  cell_xy <- xyFromCell(aspect, cell_numbers)

  # Convert degree offsets from the seamount centre into metres.
  metres_per_degree_lon <- dx / cell_res[1]
  metres_per_degree_lat <- dy / cell_res[2]

  radial_dx <- (
    cell_xy[, 1] - lon
  ) * metres_per_degree_lon

  radial_dy <- (
    cell_xy[, 2] - lat
  ) * metres_per_degree_lat

  # Direction from the seamount centre to each raster cell.
  outward_aspect_values <- (
    atan2(
      radial_dx,
      radial_dy
    ) * 180 / pi + 360
  ) %% 360

  outward_aspect <- setValues(
    rast(aspect),
    outward_aspect_values
  )

  # Difference between the actual slope direction and the direction pointing
  # away from seamount.
  aspect_difference <- abs(
    aspect - outward_aspect
  )

  aspect_difference <- ifel(
    aspect_difference > 180,
    360 - aspect_difference,
    aspect_difference
  )

  # Keep concave cells facing within 90 degrees of outward.
  candidate_mask <- (
    profile_convexity < 0 &
      aspect_difference < 90
  )

  # Candidate cost - NETVOLC Eq. 4c
  #
  # Cost = -1 / profile convexity + slope
  #
  # Eq. 4c (from link at top) includes slope to reduce selection of steep
  # internal slope breaks and favour lower-gradient breaks near the base.
  # Only concave, outward-facing candidate cells receive a cost.

  candidate_cost <- ifel(
    candidate_mask,
    (-1 / profile_convexity) + slope,
    NA
  )

  # Candidate cell table
  #
  # Keep only valid candidate cells and store the terrain
  # information needed to build the NETVOLC network.

  candidate_values <- values(candidate_mask, mat = FALSE)
  cost_values <- values(candidate_cost, mat = FALSE)
  convexity_values <- values(profile_convexity, mat = FALSE)
  slope_values <- values(slope, mat = FALSE)

  candidate_cells <- which(
    candidate_values == 1 &
      is.finite(cost_values)
  )

  # At least three points are needed to build a spatial network.
  if (length(candidate_cells) < 3) {
    stop("Too few valid NETVOLC candidate cells.")
  }

  # Get the longitude and latitude of each candidate cell.
  candidate_xy <- xyFromCell(
    gebco_crop,
    candidate_cells
  )

  candidates <- data.frame(
    cell = candidate_cells,
    lon = candidate_xy[, 1],
    lat = candidate_xy[, 2],
    profile_convexity = convexity_values[candidate_cells],
    slope_deg = slope_values[candidate_cells],
    cost = cost_values[candidate_cells]
  )

  # Candidate positions relative to the seamount centre
  #
  # Convert candidate locations to metres, then describe each one by:
  # radius = distance from the seamount centre
  # theta  = angle around the seamount centre

  candidates_sf <- st_as_sf(
    candidates,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )

  candidates_utm <- st_transform(
    candidates_sf,
    crs = utm_epsg
  )

  candidate_xy_m <- st_coordinates(candidates_utm)

  candidates$x_m <- candidate_xy_m[, "X"]
  candidates$y_m <- candidate_xy_m[, "Y"]

  # Seamount centre in metres.
  centre_m <- st_coordinates(seamount_utm)[1, ]

  # Position of each candidate relative to the centre.
  candidate_dx <- candidates$x_m - centre_m["X"]
  candidate_dy <- candidates$y_m - centre_m["Y"]

  # Radial distance from the seamount centre.
  candidates$radius_m <- sqrt(
    candidate_dx^2 + candidate_dy^2
  )

  # Angle around the seamount centre.
  candidates$theta_deg <- atan2(
    candidate_dy,
    candidate_dx
  ) * 180 / pi

  # Express radial distance in approximate GEBCO-cell units.
  cell_size_m <- sqrt(dx * dy)

  candidates$radius_px <- (
    candidates$radius_m / cell_size_m
  )

  # Build the candidate network
  #
  # Connect nearby candidate cells using Delaunay triangulation
  # based on angle around the seamount and distance from the centre.

  delaunay <- deldir(
    x = candidates$theta_deg,
    y = candidates$radius_px,
    id = seq_len(nrow(candidates))
  )

  # Extract the connections between candidate cells.
  edges <- delaunay$delsgs

  # Candidate IDs are returned as characters, so convert to integers.
  edges$ind1 <- as.integer(edges$ind1)
  edges$ind2 <- as.integer(edges$ind2)

  # NETVOLC edge cost - Eq. 5
  #
  # Edge cost = ((C1 + C2) / 2) * L^2 , where:
  # C1, C2 = costs of the two candidate cells
  # L = distance between them in polar coordinates
  # Lower-cost edges are preferred by the boundary search.

  # Attach the terrain cost at each end of the edge.
  edges$cost_1 <- candidates$cost[edges$ind1]
  edges$cost_2 <- candidates$cost[edges$ind2]

  # Calculate edge length using polar coordinates.
  edges$length <- sqrt(
    (edges$x2 - edges$x1)^2 +
      (edges$y2 - edges$y1)^2
  )

  # NETVOLC Eq. 5.
  edges$edge_cost <- (
    (edges$cost_1 + edges$cost_2) / 2
  ) * edges$length^2

  # Filter the candidate network
  #
  # NETVOLC removes unusually long edges and convex-hull edges
  # as they can create unrealistic shortcuts in the network.

  # Remove the longest 1% of edges.
  # Chose 99th percentile to define unusually long edges

  length_cutoff <- quantile(
    edges$length,
    0.99,
    na.rm = TRUE
  )

  edges <- edges[
    edges$length <= length_cutoff,
  ]

  # Remove edges forming the convex hull of the network.
  hull_nodes <- chull(
    candidates$theta_deg,
    candidates$radius_px
  )

  hull_pairs <- cbind(
    hull_nodes,
    c(hull_nodes[-1], hull_nodes[1])
  )

  # Give each edge a direction-independent ID.
  hull_keys <- paste(
    pmin(hull_pairs[, 1], hull_pairs[, 2]),
    pmax(hull_pairs[, 1], hull_pairs[, 2]),
    sep = "_"
  )

  edge_keys <- paste(
    pmin(edges$ind1, edges$ind2),
    pmax(edges$ind1, edges$ind2),
    sep = "_"
  )

  edges <- edges[
    !edge_keys %in% hull_keys,
  ]

  # Prepare the graph for the boundary search
  #
  # Polar seam:
  # Converting candidate positions to angles unwraps the seamount
  # from -180 to +180 degrees.
  # -180 and +180 represent the same direction,
  # but appear at opposite ends of the polar coordinates.
  #
  # A 5 degree band was used on either side of the polar seam,
  # following the NETVOLC v1.1 implementation.
  #
  # Since an exact candidate may not occur at +/-180 degrees, the
  # small band allows nearby candidates to act as temporary
  # seam endpoints (so the path forms a continuous closed polygon).

  # Build the graph from the filtered Delaunay edges.
  graph_base <- graph_from_data_frame(
    data.frame(
      from = as.character(edges$ind1),
      to = as.character(edges$ind2),
      edge_cost = edges$edge_cost
    ),
    directed = FALSE,
    vertices = data.frame(
      name = as.character(seq_len(nrow(candidates)))
    )
  )

  # Identify connected parts of the candidate network.
  graph_components <- components(graph_base)

  candidates$component <- as.integer(
    graph_components$membership
  )

  # Candidate cells near the -180 degree side.
  start_candidates <- which(
    candidates$theta_deg <= -175
  )

  # Candidate cells near the +180 degree side.
  finish_candidates <- which(
    candidates$theta_deg >= 175
  )

  # Find components that reach both sides of the polar seam.
  eligible_components <- intersect(
    unique(candidates$component[start_candidates]),
    unique(candidates$component[finish_candidates])
  )

  # A complete boundary requires a network spanning both sides.
  if (length(eligible_components) == 0) {
    stop("No candidate network spans the polar seam.")
  }

  # If more than one component spans the seam,
  # use the largest connected component.
  component_sizes <- table(candidates$component)

  main_component <- eligible_components[
    which.max(
      component_sizes[
        as.character(eligible_components)
      ]
    )
  ]

  # Keep only seam candidates from the selected component.
  start_nodes <- start_candidates[
    candidates$component[start_candidates] == main_component
  ]

  finish_nodes <- finish_candidates[
    candidates$component[finish_candidates] == main_component
  ]

  # Temporary start and finish nodes
  #
  # 'START' and 'FINISH' are artificial nodes used only for the
  # minimum-cost path calculation.
  # They connect to the candidate cells near each side of the seam
  # with zero cost, so they won't affect the overall terrain cost.

  start_edges <- data.frame(
    from = "START",
    to = as.character(start_nodes),
    edge_cost = 0
  )

  finish_edges <- data.frame(
    from = as.character(finish_nodes),
    to = "FINISH",
    edge_cost = 0
  )

  # Store the real candidate-network edges separately.
  real_edges <- data.frame(
    from = as.character(edges$ind1),
    to = as.character(edges$ind2),
    edge_cost = edges$edge_cost
  )

  # Run the iterative NETVOLC boundary search (F1).

  boundary_search <- run_boundary_search(
    real_edges = real_edges,
    start_edges = start_edges,
    finish_edges = finish_edges
  )

  # Candidate cells making up the lowest normalised-cost boundary.
  best_path <- boundary_search$best_solution$real_nodes

  # Stitch the polar seam
  #
  # The selected path currently runs from one side of the
  # -180 / +180 degree seam to the other.
  # Rebuild the network from the selected iteration so the
  # two ends can be joined into one continuous boundary.

  # Edges removed before the selected iteration.
  cancelled_before <- boundary_search$best_solution$cancelled_before

  # Give each real edge the same ID regardless of direction.
  real_edge_keys <- paste(
    pmin(real_edges$from, real_edges$to),
    pmax(real_edges$from, real_edges$to),
    sep = "_"
  )

  # Rebuild the network that produced the selected boundary.
  final_edges <- real_edges[
    !real_edge_keys %in% cancelled_before,
  ]

  final_graph <- graph_from_data_frame(
    final_edges,
    directed = FALSE
  )

  # First and last candidate cells in the selected path.
  start_node <- as.integer(best_path[1])
  finish_node <- as.integer(best_path[length(best_path)])

  # Exact seam points
  #
  # The seam lies at +/-180 degrees, directly west of the
  # seamount centre in the local polar-coordinate system.
  # Two possible seam radii are tested: one based on each
  # end of the selected boundary path.

  start_radius_m <- candidates$radius_m[start_node]
  finish_radius_m <- candidates$radius_m[finish_node]

  seam_points_m <- data.frame(
    solution = c("start_radius", "finish_radius"),
    radius_m = c(start_radius_m, finish_radius_m),
    x_m = c(
      centre_m["X"] - start_radius_m,
      centre_m["X"] - finish_radius_m
    ),
    y_m = c(
      centre_m["Y"],
      centre_m["Y"]
    )
  )

  # Convert the seam points back to longitude and latitude.

  seam_points_sf <- st_as_sf(
    seam_points_m,
    coords = c("x_m", "y_m"),
    crs = utm_epsg,
    remove = FALSE
  )

  seam_points_ll <- st_transform(
    seam_points_sf,
    crs = 4326
  )

  seam_lonlat <- st_coordinates(seam_points_ll)

  seam_points_m$lon <- seam_lonlat[, "X"]
  seam_points_m$lat <- seam_lonlat[, "Y"]

  # Find the candidate cells closest to each exact seam point.
  # Only candidates still present in the selected network are used.

  active_nodes <- V(final_graph)$name

  active_start_nodes <- start_nodes[
    as.character(start_nodes) %in% active_nodes
  ]

  active_finish_nodes <- finish_nodes[
    as.character(finish_nodes) %in% active_nodes
  ]

  if (length(active_start_nodes) == 0 || length(active_finish_nodes) == 0) {
    stop("The selected NETVOLC network has no active seam candidates.")
  }

  # Find candidate cells nearest to the exact seam points.
  #
  # The exact seam point usually falls between candidate cells.
  # For each possible seam radius, find the nearest active candidate
  # on the START and FINISH sides.
  #
  # This gives two possible seam closures - keep the lower cost closure

  # Closure option 1: use the radius of the START end.

  start_distances <- sqrt(
    (candidates$x_m[active_start_nodes] - seam_points_m$x_m[1])^2 +
      (candidates$y_m[active_start_nodes] - seam_points_m$y_m[1])^2
  )

  finish_distances <- sqrt(
    (candidates$x_m[active_finish_nodes] - seam_points_m$x_m[1])^2 +
      (candidates$y_m[active_finish_nodes] - seam_points_m$y_m[1])^2
  )

  closure1_start <- active_start_nodes[which.min(start_distances)]
  closure1_finish <- active_finish_nodes[which.min(finish_distances)]

  # Closure option 2: use the radius of the FINISH end.

  start_distances <- sqrt(
    (candidates$x_m[active_start_nodes] - seam_points_m$x_m[2])^2 +
      (candidates$y_m[active_start_nodes] - seam_points_m$y_m[2])^2
  )

  finish_distances <- sqrt(
    (candidates$x_m[active_finish_nodes] - seam_points_m$x_m[2])^2 +
      (candidates$y_m[active_finish_nodes] - seam_points_m$y_m[2])^2
  )

  closure2_start <- active_start_nodes[which.min(start_distances)]
  closure2_finish <- active_finish_nodes[which.min(finish_distances)]

  # Test the two possible seam closures.
  #
  # For each option, find the minimum-cost path between the
  # nearest START- and FINISH-side candidates.

  closure1 <- shortest_paths(
    final_graph,
    from = as.character(closure1_start),
    to = as.character(closure1_finish),
    weights = E(final_graph)$edge_cost,
    output = "both"
  )

  closure2 <- shortest_paths(
    final_graph,
    from = as.character(closure2_start),
    to = as.character(closure2_finish),
    weights = E(final_graph)$edge_cost,
    output = "both"
  )

  # Check whether each closure produced a path.
  closure1_exists <- length(closure1$vpath[[1]]) > 0
  closure2_exists <- length(closure2$vpath[[1]]) > 0

  if (!closure1_exists && !closure2_exists) {
    stop("Neither seam closure could be found.")
  }

  # Calculate the normalised cost of each available closure.
  closure1_cost <- Inf
  closure2_cost <- Inf

  if (closure1_exists) {

    # Extract the candidate cells used by this closure.
    closure1_nodes <- V(final_graph)$name[
      as.integer(closure1$vpath[[1]])
    ]

    closure1_cost <- sum(
      E(final_graph)[closure1$epath[[1]]]$edge_cost
    ) / length(closure1_nodes)
  }

  if (closure2_exists) {

    # Extract the candidate cells used by this closure.
    closure2_nodes <- V(final_graph)$name[
      as.integer(closure2$vpath[[1]])
    ]

    closure2_cost <- sum(
      E(final_graph)[closure2$epath[[1]]]$edge_cost
    ) / length(closure2_nodes)
  }

  # Keep the lower-cost seam closure.

  if (closure1_cost <= closure2_cost) {

    seam_path <- closure1_nodes
    seam_point <- seam_points_m[1, ]
    final_normalised_cost <- closure1_cost

  } else {

    seam_path <- closure2_nodes
    seam_point <- seam_points_m[2, ]
    final_normalised_cost <- closure2_cost
  }

  # Final boundary coordinates
  #
  # seam_path is the final ordered boundary after the
  # polar-seam continuity check.

  final_nodes <- as.integer(seam_path)

  final_xy <- candidates[
    final_nodes,
    c("lon", "lat")
  ]

  # Add the exact seam point, then repeat the first point
  # to close the boundary into a ring.

  ring_xy <- rbind(
    as.matrix(final_xy),
    c(seam_point$lon, seam_point$lat),
    as.numeric(final_xy[1, ])
  )

  # Convert the closed boundary into a 2-D sf polygon.

  footprint_polygon <- st_sfc(
    st_polygon(
      list(ring_xy)
    ),
    crs = 4326
  )

  # Check that the final boundary forms a valid polygon.

  if (!st_is_valid(footprint_polygon)) {
    stop("Final basal-footprint polygon is invalid.")
  }

  # Project the polygon into metres for area and perimeter calculations.

  footprint_utm <- st_transform(
    footprint_polygon,
    crs = utm_epsg
  )

  # Footprint area in km2.
  area_km2 <- as.numeric(
    st_area(footprint_utm)
  ) / 1e6

  # Footprint perimeter in km.
  perimeter_km <- as.numeric(
    st_length(st_boundary(footprint_utm))
  ) / 1000

  # Store the final basal footprint and main method information.

  footprint <- st_sf(
    node_id = seamount$node_id,
    name = seamount$name,
    chain = seamount$chain,
    crop_buffer_deg = crop_buffer_deg,
    area_km2 = area_km2,
    perimeter_km = perimeter_km,
    best_iteration = boundary_search$best_iteration,
    seam_solution = seam_point$solution,
    normalised_cost = final_normalised_cost,
    geometry = footprint_polygon
  )

  # Return the final footprint and useful information for checking in script 04.

  return(
    list(
      footprint = footprint,
      candidates = candidates,
      best_path = best_path,
      seam_path = seam_path,
      boundary_search = boundary_search
    )
  )

}

# Apply NETVOLC to all study features
#
# Apply the single-seamount delineation function to each
# seamount using the same crop size.

delineate_all_seamounts <- function(
    seamounts,
    gebco,
    crop_buffer_deg) {

  results <- vector(
    "list",
    nrow(seamounts)
  )

  for (i in seq_len(nrow(seamounts))) {

    seamount_i <- seamounts[
      i,
      ,
      drop = FALSE
    ]

    # Keep individual failures from stopping the full batch.
    results[[i]] <- tryCatch({
      result_i <- delineate_seamount(
        seamount = seamount_i,
        gebco = gebco,
        crop_buffer_deg = crop_buffer_deg
      )
      # Only keep the final footprint.
      footprint_i <- result_i$footprint
      rm(result_i)
      footprint_i
    }, error = function(e) {
      list(
        node_id = seamount_i$node_id,
        name = seamount_i$name,
        error = conditionMessage(e)
      )
    })

    gc()
  }

  # Identify successful delineations.
  successful <- vapply(
    results,
    function(x) inherits(x, "sf"),
    logical(1)
  )

  # Combine successful footprints.
  if (any(successful)) {

    footprints <- do.call(
      rbind,
      results[successful]
    )

  } else {

    footprints <- NULL
  }

  return(
    list(
      results = results,
      footprints = footprints,
      successful = successful
    )
  )
}

# Delineate all features with a 0.75 degree search radius.
# Every study node is required downstream, so failed delineations are treated as
# an error rather than being silently omitted.

results_075 <- delineate_all_seamounts(
  seamounts = candidate_nodes_sf,
  gebco = gebco,
  crop_buffer_deg = 0.75
)

failed_075 <- candidate_nodes_sf$node_id[
  !results_075$successful
]

if (length(failed_075) > 0) {
  stop(
    "NETVOLC delineation failed for the 0.75 degree search window: ",
    paste(failed_075, collapse = ", ")
  )
}

saveRDS(
  results_075,
  file.path(
    "Outputs",
    "results_075.rds"
  )
)

# Repeat the delineation with a 1.0 degree search radius.
# Comparing two search windows helps identify boundaries that are sensitive to
# the spatial extent supplied to NETVOLC.

results_100 <- delineate_all_seamounts(
  seamounts = candidate_nodes_sf,
  gebco = gebco,
  crop_buffer_deg = 1.0
)

failed_100 <- candidate_nodes_sf$node_id[
  !results_100$successful
]

if (length(failed_100) > 0) {
  stop(
    "NETVOLC delineation failed for the 1.0 degree search window: ",
    paste(failed_100, collapse = ", ")
  )
}

saveRDS(
  results_100,
  file.path(
    "Outputs",
    "results_100.rds"
  )
)

# Compare the two search-window solutions.

# Extract the two complete footprint sets.

footprints_075 <- results_075$footprints
footprints_100 <- results_100$footprints

# Compare footprint area and perimeter before assessing spatial overlap.

comparison <- footprints_075 %>%
  st_drop_geometry() %>%
  select(
    node_id,
    name,
    chain,
    area_075 = area_km2,
    perimeter_075 = perimeter_km
  ) %>%
  left_join(
    footprints_100 %>%
      st_drop_geometry() %>%
      select(
        node_id,
        area_100 = area_km2,
        perimeter_100 = perimeter_km
      ),
    by = "node_id"
  ) %>%
  mutate(
    area_difference_km2 =
      area_100 - area_075,
    area_difference_percent =
      100 * (area_100 - area_075) / area_075,
    perimeter_difference_km =
      perimeter_100 - perimeter_075,
    perimeter_difference_percent =
      100 * (perimeter_100 - perimeter_075) / perimeter_075
  )

# Spatial overlap between the two solutions

overlap <- vector(
  "list",
  nrow(footprints_075)
)

for (i in seq_len(nrow(footprints_075))) {

  node <- footprints_075$node_id[i]

  polygon_075 <- footprints_075[
    footprints_075$node_id == node,
  ]

  polygon_100 <- footprints_100[
    footprints_100$node_id == node,
  ]

  # Choose the local UTM zone from the footprint location.

  bbox <- st_bbox(polygon_075)

  longitude <- mean(
    c(
      bbox["xmin"],
      bbox["xmax"]
    )
  )

  latitude <- mean(
    c(
      bbox["ymin"],
      bbox["ymax"]
    )
  )

  utm_zone <- floor(
    (longitude + 180) / 6
  ) + 1

  if (latitude < 0) {
    utm_epsg <- 32700 + utm_zone
  } else {
    utm_epsg <- 32600 + utm_zone
  }

  # Project to metres and repair any small geometry problems,
  # such as duplicate vertices.

  polygon_075_m <- polygon_075 %>%
    st_transform(utm_epsg) %>%
    st_make_valid()

  polygon_100_m <- polygon_100 %>%
    st_transform(utm_epsg) %>%
    st_make_valid()

  # Calculate intersection and union area.

  intersection <- st_intersection(
    st_geometry(polygon_075_m),
    st_geometry(polygon_100_m)
  )

  area_075 <- sum(
    as.numeric(
      st_area(polygon_075_m)
    )
  )

  area_100 <- sum(
    as.numeric(
      st_area(polygon_100_m)
    )
  )

  intersection_area <- sum(
    as.numeric(
      st_area(intersection)
    )
  )

  union_area <-
    area_075 +
    area_100 -
    intersection_area

  overlap[[i]] <- tibble(
    node_id = node,
    overlap_percent =
      100 * intersection_area / union_area
  )
}

overlap <- bind_rows(overlap)

comparison <- comparison %>%
  left_join(
    overlap,
    by = "node_id"
  )

write.csv(
  comparison,
  file.path(
    "Outputs",
    "netvolc_075_vs_100_comparison.csv"
  ),
  row.names = FALSE
)

# NETVOLC search-window comparison
#
# Compare the 0.75 degree and 1.0 degree footprint solutions.
# The numerical comparison is used to identify seamounts where
# the two search-window sizes produce meaningfully different boundaries.
# Plot 0.75 and 1.0 degree NETVOLC footprints
#
# Both solutions are plotted over GEBCO so differences in
# the estimated basal boundary can be visually checked.

plot_footprint_comparison <- function(node) {

  # Select the two footprint solutions.

  polygon_075 <- footprints_075 %>%
    filter(node_id == node)

  polygon_100 <- footprints_100 %>%
    filter(node_id == node)

  # Find an extent containing both footprints.

  bbox_075 <- st_bbox(polygon_075)
  bbox_100 <- st_bbox(polygon_100)

  comparison_extent <- terra::ext(
    min(bbox_075["xmin"], bbox_100["xmin"]),
    max(bbox_075["xmax"], bbox_100["xmax"]),
    min(bbox_075["ymin"], bbox_100["ymin"]),
    max(bbox_075["ymax"], bbox_100["ymax"])
  )

  # Crop GEBCO to the footprint area.

  gebco_crop <- terra::crop(
    gebco,
    comparison_extent
  )

  # Plot GEBCO bathymetry.

  plot(
    gebco_crop,
    main = paste(
      node,
      polygon_075$name
    )
  )

  # Solid line = 0.75 degree solution.

  plot(
    st_geometry(polygon_075),
    add = TRUE,
    border = "black",
    lwd = 2
  )

  # Dashed line = 1.0 degree solution.

  plot(
    st_geometry(polygon_100),
    add = TRUE,
    border = "black",
    lwd = 2,
    lty = 2
  )

  legend(
    "topright",
    legend = c(
      "0.75 degree",
      "1.0 degree"
    ),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}

# Order seamounts from lowest to highest polygon overlap.
#
# Lower overlap means the two NETVOLC solutions differ more,
# so these sites should be visually checked first.

qc_table <- comparison %>%
  arrange(overlap_percent) %>%
  select(
    node_id,
    name,
    chain,
    area_075,
    area_100,
    area_difference_percent,
    perimeter_difference_percent,
    overlap_percent
  )

# Create a table for recording the final footprint choice.
#
# Each seamount will ultimately use either the 0.75 degree
# or 1.0 degree NETVOLC solution.

# Begin visual quality control with the feature showing the lowest spatial
# agreement between search-window solutions.

lowest_overlap_node <- qc_table$node_id[1]

plot_footprint_comparison(
  lowest_overlap_node
)

# Use the 0.75 degree solution by default.

footprint_choices <- qc_table %>%
  select(
    node_id,
    name,
    chain,
    overlap_percent
  ) %>%
  mutate(
    final_choice = "0.75",
    qc_notes = "Footprints identical or very similar"
  )

# Chesterfield Plateau (L02) uses the 1.0 degree solution because visual
# inspection against GEBCO indicates that it follows the basal break more closely.

footprint_choices$final_choice[
  footprint_choices$node_id == "L02"
] <- "1.0"

footprint_choices$qc_notes[
  footprint_choices$node_id == "L02"
] <- "1.0 degree footprint follows basal boundary better"

# Britannia requires separate treatment because NETVOLC identifies T15, T16
# and T17 as one continuous basal complex. Inspect the shared boundary and the
# three node locations before partitioning the complex.

plot_footprint_comparison("T15")
plot_footprint_comparison("T16")
plot_footprint_comparison("T17")

# Use the T16 0.75 degree NETVOLC solution as the shared outer boundary of the
# Britannia complex, then partition that boundary among the three node centres.

britannia_parent <- footprints_075 %>%
  filter(node_id == "T16")

# Get the three Britannia node centres.

britannia_nodes <- candidate_nodes_sf %>%
  filter(node_id %in% c("T15", "T16", "T17"))

# Work in metres rather than longitude/latitude.
# Britannia is in UTM Zone 56S.

britannia_parent_m <- st_transform(
  britannia_parent,
  32756
)

britannia_nodes_m <- st_transform(
  britannia_nodes,
  32756
)

# Create one geometry containing the three node points.

britannia_points <- st_union(
  st_geometry(britannia_nodes_m)
)

# Make Voronoi regions around the three nodes.
# A large envelope ensures that all three Voronoi regions fully cover the
# parent footprint before clipping.

voronoi <- st_voronoi(
  britannia_points,
  envelope = st_as_sfc(
    st_bbox(
      st_buffer(
        britannia_parent_m,
        10000
      )
    )
  )
)

# Extract the individual polygons.

voronoi_polygons <- st_collection_extract(
  voronoi,
  "POLYGON"
)

voronoi_sf <- st_sf(
  geometry = voronoi_polygons
)

# Assign each Voronoi polygon to the node contained inside it.

voronoi_sf <- st_join(
  voronoi_sf,
  britannia_nodes_m %>%
    select(node_id, name, chain),
  join = st_intersects
)

# Clip the Voronoi regions to the actual NETVOLC parent footprint.

britannia_split_m <- st_intersection(
  voronoi_sf,
  st_geometry(britannia_parent_m)
)

# Calculate area of each resulting node polygon.

britannia_split_m$area_km2 <- as.numeric(
  st_area(britannia_split_m)
) / 1e6

# Convert back to longitude/latitude.

britannia_split <- st_transform(
  britannia_split_m,
  4326
)

# Plot the result.

plot(
  st_geometry(britannia_parent),
  border = "black",
  lwd = 3,
  main = "Britannia complex split among T15-T17"
)

plot(
  st_geometry(britannia_split),
  add = TRUE,
  border = "purple",
  lwd = 2
)

plot(
  st_geometry(britannia_nodes),
  add = TRUE,
  pch = 19
)

text(
  st_coordinates(britannia_nodes),
  labels = britannia_nodes$node_id,
  pos = 3
)

# Plot the Britannia partition over GEBCO bathymetry for visual quality control.

# Get extent of the shared Britannia footprint
britannia_bbox <- st_bbox(britannia_parent)

# Add a small amount of space around it
britannia_extent <- terra::ext(
  britannia_bbox["xmin"] - 0.1,
  britannia_bbox["xmax"] + 0.1,
  britannia_bbox["ymin"] - 0.1,
  britannia_bbox["ymax"] + 0.1
)

# Crop GEBCO to the Britannia area
britannia_gebco <- terra::crop(
  gebco,
  britannia_extent
)

# Plot bathymetry
plot(
  britannia_gebco,
  main = "Britannia complex: GEBCO and node partitions"
)

# Add GEBCO depth contours
contour(
  britannia_gebco,
  add = TRUE,
  drawlabels = TRUE
)

# Add outer NETVOLC footprint
plot(
  st_geometry(britannia_parent),
  add = TRUE,
  border = "black",
  lwd = 3
)

# Add the three Voronoi-split polygons
plot(
  st_geometry(britannia_split),
  add = TRUE,
  border = "purple",
  lwd = 2
)

# Add the three Britannia node centres
plot(
  st_geometry(britannia_nodes),
  add = TRUE,
  pch = 19
)

# Label the nodes
text(
  st_coordinates(britannia_nodes),
  labels = britannia_nodes$node_id,
  pos = 3
)

# Build the final set of seamount footprints.
#
# Most nodes use the 0.75 degree NETVOLC solution.
# Chesterfield Plateau (L02) uses the 1.0 degree solution because
# visual inspection showed that it followed the basal boundary better.
#
# North, Central and South Britannia (T15-T17) are treated separately
# because NETVOLC identified the three guyots as one shared basal complex.
# Their shared outer footprint is therefore partitioned among the three
# predefined node centres using the Voronoi split created above.

# Start with the 0.75 degree NETVOLC footprints for all nodes.

final_footprints <- footprints_075

# Replace Chesterfield Plateau with the preferred 1.0 degree footprint.

final_footprints <- final_footprints %>%
  filter(node_id != "L02")

chesterfield_final <- footprints_100 %>%
  filter(node_id == "L02")

final_footprints <- bind_rows(
  final_footprints,
  chesterfield_final
)

# Remove the original NETVOLC footprints for the three Britannia nodes.
#
# These polygons overlap because NETVOLC identifies the Britannia guyots
# as part of one continuous basal complex. They will be replaced below
# with the three non-overlapping Voronoi-partitioned polygons.

final_footprints <- final_footprints %>%
  filter(
    !node_id %in% c("T15", "T16", "T17")
  )

# Prepare the three Britannia polygons for inclusion in the final dataset.
#
# These polygons share the NETVOLC-derived outer boundary of the
# Britannia complex, but the internal boundaries are nearest-node
# allocation boundaries rather than independently derived basal breaks.

britannia_final <- britannia_split %>%
  mutate(
    crop_buffer_deg = NA_real_,
    perimeter_km = NA_real_,
    best_iteration = NA_integer_,
    seam_solution = "Voronoi partition",
    normalised_cost = NA_real_
  )

# Calculate area and perimeter for each Britannia partition.
#
# UTM Zone 56S is used so that distances and areas are measured in metres.

britannia_final_m <- st_transform(
  britannia_final,
  32756
)

britannia_final$area_km2 <- as.numeric(
  st_area(britannia_final_m)
) / 1e6

britannia_final$perimeter_km <- as.numeric(
  st_length(
    st_boundary(britannia_final_m)
  )
) / 1000

# Keep the same columns used by the standard NETVOLC footprints.
#
# This allows the Britannia polygons to be combined with the other nodes
# without changing the structure of the final spatial dataset.

britannia_final <- britannia_final %>%
  select(
    node_id,
    name,
    chain,
    crop_buffer_deg,
    area_km2,
    perimeter_km,
    best_iteration,
    seam_solution,
    normalised_cost,
    geometry
  )

# Add the three Britannia partitions back into the final footprint dataset.

final_footprints <- bind_rows(
  final_footprints,
  britannia_final
)

# Validate the final footprint inventory.
# Each study node must have exactly one polygon before the footprints are used
# to construct LLC4320v2 benthic habitat surfaces.

expected_node_ids <- sort(candidate_nodes_sf$node_id)
final_node_ids <- sort(final_footprints$node_id)

missing_final_nodes <- setdiff(
  expected_node_ids,
  final_node_ids
)

unexpected_final_nodes <- setdiff(
  final_node_ids,
  expected_node_ids
)

if (length(missing_final_nodes) > 0) {
  stop(
    "Final footprint dataset is missing node(s): ",
    paste(missing_final_nodes, collapse = ", ")
  )
}

if (length(unexpected_final_nodes) > 0) {
  stop(
    "Final footprint dataset contains unexpected node(s): ",
    paste(unexpected_final_nodes, collapse = ", ")
  )
}

if (anyDuplicated(final_footprints$node_id) > 0) {
  stop("Final footprint dataset contains duplicate node IDs.")
}

# Create a non-spatial summary table documenting the final footprint geometry
# and the method used to obtain each polygon.
#
# The footprint_method column records how each final polygon was obtained,
# which provides a clear audit trail for the few nodes requiring
# special treatment.

final_list <- final_footprints %>%
  st_drop_geometry() %>%
  arrange(chain, node_id) %>%
  select(
    node_id,
    name,
    chain,
    area_km2,
    perimeter_km
  ) %>%
  mutate(
    footprint_method = case_when(
      node_id == "L02" ~
        "NETVOLC - 1.0 degree search window",
      node_id %in% c("T15", "T16", "T17") ~
        "Shared Britannia NETVOLC footprint + nearest-node partition",
      TRUE ~
        "NETVOLC - 0.75 degree search window"
    )
  )

# Save

st_write(
  final_footprints,
  file.path(
    "Outputs",
    "final_seamount_footprints.gpkg"
  ),
  layer = "seamount_footprints",
  delete_dsn = TRUE
)

write.csv(
  final_list,
  "Outputs/final_seamount_footprint_list.csv",
  row.names = FALSE
)

saveRDS(
  final_footprints,
  "Outputs/final_seamount_footprints.rds"
)

write.csv(
  footprint_choices,
  "Outputs/netvolc_final_footprint_choices.csv",
  row.names = FALSE
)

