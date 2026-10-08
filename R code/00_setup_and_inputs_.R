# Taupo Bank Analysis pt. 1 :)
# Data cleaning and setup
#
# Purpose:
# Load and clean the core project inputs used throughout the analysis.

library(terra)
library(sf)
library(dplyr)
library(readr)
library(ncdf4)


# File paths

gebco_file <- paste0(
  "GEBCO_23_Sep_2026_4a55cbf62670/",
  "gebco_2026_n-16.0_s-38.0_w154.0_e161.0.nc"
)

particle_file <- "tasman_seamount_2023_merged.nc"


# Read node tables

tasmantid_nodes <- read_csv(
  "tasmantid_nodes.csv",
  show_col_types = FALSE
)

lord_howe_nodes <- read_csv(
  "lord_howe_nodes.csv",
  show_col_types = FALSE
)


# Combine node tables

candidate_nodes <- bind_rows(
  tasmantid_nodes,
  lord_howe_nodes
)


# Clean node-table text fields

candidate_nodes <- candidate_nodes %>%
  mutate(
    node_id = trimws(node_id),
    name = trimws(name),
    chain = trimws(chain),
    feature = stringr::str_to_sentence(
      trimws(feature)
    ),
    status = tolower(
      trimws(status)
    ),
    notes = as.character(notes)
  )


# Standardise Middleton Reef node ID.
# Original CSV uses L9; analysis uses L09.

candidate_nodes <- candidate_nodes %>%
  mutate(
    node_id = if_else(
      node_id == "L9",
      "L09",
      node_id
    )
  )


# Standardise Recorder Guyot names.

candidate_nodes <- candidate_nodes %>%
  mutate(
    name = case_when(
      node_id == "T06" ~ "North Recorder Guyot",
      node_id == "T07" ~ "South Recorder Guyot",
      TRUE ~ name
    )
  )


# Update the Chesterfield Plateau working centre.
# The original QGIS point was visually checked against GEBCO and
# repositioned to the centre of the main plateau body.

candidate_nodes <- candidate_nodes %>%
  mutate(
    lon = if_else(
      node_id == "L02",
      158.600400,
      lon
    ),
    lat = if_else(
      node_id == "L02",
      -19.261420,
      lat
    )
  )


# Keep the original feature description
# and create a standardised working feature category.

candidate_nodes <- candidate_nodes %>%
  rename(
    feature_original = feature
  ) %>%
  mutate(
    feature_class = case_when(

      # Banks
      node_id %in% c(
        "L03", "L04", "L05", "L06",
        "T21", "T22"
      ) ~ "Bank",

      # Guyots
      node_id %in% c(
        "L07",
        "T06", "T07",
        "T12", "T13",
        "T14", "T15", "T16", "T17",
        "T20"
      ) ~ "Guyot",

      # Reef-capped features
      node_id %in% c(
        "T01", "T02", "T03",
        "L09", "L10"
      ) ~ "Reef-capped feature",

      # Emergent volcanic edifices
      node_id %in% c(
        "L11", "L12"
      ) ~ "Island-capped volcanic edifice",

      # Chesterfield Plateau
      node_id == "L02" ~
        "Volcanic edifice on plateau",

      # Remaining submarine volcanic features
      TRUE ~ "Seamount"
    )
  )


# Check node-table integrity.

if (any(is.na(candidate_nodes$node_id))) {
  stop("Missing node IDs found.")
}

if (
  any(is.na(candidate_nodes$lon)) ||
  any(is.na(candidate_nodes$lat))
) {
  stop("Missing node coordinates found.")
}

if (anyDuplicated(candidate_nodes$node_id) > 0) {
  stop("Duplicate node IDs found.")
}


# Convert cleaned node coordinates to sf points.

candidate_nodes_sf <- st_as_sf(
  candidate_nodes,
  coords = c(
    "lon",
    "lat"
  ),
  crs = 4326,
  remove = FALSE
)


# Load GEBCO bathymetry.

gebco <- rast(
  gebco_file,
  subds = "elevation"
)


