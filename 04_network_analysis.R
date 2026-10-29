# 04 Depth-resolved connectivity and network analysis
#
# Purpose:
# Implements all network analysis proposed in Section 3.2
#
# Core analysis:
# - depth-resolved directed connectivity among seamount-depth nodes;
# - external degree and strength, source-sink index, retention and isolation;
# - weighted directed betweenness and Infomap subnetworks;
# - source-normalised depth transitions and chain-to-chain exchange;
# - depth-integrated connectivity using equal weighting across represented
#   source depths;
# - comparison of depth-resolved and depth-integrated node roles;
# - removal of cross-depth edges to quantify their structural contribution; and
# - sensitivity to PLD, pre-competency, release timing and optional edge
#   thresholds.
#
# Interpretation:
# The Taupo pilot contains releases from only one source seamount, so source/sink
# and whole-network role metrics are workflow demonstrations rather than final
# ecological inferences. They become fully interpretable once releases are
# available from all study seamounts.
#

library(dplyr)
library(tibble)
library(igraph)

# Depth strata

depth_levels <- c(
  "200-700 m",
  "700-1500 m",
  "1500-2500 m",
  ">2500 m"
)

# Inputs from habitat-encounter analysis (03)

settlement_events <- readRDS(
  file.path(
    "Outputs",
    "settlement_events.rds"
  )
)

release_counts <- readRDS(
  file.path(
    "Outputs",
    "release_counts.rds"
  )
)

network_nodes <- readRDS(
  file.path(
    "Outputs",
    "network_nodes.rds"
  )
)

# Site information
#
# Reconstruct site-level metadata from the saved network-node inventory.
# network_nodes was created in 03 from the mapped benthic surfaces, with
# one node for each seamount-depth combination containing eligible habitat.
# This keeps the network analysis self-contained once the saved node inventory exists.

candidate_nodes <- network_nodes %>%
  select(
    node_id,
    name,
    chain
  ) %>%
  distinct()

# Helpers

make_multilayer_node_id <- function(
    node_id,
    depth_stratum) {

  paste(
    node_id,
    depth_stratum,
    sep = "__"
  )
}

make_period_label <- function(
    release_id) {

  paste0(
    "release_",
    release_id
  )
}

# Austral-season helper
#
# The full study releases particles monthly and then groups those releases
# into austral seasons for the release-timing sensitivity analysis.
#
# The Taupo pilot contains only June-August releases, so all pilot
# releases belong to winter. We therefore leave release_calendar as NULL for
# the pilot and retain release-to-release comparisons.
#
# When the full 12-month simulation is available, populate release_calendar
# with the release_id used by Dustin's NetCDF file and its calendar month.
# Seasonal pooled analyses will then be created automatically.
#
# Example only:
#
# release_calendar <- tibble(
#   release_id = 1:12,
#   month_number = 1:12
# )
#
# Do not assume that Dustin's release_id values equal calendar month numbers;
# use the actual release metadata when the full trajectory file is supplied.

month_to_austral_season <- function(
    month_number) {

  case_when(
    month_number %in%
      c(
        12,
        1,
        2
      ) ~
      "summer",
    month_number %in%
      c(
        3,
        4,
        5
      ) ~
      "autumn",
    month_number %in%
      c(
        6,
        7,
        8
      ) ~
      "winter",
    month_number %in%
      c(
        9,
        10,
        11
      ) ~
      "spring",
    TRUE ~
      NA_character_
  )
}

# Pilot default.
#
# Leave as NULL until a verified release_id-to-month lookup is available.

release_calendar <- NULL

# Analysis settings
#
# Infomap uses repeated optimisation trials. A fixed seed makes the workflow
# reproducible; both settings are saved with the analysis outputs.

infomap_trials <- 100L
infomap_seed <- 20261029L

# Minimum dispersal-probability thresholds for threshold sensitivity.
#
# Leave empty for the Taupo pilot. Populate after the full particle ensemble is
# available and the threshold range has been selected, as described in S4.

connectivity_thresholds <- numeric(0)

make_scenario_name <- function(
    analysis_period,
    pld_days,
    precompetency_days) {

  paste0(
    analysis_period,
    "__PLD",
    pld_days,
    "__PC",
    precompetency_days,
    "d"
  )
}

# Release periods, seasonal periods and pooled data
#
# Three time aggregations can be carried through the same analysis:
#
# 1. individual releases ("release") for the Taupo pilot and diagnostic checks;
# 2. austral seasons ("season") when a verified release calendar is supplied;
# 3. all releases pooled ("pooled") for the overall connectivity estimate.
#
# Connectivity is always recalculated from particle counts after aggregation.
# Probabilities are never averaged across months or seasons.

prepare_analysis_periods <- function(
    settlement_events,
    release_counts,
    release_calendar = NULL) {

  # Individual release periods

  release_settlements <- settlement_events %>%
    mutate(
      analysis_period =
        make_period_label(
          release_id
        ),
      period_type =
        "release"
    )

  release_releases <- release_counts %>%
    mutate(
      analysis_period =
        make_period_label(
          release_id
        ),
      period_type =
        "release"
    )

  # Overall pooled period

  pooled_settlements <- settlement_events %>%
    mutate(
      analysis_period =
        "pooled",
      period_type =
        "pooled"
    )

  pooled_releases <- release_counts %>%
    group_by(
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    summarise(
      N_released =
        sum(
          N_released
        ),
      .groups =
        "drop"
    ) %>%
    mutate(
      release_id =
        NA_integer_,
      analysis_period =
        "pooled",
      period_type =
        "pooled"
    )

  # Optional seasonal pooling
  #
  # Seasonal rows are created only when a verified release calendar is
  # supplied. This prevents the pilot code from guessing which release_id
  # corresponds to which month.

  seasonal_settlements <-
    tibble()

  seasonal_releases <-
    tibble()

  if (!is.null(
    release_calendar
  )) {

    required_calendar_columns <-
      c(
        "release_id",
        "month_number"
      )

    if (!all(
      required_calendar_columns %in%
        names(
          release_calendar
        )
    )) {

      stop(
        paste(
          "release_calendar must contain:",
          paste(
            required_calendar_columns,
            collapse = ", "
          )
        )
      )
    }

    release_calendar_checked <- release_calendar %>%
      distinct(
        release_id,
        month_number
      ) %>%
      mutate(
        austral_season =
          month_to_austral_season(
            month_number
          )
      )

    if (any(
      is.na(
        release_calendar_checked$austral_season
      )
    )) {

      stop(
        "release_calendar contains a month_number outside 1-12."
      )
    }

    # Every release represented in the analysis must have a calendar entry.

    missing_release_ids <-
      setdiff(
        unique(
          release_counts$release_id
        ),
        release_calendar_checked$release_id
      )

    if (length(
      missing_release_ids
    ) > 0) {

      stop(
        paste(
          "No calendar month supplied for release_id:",
          paste(
            missing_release_ids,
            collapse = ", "
          )
        )
      )
    }

    seasonal_settlements <- settlement_events %>%
      left_join(
        release_calendar_checked,
        by =
          "release_id"
      ) %>%
      mutate(
        analysis_period =
          paste0(
            "season_",
            austral_season
          ),
        period_type =
          "season"
      ) %>%
      select(
        -month_number,
        -austral_season
      )

    # Sum release counts before calculating seasonal probabilities.
    #
    # This is important: a seasonal connectivity probability is based on the
    # total number of successful events divided by the total number released
    # in that season, not the mean of monthly connectivity probabilities.

    seasonal_releases <- release_counts %>%
      left_join(
        release_calendar_checked,
        by =
          "release_id"
      ) %>%
      group_by(
        austral_season,
        pld_days,
        precompetency_days,
        source_node,
        source_depth_stratum
      ) %>%
      summarise(
        N_released =
          sum(
            N_released
          ),
        .groups =
          "drop"
      ) %>%
      mutate(
        release_id =
          NA_integer_,
        analysis_period =
          paste0(
            "season_",
            austral_season
          ),
        period_type =
          "season"
      ) %>%
      select(
        -austral_season
      )
  }

  list(
    settlements =
      bind_rows(
        release_settlements,
        seasonal_settlements,
        pooled_settlements
      ),
    releases =
      bind_rows(
        release_releases,
        seasonal_releases,
        pooled_releases
      )
  )
}

# Depth-resolved connectivity
#
# P_(i,alpha -> j,beta)
# =
# N_(i,alpha -> j,beta) / N_(i,alpha)
#
# Connectivity is the proportion of particles released from source
# seamount-depth node (i,alpha) whose first qualifying potential habitat-encounter
# event occurs at destination seamount-depth node (j,beta).
#
# Keep zero-connectivity destinations so a zero means the source was
# simulated but no particles produced a potential habitat-encounter event there.

calculate_connectivity <- function(
    settlements_period,
    releases_period,
    network_nodes) {

  settlement_counts <- settlements_period %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum,
      destination_node,
      destination_depth_stratum,
      name = "N_settled"
    )

  destination_nodes <- network_nodes %>%
    transmute(
      destination_node =
        node_id,
      destination_depth_stratum =
        depth_stratum
    ) %>%
    distinct()

  dplyr::cross_join(
    releases_period,
    destination_nodes
  ) %>%
    left_join(
      settlement_counts,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node",
        "source_depth_stratum",
        "destination_node",
        "destination_depth_stratum"
      )
    ) %>%
    mutate(
      N_settled =
        coalesce(
          N_settled,
          0L
        ),
      connectivity =
        N_settled /
        N_released
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days,
      source_node,
      source_depth_stratum,
      destination_node,
      destination_depth_stratum
    )
}

# Total potential habitat-encounter probability
#
# Because 03 keeps only the first qualifying potential habitat-encounter event per
# particle and biological scenario, summing outgoing connectivity gives the
# proportion of released particles with a qualifying event somewhere in the network.

calculate_settlement_summary <- function(
    connectivity) {

  connectivity %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    summarise(
      total_settlement_probability =
        sum(
          connectivity
        ),
      unsettled_probability =
        1 -
        total_settlement_probability,
      .groups =
        "drop"
    )
}

# Directionality
#
# Connectivity is directional, so P_(i,alpha -> j,beta) does not have
# to equal P_(j,beta -> i,alpha).
#
# Reverse connectivity remains NA when the reverse source was not
# represented by releases.

calculate_directionality <- function(
    connectivity) {

  reverse_table <- connectivity %>%
    transmute(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node =
        destination_node,
      source_depth_stratum =
        destination_depth_stratum,
      destination_node =
        source_node,
      destination_depth_stratum =
        source_depth_stratum,
      reverse_connectivity =
        connectivity
    )

  connectivity %>%
    left_join(
      reverse_table,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node",
        "source_depth_stratum",
        "destination_node",
        "destination_depth_stratum"
      )
    ) %>%
    mutate(
      reverse_available =
        !is.na(
          reverse_connectivity
        ),
      directional_difference =
        connectivity -
        reverse_connectivity,
      absolute_directional_difference =
        abs(
          directional_difference
        )
    )
}

# Connectivity roles
#
# External out-strength:
# s_OUT(i,alpha) = sum P_(i,alpha -> j,beta), j != i
#
# External in-strength:
# s_IN(i,alpha) = sum P_(j,beta -> i,alpha), j != i
#
# Same-seamount potential habitat encounter is kept separate as retention.

calculate_strength_retention_isolation <- function(
    connectivity,
    releases_period,
    network_nodes) {

  # External connectivity excludes all same-seamount connections.
  #
  # This keeps export/import to the wider seamount network separate from
  # same-site retention, including cross-depth return on the source seamount.

  external_edges <- connectivity %>%
    filter(
      source_node !=
        destination_node
    )

  # External out-degree
  #
  # Number of non-zero outgoing connections from a focal site-depth node
  # to site-depth habitats on other seamounts.

  out_degree <- external_edges %>%
    filter(
      connectivity > 0
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    summarise(
      k_out =
        n(),
      .groups =
        "drop"
    ) %>%
    rename(
      node_id =
        source_node,
      depth_stratum =
        source_depth_stratum
    )

  # External in-degree
  #
  # Number of non-zero incoming connections from site-depth habitats
  # on other seamounts.

  in_degree <- external_edges %>%
    filter(
      connectivity > 0
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      destination_node,
      destination_depth_stratum
    ) %>%
    summarise(
      k_in =
        n(),
      .groups =
        "drop"
    ) %>%
    rename(
      node_id =
        destination_node,
      depth_stratum =
        destination_depth_stratum
    )

  # External out-strength
  #
  # s_OUT(i,alpha) = sum P_(i,alpha -> j,beta), j != i

  out_strength <- external_edges %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    summarise(
      s_out =
        sum(
          connectivity
        ),
      .groups =
        "drop"
    ) %>%
    rename(
      node_id =
        source_node,
      depth_stratum =
        source_depth_stratum
    )

  # External in-strength
  #
  # s_IN(i,alpha) = sum P_(j,beta -> i,alpha), j != i

  in_strength <- external_edges %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      destination_node,
      destination_depth_stratum
    ) %>%
    summarise(
      s_in =
        sum(
          connectivity
        ),
      .groups =
        "drop"
    ) %>%
    rename(
      node_id =
        destination_node,
      depth_stratum =
        destination_depth_stratum
    )

  # Same-site retention
  #
  # R_(i,alpha) = sum_beta P_(i,alpha -> i,beta)
  #
  # This includes return to the same depth stratum and cross-depth return
  # to another depth stratum on the same seamount.

  retention <- connectivity %>%
    filter(
      source_node ==
        destination_node
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    summarise(
      retention =
        sum(
          connectivity
        ),
      .groups =
        "drop"
    ) %>%
    rename(
      node_id =
        source_node,
      depth_stratum =
        source_depth_stratum
    )

  # Nodes represented by releases
  #
  # The full study releases particles from all sampled source habitats.
  # In the Taupo pilot only Taupo source-depth nodes are represented by releases,
  # so node-role comparisons among seamounts are not interpreted from the pilot.

  role_nodes <- releases_period %>%
    transmute(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      node_id =
        source_node,
      depth_stratum =
        source_depth_stratum,
      N_released
    ) %>%
    distinct() %>%
    left_join(
      network_nodes %>%
        select(
          node_id,
          name,
          chain,
          depth_stratum
        ) %>%
        distinct(),
      by = c(
        "node_id",
        "depth_stratum"
      )
    )

  role_nodes %>%
    left_join(
      out_degree,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    left_join(
      in_degree,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    left_join(
      out_strength,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    left_join(
      in_strength,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    left_join(
      retention,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    mutate(
      k_out =
        coalesce(
          k_out,
          0L
        ),
      k_in =
        coalesce(
          k_in,
          0L
        ),
      s_out =
        coalesce(
          s_out,
          0
        ),
      s_in =
        coalesce(
          s_in,
          0
        ),
      retention =
        coalesce(
          retention,
          0
        ),
      source_sink_index =
        if_else(
          s_out + s_in > 0,
          (s_out - s_in) /
            (s_out + s_in),
          NA_real_
        ),
      # Total external strength is used as the inverse indicator of
      # relative isolation: lower values indicate greater isolation.
      s_external =
        s_out +
        s_in
    )
}

# Weighted shortest-path cost
#
# Larger connectivity probabilities represent stronger connections.
# Convert probabilities into additive costs for shortest-path analysis.

connectivity_to_cost <- function(
    connectivity) {

  -log(
    connectivity
  )
}

# Betweenness centrality

calculate_graph_betweenness <- function(
    graph) {

  if (ecount(graph) == 0) {

    return(
      tibble(
        graph_node =
          V(graph)$name,
        betweenness =
          0,
        betweenness_normalised =
          0
      )
    )
  }

  tibble(
    graph_node =
      V(graph)$name,
    betweenness =
      igraph::betweenness(
        graph,
        directed = TRUE,
        weights =
          E(graph)$path_cost,
        normalized = FALSE
      ),
    betweenness_normalised =
      igraph::betweenness(
        graph,
        directed = TRUE,
        weights =
          E(graph)$path_cost,
        normalized = TRUE
      )
  )
}

# Multilayer graph
#
# Each vertex is one seamount-depth node.
# Within-depth and cross-depth edges are retained.
#
# Exact self-loops are removed for path-based metrics, but
# same-seamount cross-depth connections remain.

build_multilayer_graph <- function(
    connectivity,
    network_nodes,
    analysis_period_value,
    pld_value,
    precompetency_value) {

  vertices <- network_nodes %>%
    distinct(
      node_id,
      name,
      chain,
      depth_stratum
    ) %>%
    mutate(
      graph_node =
        make_multilayer_node_id(
          node_id,
          depth_stratum
        )
    ) %>%
    rename(
      seamount_name =
        name
    ) %>%
    transmute(
      name =
        graph_node,
      node_id,
      seamount_name,
      chain,
      depth_stratum
    )

  edges <- connectivity %>%
    filter(
      analysis_period ==
        analysis_period_value,
      pld_days ==
        pld_value,
      precompetency_days ==
        precompetency_value,
      connectivity >
        0
    ) %>%
    mutate(
      from =
        make_multilayer_node_id(
          source_node,
          source_depth_stratum
        ),
      to =
        make_multilayer_node_id(
          destination_node,
          destination_depth_stratum
        ),
      path_cost =
        pmax(
          connectivity_to_cost(
            connectivity
          ),
          .Machine$double.eps
        )
    ) %>%
    filter(
      from !=
        to
    ) %>%
    select(
      from,
      to,
      connectivity,
      path_cost
    )

  graph_from_data_frame(
    d =
      edges,
    directed =
      TRUE,
    vertices =
      vertices
  )
}

# Betweenness for every release-specific and pooled scenario

calculate_all_betweenness <- function(
    connectivity,
    network_nodes) {

  scenarios <- connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days
    )

  output <- vector(
    "list",
    nrow(
      scenarios
    )
  )

  for (i in seq_len(
    nrow(
      scenarios
    )
  )) {

    scenario <-
      scenarios[i, ]

    graph <- build_multilayer_graph(
      connectivity =
        connectivity,
      network_nodes =
        network_nodes,
      analysis_period_value =
        scenario$analysis_period,
      pld_value =
        scenario$pld_days,
      precompetency_value =
        scenario$precompetency_days
    )

    vertices <- igraph::as_data_frame(
      graph,
      what = "vertices"
    )

    betweenness_table <-
      calculate_graph_betweenness(
        graph
      )

    output[[i]] <- vertices %>%
      left_join(
        betweenness_table,
        by = c(
          "name" =
            "graph_node"
        )
      ) %>%
      transmute(
        analysis_period =
          scenario$analysis_period,
        period_type =
          scenario$period_type,
        pld_days =
          scenario$pld_days,
        precompetency_days =
          scenario$precompetency_days,
        node_id,
        name =
          seamount_name,
        chain,
        depth_stratum,
        betweenness,
        betweenness_normalised
      )
  }

  bind_rows(
    output
  )
}

# Combine strengths, retention, SSI, isolation and betweenness

calculate_network_roles <- function(
    connectivity,
    releases_period,
    network_nodes) {

  roles <-
    calculate_strength_retention_isolation(
      connectivity =
        connectivity,
      releases_period =
        releases_period,
      network_nodes =
        network_nodes
    )

  betweenness <-
    calculate_all_betweenness(
      connectivity =
        connectivity,
      network_nodes =
        network_nodes
    )

  roles %>%
    left_join(
      betweenness %>%
        select(
          analysis_period,
          period_type,
          pld_days,
          precompetency_days,
          node_id,
          depth_stratum,
          betweenness,
          betweenness_normalised
        ),
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    mutate(
      betweenness =
        coalesce(
          betweenness,
          0
        ),
      betweenness_normalised =
        coalesce(
          betweenness_normalised,
          0
        )
    )
}

# Multilayer connectivity matrix
#
# Rows = source seamount-depth nodes
# Columns = destination seamount-depth nodes
#
# 0 = source simulated but no qualifying habitat encounter
# NA = source not released in that scenario

make_multilayer_matrix <- function(
    connectivity,
    network_nodes,
    analysis_period_value,
    pld_value,
    precompetency_value) {

  node_ids <- network_nodes %>%
    distinct(
      node_id,
      depth_stratum
    ) %>%
    mutate(
      depth_order =
        match(
          depth_stratum,
          depth_levels
        )
    ) %>%
    arrange(
      node_id,
      depth_order
    ) %>%
    transmute(
      graph_node =
        make_multilayer_node_id(
          node_id,
          depth_stratum
        )
    ) %>%
    pull(
      graph_node
    )

  matrix_out <- matrix(
    NA_real_,
    nrow =
      length(
        node_ids
      ),
    ncol =
      length(
        node_ids
      ),
    dimnames = list(
      source =
        node_ids,
      destination =
        node_ids
    )
  )

  dat <- connectivity %>%
    filter(
      analysis_period ==
        analysis_period_value,
      pld_days ==
        pld_value,
      precompetency_days ==
        precompetency_value
    ) %>%
    mutate(
      source_graph_node =
        make_multilayer_node_id(
          source_node,
          source_depth_stratum
        ),
      destination_graph_node =
        make_multilayer_node_id(
          destination_node,
          destination_depth_stratum
        )
    )

  released_sources <-
    unique(
      dat$source_graph_node
    )

  matrix_out[
    released_sources,
    node_ids
  ] <- 0

  for (i in seq_len(
    nrow(
      dat
    )
  )) {

    matrix_out[
      dat$source_graph_node[i],
      dat$destination_graph_node[i]
    ] <-
      dat$connectivity[i]
  }

  matrix_out
}

# Depth-block matrices
#
# One seamount-to-seamount matrix for each source-depth to
# destination-depth combination.
#
# alpha = beta is within-depth.
# alpha != beta is cross-depth.

make_depth_blocks <- function(
    connectivity,
    network_nodes,
    candidate_nodes,
    analysis_period_value,
    pld_value,
    precompetency_value) {

  seamount_ids <-
    candidate_nodes$node_id

  blocks <-
    list()

  scenario_data <- connectivity %>%
    filter(
      analysis_period ==
        analysis_period_value,
      pld_days ==
        pld_value,
      precompetency_days ==
        precompetency_value
    )

  for (alpha in depth_levels) {

    for (beta in depth_levels) {

      block <- matrix(
        NA_real_,
        nrow =
          length(
            seamount_ids
          ),
        ncol =
          length(
            seamount_ids
          ),
        dimnames = list(
          source =
            seamount_ids,
          destination =
            seamount_ids
        )
      )

      valid_destinations <- network_nodes %>%
        filter(
          depth_stratum ==
            beta
        ) %>%
        pull(
          node_id
        )

      dat <- scenario_data %>%
        filter(
          source_depth_stratum ==
            alpha,
          destination_depth_stratum ==
            beta
        )

      released_sources <-
        unique(
          dat$source_node
        )

      if (
        length(
          released_sources
        ) > 0 &&
        length(
          valid_destinations
        ) > 0
      ) {

        block[
          released_sources,
          valid_destinations
        ] <- 0
      }

      if (nrow(dat) > 0) {

        for (i in seq_len(
          nrow(
            dat
          )
        )) {

          block[
            dat$source_node[i],
            dat$destination_node[i]
          ] <-
            dat$connectivity[i]
        }
      }

      blocks[[
        paste(
          alpha,
          "to",
          beta,
          sep = "__"
        )
      ]] <-
        block
    }
  }

  blocks
}

# Build multilayer matrices and all 16 depth blocks
# for every analysis scenario.

build_all_matrices <- function(
    connectivity,
    network_nodes,
    candidate_nodes) {

  scenarios <- connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days
    )

  multilayer <-
    list()

  depth_blocks <-
    list()

  for (i in seq_len(
    nrow(
      scenarios
    )
  )) {

    scenario <-
      scenarios[i, ]

    scenario_name <- make_scenario_name(
      scenario$analysis_period,
      scenario$pld_days,
      scenario$precompetency_days
    )

    multilayer[[scenario_name]] <-
      make_multilayer_matrix(
        connectivity =
          connectivity,
        network_nodes =
          network_nodes,
        analysis_period_value =
          scenario$analysis_period,
        pld_value =
          scenario$pld_days,
        precompetency_value =
          scenario$precompetency_days
      )

    depth_blocks[[scenario_name]] <-
      make_depth_blocks(
        connectivity =
          connectivity,
        network_nodes =
          network_nodes,
        candidate_nodes =
          candidate_nodes,
        analysis_period_value =
          scenario$analysis_period,
        pld_value =
          scenario$pld_days,
        precompetency_value =
          scenario$precompetency_days
      )
  }

  list(
    scenario_index =
      scenarios %>%
      mutate(
        scenario_name =
          make_scenario_name(
            analysis_period,
            pld_days,
            precompetency_days
          )
      ),
    multilayer =
      multilayer,
    depth_blocks =
      depth_blocks
  )
}

# Connectivity among depth strata
#
# T_(alpha -> beta) is calculated from source-normalised dispersal
# probabilities, following Supplementary Methods S5.
#
# For each source seamount represented in source stratum alpha, connectivity to
# all destination seamounts in destination stratum beta is summed. These
# source-seamount values are then averaged so each represented source
# site-depth habitat contributes equal weight irrespective of particle number.

calculate_depth_connectivity <- function(
    connectivity) {

  source_nodes <- connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    )

  source_to_depth <- connectivity %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum,
      destination_depth_stratum
    ) %>%
    summarise(
      source_depth_connectivity =
        sum(
          connectivity
        ),
      .groups =
        "drop"
    )

  dplyr::cross_join(
    source_nodes,
    tibble(
      destination_depth_stratum =
        depth_levels
    )
  ) %>%
    left_join(
      source_to_depth,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node",
        "source_depth_stratum",
        "destination_depth_stratum"
      )
    ) %>%
    mutate(
      source_depth_connectivity =
        coalesce(
          source_depth_connectivity,
          0
        )
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_depth_stratum,
      destination_depth_stratum
    ) %>%
    summarise(
      depth_connectivity =
        mean(
          source_depth_connectivity
        ),
      n_source_seamounts =
        n_distinct(
          source_node
        ),
      .groups =
        "drop"
    ) %>%
    mutate(
      depth_transition =
        if_else(
          source_depth_stratum == destination_depth_stratum,
          "within-depth",
          "between-depth"
        ),
      source_depth_stratum =
        factor(
          source_depth_stratum,
          levels = depth_levels
        ),
      destination_depth_stratum =
        factor(
          destination_depth_stratum,
          levels = depth_levels
        )
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days,
      precompetency_days,
      source_depth_stratum,
      destination_depth_stratum
    )
}

# Overall within-depth versus between-depth potential habitat encounter
#
# Both are calculated relative to all particles released
# in each scenario.

calculate_within_between_depth <- function(
    settlements_period,
    releases_period) {

  total_releases <- releases_period %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    summarise(
      N_released =
        sum(
          N_released
        ),
      .groups =
        "drop"
    )

  settlement_counts <- settlements_period %>%
    mutate(
      depth_transition =
        if_else(
          source_depth_stratum ==
            destination_depth_stratum,
          "within-depth",
          "between-depth"
        )
    ) %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      depth_transition,
      name = "N_settled"
    )

  dplyr::cross_join(
    total_releases,
    tibble(
      depth_transition =
        c(
          "within-depth",
          "between-depth"
        )
    )
  ) %>%
    left_join(
      settlement_counts,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "depth_transition"
      )
    ) %>%
    mutate(
      N_settled =
        coalesce(
          N_settled,
          0L
        ),
      settlement_fraction =
        N_settled /
        N_released
    )
}

# Chain-level connectivity
#
# Chain aggregation follows the two-stage source-normalised procedure in S4.
# First, connectivity from each source seamount to each receiving chain is
# averaged across the source depth strata represented on that seamount. Second,
# these seamount-level values are averaged across source seamounts in the source
# chain. Same-seamount retention is excluded throughout.

calculate_chain_connectivity <- function(
    connectivity,
    candidate_nodes) {

  chain_lookup <- candidate_nodes %>%
    select(
      node_id,
      chain
    ) %>%
    distinct()

  destination_chains <- tibble(
    destination_chain =
      sort(
        unique(
          chain_lookup$chain
        )
      )
  )

  represented_source_depths <- connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    ) %>%
    left_join(
      chain_lookup %>%
        rename(
          source_chain = chain
        ),
      by = c(
        "source_node" = "node_id"
      )
    )

  source_depth_to_chain <- connectivity %>%
    filter(
      source_node != destination_node
    ) %>%
    left_join(
      chain_lookup %>%
        rename(
          destination_chain = chain
        ),
      by = c(
        "destination_node" = "node_id"
      )
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum,
      destination_chain
    ) %>%
    summarise(
      source_depth_to_chain =
        sum(
          connectivity
        ),
      .groups = "drop"
    )

  source_depth_complete <- dplyr::cross_join(
    represented_source_depths,
    destination_chains
  ) %>%
    left_join(
      source_depth_to_chain,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node",
        "source_depth_stratum",
        "destination_chain"
      )
    ) %>%
    mutate(
      source_depth_to_chain =
        coalesce(
          source_depth_to_chain,
          0
        )
    )

  source_seamount_to_chain <- source_depth_complete %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_chain,
      destination_chain
    ) %>%
    summarise(
      source_seamount_connectivity =
        mean(
          source_depth_to_chain
        ),
      n_source_depth_strata =
        n_distinct(
          source_depth_stratum
        ),
      .groups = "drop"
    )

  source_seamount_to_chain %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_chain,
      destination_chain
    ) %>%
    summarise(
      chain_connectivity =
        mean(
          source_seamount_connectivity
        ),
      n_source_seamounts =
        n_distinct(
          source_node
        ),
      .groups = "drop"
    ) %>%
    mutate(
      chain_transition =
        if_else(
          source_chain == destination_chain,
          "within-chain",
          "between-chain"
        )
    )
}

# Chain directionality
#
# Compares P_(A -> B) with P_(B -> A).
#
# The reverse value remains NA when the reverse source chain was not
# represented by releases.

calculate_chain_directionality <- function(
    chain_connectivity) {

  reverse_table <- chain_connectivity %>%
    transmute(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_chain =
        destination_chain,
      destination_chain =
        source_chain,
      reverse_chain_connectivity =
        chain_connectivity
    )

  chain_connectivity %>%
    left_join(
      reverse_table,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_chain",
        "destination_chain"
      )
    ) %>%
    mutate(
      reverse_available =
        !is.na(
          reverse_chain_connectivity
        ),
      directional_difference =
        chain_connectivity -
        reverse_chain_connectivity,
      absolute_directional_difference =
        abs(
          directional_difference
        )
    )
}

# Depth-integrated connectivity
#
# C^I_(i -> j) is calculated by first summing across destination depths for each
# represented source depth, then averaging those source-depth probabilities.
# This gives every represented source site-depth habitat equal weight and avoids
# giving greater influence to depths with more simulated particles.
#
# N_released is retained below as metadata describing release effort; it is not used as
# the denominator of the depth-integrated connectivity value.

calculate_collapsed_connectivity <- function(
    connectivity,
    releases_period,
    candidate_nodes) {

  represented_source_depths <- releases_period %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum
    )

  destination_nodes <- candidate_nodes %>%
    distinct(
      destination_node = node_id
    )

  source_depth_to_destination <- connectivity %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      source_depth_stratum,
      destination_node
    ) %>%
    summarise(
      source_depth_connectivity =
        sum(
          connectivity
        ),
      .groups = "drop"
    )

  integrated <- dplyr::cross_join(
    represented_source_depths,
    destination_nodes
  ) %>%
    left_join(
      source_depth_to_destination,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node",
        "source_depth_stratum",
        "destination_node"
      )
    ) %>%
    mutate(
      source_depth_connectivity =
        coalesce(
          source_depth_connectivity,
          0
        )
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      destination_node
    ) %>%
    summarise(
      connectivity =
        mean(
          source_depth_connectivity
        ),
      n_source_depth_strata =
        n_distinct(
          source_depth_stratum
        ),
      .groups = "drop"
    )

  release_metadata <- releases_period %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node
    ) %>%
    summarise(
      N_released =
        sum(
          N_released
        ),
      .groups = "drop"
    )

  integrated %>%
    left_join(
      release_metadata,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node"
      )
    )
}

# Directionality in the depth-collapsed network

calculate_collapsed_directionality <- function(
    collapsed_connectivity) {

  reverse_table <- collapsed_connectivity %>%
    transmute(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node =
        destination_node,
      destination_node =
        source_node,
      reverse_connectivity =
        connectivity
    )

  collapsed_connectivity %>%
    left_join(
      reverse_table,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_node",
        "destination_node"
      )
    ) %>%
    mutate(
      reverse_available =
        !is.na(
          reverse_connectivity
        ),
      directional_difference =
        connectivity -
        reverse_connectivity,
      absolute_directional_difference =
        abs(
          directional_difference
        )
    )
}

# Network roles in the depth-collapsed comparison
#
# The same definitions are used as in the multilayer network:
#
# - external strengths exclude i = j
# - retention contains i = j
# - betweenness uses -log(P) edge costs

calculate_collapsed_roles <- function(
    collapsed_connectivity,
    candidate_nodes) {

  scenarios <- collapsed_connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days
    )

  output <- vector(
    "list",
    nrow(
      scenarios
    )
  )

  for (k in seq_len(
    nrow(
      scenarios
    )
  )) {

    scenario <-
      scenarios[k, ]

    dat <- collapsed_connectivity %>%
      filter(
        analysis_period ==
          scenario$analysis_period,
        pld_days ==
          scenario$pld_days,
        precompetency_days ==
          scenario$precompetency_days
      )

    external <- dat %>%
      filter(
        source_node !=
          destination_node
      )

    out_strength <- external %>%
      group_by(
        node_id =
          source_node
      ) %>%
      summarise(
        s_out =
          sum(
            connectivity
          ),
        .groups =
          "drop"
      )

    in_strength <- external %>%
      group_by(
        node_id =
          destination_node
      ) %>%
      summarise(
        s_in =
          sum(
            connectivity
          ),
        .groups =
          "drop"
      )

    retention <- dat %>%
      filter(
        source_node ==
          destination_node
      ) %>%
      transmute(
        node_id =
          source_node,
        retention =
          connectivity
      )

    vertices <- candidate_nodes %>%
      distinct(
        node_id,
        name,
        chain
      ) %>%
      rename(
        seamount_name =
          name
      ) %>%
      transmute(
        name =
          node_id,
        seamount_name,
        chain
      )

    edges <- dat %>%
      filter(
        connectivity >
          0,
        source_node !=
          destination_node
      ) %>%
      transmute(
        from =
          source_node,
        to =
          destination_node,
        connectivity,
        path_cost =
          pmax(
            connectivity_to_cost(
              connectivity
            ),
            .Machine$double.eps
          )
      )

    graph <- graph_from_data_frame(
      d =
        edges,
      directed =
        TRUE,
      vertices =
        vertices
    )

    betweenness_table <-
      calculate_graph_betweenness(
        graph
      ) %>%
      rename(
        node_id =
          graph_node
      )

    output[[k]] <- dat %>%
      distinct(
        node_id =
          source_node,
        N_released
      ) %>%
      left_join(
        candidate_nodes %>%
          distinct(
            node_id,
            name,
            chain
          ),
        by =
          "node_id"
      ) %>%
      left_join(
        out_strength,
        by =
          "node_id"
      ) %>%
      left_join(
        in_strength,
        by =
          "node_id"
      ) %>%
      left_join(
        retention,
        by =
          "node_id"
      ) %>%
      left_join(
        betweenness_table,
        by =
          "node_id"
      ) %>%
      mutate(
        analysis_period =
          scenario$analysis_period,
        period_type =
          scenario$period_type,
        pld_days =
          scenario$pld_days,
        precompetency_days =
          scenario$precompetency_days,
        s_out =
          coalesce(
            s_out,
            0
          ),
        s_in =
          coalesce(
            s_in,
            0
          ),
        retention =
          coalesce(
            retention,
            0
          ),
        s_external =
          s_out +
          s_in,
        source_sink_index =
          if_else(
            s_external > 0,
            (s_out - s_in) /
              s_external,
            NA_real_
          ),
        betweenness =
          coalesce(
            betweenness,
            0
          ),
        betweenness_normalised =
          coalesce(
            betweenness_normalised,
            0
          )
      )
  }

  bind_rows(
    output
  )
}

# Depth-collapsed matrices
#
# Rows and columns represent seamounts.
#
# Zero indicates a released source with no qualifying habitat encounter at that
# destination.
#
# Unreleased source rows remain NA.

make_collapsed_matrices <- function(
    collapsed_connectivity,
    candidate_nodes) {

  scenarios <- collapsed_connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days
    )

  node_ids <-
    candidate_nodes$node_id

  matrices <-
    list()

  for (i in seq_len(
    nrow(
      scenarios
    )
  )) {

    scenario <-
      scenarios[i, ]

    dat <- collapsed_connectivity %>%
      filter(
        analysis_period ==
          scenario$analysis_period,
        pld_days ==
          scenario$pld_days,
        precompetency_days ==
          scenario$precompetency_days
      )

    matrix_out <- matrix(
      NA_real_,
      nrow =
        length(
          node_ids
        ),
      ncol =
        length(
          node_ids
        ),
      dimnames = list(
        source =
          node_ids,
        destination =
          node_ids
      )
    )

    released_sources <-
      unique(
        dat$source_node
      )

    matrix_out[
      released_sources,
      node_ids
    ] <- 0

    for (j in seq_len(
      nrow(
        dat
      )
    )) {

      matrix_out[
        dat$source_node[j],
        dat$destination_node[j]
      ] <-
        dat$connectivity[j]
    }

    matrices[[
      make_scenario_name(
        scenario$analysis_period,
        scenario$pld_days,
        scenario$precompetency_days
      )
    ]] <-
      matrix_out
  }

  matrices
}

# Depth-resolved versus depth-integrated role comparison
#
# The proposal compares ecological roles between:
#
# 1. site-depth nodes in the multilayer network; and
# 2. the corresponding seamount in the depth-integrated network.
#
# Raw metric magnitudes are not compared directly because the two networks
# differ in node definition and size. Instead, comparisons use within-network
# ranks for metrics with a comparable ecological interpretation:
#
# - external out-strength
# - external in-strength
# - source-sink index
# - same-site retention
# - total external strength / relative isolation
# - betweenness centrality
#
# External degree is intentionally excluded because a connection is to a
# site-depth habitat in the multilayer network but to a whole seamount in the
# depth-integrated network.

safe_spearman <- function(
    x,
    y,
    minimum_n = 3) {

  keep <-
    is.finite(x) &
    is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < minimum_n ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {

    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x,
      y,
      method =
        "spearman"
    )
  )
}

make_depth_role_metric_table <- function(
    network_roles) {

  bind_rows(
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        name,
        chain,
        depth_stratum,
        metric =
          "external_out_strength",
        metric_value =
          s_out
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        name,
        chain,
        depth_stratum,
        metric =
          "external_in_strength",
        metric_value =
          s_in
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        name,
        chain,
        depth_stratum,
        metric =
          "source_sink_index",
        metric_value =
          source_sink_index
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        name,
        chain,
        depth_stratum,
        metric =
          "same_site_retention",
        metric_value =
          retention
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        name,
        chain,
        depth_stratum,
        metric =
          "total_external_strength",
        metric_value =
          s_external
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        name,
        chain,
        depth_stratum,
        metric =
          "betweenness",
        metric_value =
          betweenness_normalised
      )
  )
}

make_collapsed_role_metric_table <- function(
    collapsed_roles) {

  bind_rows(
    collapsed_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        metric =
          "external_out_strength",
        collapsed_value =
          s_out
      ),
    collapsed_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        metric =
          "external_in_strength",
        collapsed_value =
          s_in
      ),
    collapsed_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        metric =
          "source_sink_index",
        collapsed_value =
          source_sink_index
      ),
    collapsed_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        metric =
          "same_site_retention",
        collapsed_value =
          retention
      ),
    collapsed_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        metric =
          "total_external_strength",
        collapsed_value =
          s_external
      ),
    collapsed_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        metric =
          "betweenness",
        collapsed_value =
          betweenness_normalised
      )
  )
}

compare_depth_resolved_collapsed_roles <- function(
    network_roles,
    collapsed_roles) {

  depth_metrics <-
    make_depth_role_metric_table(
      network_roles
    )

  collapsed_metrics <-
    make_collapsed_role_metric_table(
      collapsed_roles
    )

  # Join the corresponding seamount between the two representations.
  #
  # Comparisons are then made separately within each depth stratum, so the
  # collapsed-network ranks are recalculated using exactly the seamounts
  # represented in that depth stratum.

  comparison_detail <- depth_metrics %>%
    inner_join(
      collapsed_metrics,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "metric"
      )
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      depth_stratum,
      metric
    ) %>%
    mutate(
      depth_resolved_rank =
        rank(
          -metric_value,
          ties.method =
            "average",
          na.last =
            "keep"
        ),
      depth_integrated_rank =
        rank(
          -collapsed_value,
          ties.method =
            "average",
          na.last =
            "keep"
        ),
      # Positive values mean the seamount has a numerically larger rank
      # in the depth-resolved network than in the integrated comparison.
      rank_difference =
        depth_resolved_rank -
        depth_integrated_rank,
      ssi_role_depth =
        case_when(
          metric != "source_sink_index" ~
            NA_character_,
          metric_value > 0 ~
            "source-like",
          metric_value < 0 ~
            "sink-like",
          TRUE ~
            "balanced"
        ),
      ssi_role_integrated =
        case_when(
          metric != "source_sink_index" ~
            NA_character_,
          collapsed_value > 0 ~
            "source-like",
          collapsed_value < 0 ~
            "sink-like",
          TRUE ~
            "balanced"
        ),
      ssi_role_change =
        if_else(
          metric ==
            "source_sink_index",
          ssi_role_depth !=
            ssi_role_integrated,
          NA
        )
    ) %>%
    ungroup()

  comparison_summary <- comparison_detail %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      depth_stratum,
      metric
    ) %>%
    summarise(
      n_seamounts =
        sum(
          is.finite(metric_value) &
            is.finite(collapsed_value)
        ),
      spearman_rho =
        safe_spearman(
          metric_value,
          collapsed_value
        ),
      n_ssi_role_changes =
        if (
          first(metric) ==
            "source_sink_index"
        ) {
          sum(
            ssi_role_change,
            na.rm =
              TRUE
          )
        } else {
          NA_integer_
        },
      .groups =
        "drop"
    )

  list(
    detail =
      comparison_detail,
    summary =
      comparison_summary
  )
}

# Pairwise sensitivity analysis of node-level network roles
#
# PLD, pre-competency and release timing are varied one at a time while the
# remaining assumptions are held constant. Spearman rank correlation is used
# to quantify whether relative site-depth-node roles are preserved.
#
# For the source-sink index, changes between source-like, sink-like and
# balanced classifications are counted explicitly.

make_sensitivity_role_metric_table <- function(
    network_roles) {

  bind_rows(
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "external_out_degree",
        metric_value =
          as.numeric(
            k_out
          )
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "external_in_degree",
        metric_value =
          as.numeric(
            k_in
          )
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "external_out_strength",
        metric_value =
          s_out
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "external_in_strength",
        metric_value =
          s_in
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "source_sink_index",
        metric_value =
          source_sink_index
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "same_site_retention",
        metric_value =
          retention
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "total_external_strength",
        metric_value =
          s_external
      ),
    network_roles %>%
      transmute(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days,
        node_id,
        depth_stratum,
        metric =
          "betweenness",
        metric_value =
          betweenness_normalised
      )
  )
}

compare_two_role_scenarios <- function(
    scenario_a,
    scenario_b) {

  joined <- scenario_a %>%
    select(
      node_id,
      depth_stratum,
      metric,
      value_a =
        metric_value
    ) %>%
    inner_join(
      scenario_b %>%
        select(
          node_id,
          depth_stratum,
          metric,
          value_b =
            metric_value
        ),
      by = c(
        "node_id",
        "depth_stratum",
        "metric"
      )
    )

  joined %>%
    group_by(
      metric
    ) %>%
    summarise(
      n_common_nodes =
        sum(
          is.finite(value_a) &
            is.finite(value_b)
        ),
      spearman_rho =
        safe_spearman(
          value_a,
          value_b
        ),
      n_ssi_role_changes =
        if (
          first(metric) ==
            "source_sink_index"
        ) {
          sum(
            case_when(
              !is.finite(value_a) |
                !is.finite(value_b) ~
                FALSE,
              value_a > 0 &
                value_b <= 0 ~
                TRUE,
              value_a < 0 &
                value_b >= 0 ~
                TRUE,
              value_a == 0 &
                value_b != 0 ~
                TRUE,
              TRUE ~
                FALSE
            ),
            na.rm =
              TRUE
          )
        } else {
          NA_integer_
        },
      .groups =
        "drop"
    )
}

calculate_role_sensitivity <- function(
    network_roles) {

  metric_table <-
    make_sensitivity_role_metric_table(
      network_roles
    )

  output <-
    list()

  output_index <-
    1L

  # PLD sensitivity
  #
  # Compare every pair of PLD values while holding analysis period and
  # pre-competency constant.

  period_values <-
    unique(
      metric_table$analysis_period
    )

  for (period_value in period_values) {

    period_data <- metric_table %>%
      filter(
        analysis_period ==
          period_value
      )

    period_type_value <-
      unique(
        period_data$period_type
      )[1]

    pc_values <-
      sort(
        unique(
          period_data$precompetency_days
        )
      )

    for (pc_value in pc_values) {

      fixed_data <- period_data %>%
        filter(
          precompetency_days ==
            pc_value
        )

      pld_values <-
        sort(
          unique(
            fixed_data$pld_days
          )
        )

      if (length(pld_values) < 2) {
        next
      }

      pld_pairs <-
        combn(
          pld_values,
          2
        )

      for (pair_index in seq_len(
        ncol(
          pld_pairs
        )
      )) {

        value_a <-
          pld_pairs[1, pair_index]

        value_b <-
          pld_pairs[2, pair_index]

        comparison <-
          compare_two_role_scenarios(
            fixed_data %>%
              filter(
                pld_days ==
                  value_a
              ),
            fixed_data %>%
              filter(
                pld_days ==
                  value_b
              )
          ) %>%
          mutate(
            parameter =
              "PLD",
            scenario_a =
              as.character(
                value_a
              ),
            scenario_b =
              as.character(
                value_b
              ),
            fixed_analysis_period =
              period_value,
            fixed_period_type =
              period_type_value,
            fixed_pld_days =
              NA_real_,
            fixed_precompetency_days =
              pc_value,
            .before =
              1
          )

        output[[output_index]] <-
          comparison

        output_index <-
          output_index +
          1L
      }
    }
  }

  # Pre-competency sensitivity
  #
  # Compare every pair of pre-competency values while holding analysis
  # period and PLD constant.

  for (period_value in period_values) {

    period_data <- metric_table %>%
      filter(
        analysis_period ==
          period_value
      )

    period_type_value <-
      unique(
        period_data$period_type
      )[1]

    pld_values <-
      sort(
        unique(
          period_data$pld_days
        )
      )

    for (pld_value in pld_values) {

      fixed_data <- period_data %>%
        filter(
          pld_days ==
            pld_value
        )

      pc_values <-
        sort(
          unique(
            fixed_data$precompetency_days
          )
        )

      if (length(pc_values) < 2) {
        next
      }

      pc_pairs <-
        combn(
          pc_values,
          2
        )

      for (pair_index in seq_len(
        ncol(
          pc_pairs
        )
      )) {

        value_a <-
          pc_pairs[1, pair_index]

        value_b <-
          pc_pairs[2, pair_index]

        comparison <-
          compare_two_role_scenarios(
            fixed_data %>%
              filter(
                precompetency_days ==
                  value_a
              ),
            fixed_data %>%
              filter(
                precompetency_days ==
                  value_b
              )
          ) %>%
          mutate(
            parameter =
              "precompetency",
            scenario_a =
              as.character(
                value_a
              ),
            scenario_b =
              as.character(
                value_b
              ),
            fixed_analysis_period =
              period_value,
            fixed_period_type =
              period_type_value,
            fixed_pld_days =
              pld_value,
            fixed_precompetency_days =
              NA_real_,
            .before =
              1
          )

        output[[output_index]] <-
          comparison

        output_index <-
          output_index +
          1L
      }
    }
  }

  # Release-timing sensitivity
  #
  # For the full study, seasonal comparisons are preferred because the
  # proposal groups monthly releases into austral seasons.
  #
  # If seasonal periods have not been created (as in the Taupo
  # pilot), the function falls back to comparing individual releases. This
  # keeps the pilot useful without pretending that June, July and August
  # represent different seasons.

  timing_period_type <-
    if (
      any(
        metric_table$period_type ==
          "season"
      )
    ) {

      "season"

    } else {

      "release"
    }

  timing_data <- metric_table %>%
    filter(
      period_type ==
        timing_period_type
    )

  pld_values <-
    sort(
      unique(
        timing_data$pld_days
      )
    )

  pc_values <-
    sort(
      unique(
        timing_data$precompetency_days
      )
    )

  for (pld_value in pld_values) {

    for (pc_value in pc_values) {

      fixed_data <- timing_data %>%
        filter(
          pld_days ==
            pld_value,
          precompetency_days ==
            pc_value
        )

      timing_values <-
        sort(
          unique(
            fixed_data$analysis_period
          )
        )

      if (length(timing_values) < 2) {
        next
      }

      timing_pairs <-
        combn(
          timing_values,
          2
        )

      for (pair_index in seq_len(
        ncol(
          timing_pairs
        )
      )) {

        value_a <-
          timing_pairs[1, pair_index]

        value_b <-
          timing_pairs[2, pair_index]

        comparison <-
          compare_two_role_scenarios(
            fixed_data %>%
              filter(
                analysis_period ==
                  value_a
              ),
            fixed_data %>%
              filter(
                analysis_period ==
                  value_b
              )
          ) %>%
          mutate(
            parameter =
              "release_timing",
            scenario_a =
              value_a,
            scenario_b =
              value_b,
            fixed_analysis_period =
              NA_character_,
            fixed_period_type =
              timing_period_type,
            fixed_pld_days =
              pld_value,
            fixed_precompetency_days =
              pc_value,
            .before =
              1
          )

        output[[output_index]] <-
          comparison

        output_index <-
          output_index +
          1L
      }
    }
  }

  bind_rows(
    output
  )
}

# Scenario pairs used by the sensitivity summaries
#
# The same "one factor at a time" comparison design is used for:
#
# - individual dispersal probabilities;
# - depth-exchange summaries; and
# - Infomap subnetwork composition.
#
# PLD comparisons hold pre-competency and analysis period constant.
# Pre-competency comparisons hold PLD and analysis period constant.
# Release-timing comparisons hold PLD and pre-competency constant.
#
# If seasonal periods are available, timing comparisons use seasons. If not,
# the Taupo pilot falls back to individual releases.

build_sensitivity_pair_index <- function(
    scenario_table) {

  scenarios <- scenario_table %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    mutate(
      scenario_key =
        make_scenario_name(
          analysis_period,
          pld_days,
          precompetency_days
        )
    )

  output <-
    list()

  output_index <-
    1L

  # PLD pairs

  period_pc_groups <- scenarios %>%
    distinct(
      analysis_period,
      period_type,
      precompetency_days
    )

  for (i in seq_len(
    nrow(
      period_pc_groups
    )
  )) {

    group_row <-
      period_pc_groups[i, ]

    dat <- scenarios %>%
      filter(
        analysis_period ==
          group_row$analysis_period,
        period_type ==
          group_row$period_type,
        precompetency_days ==
          group_row$precompetency_days
      ) %>%
      arrange(
        pld_days
      )

    if (nrow(dat) < 2) {
      next
    }

    pair_columns <-
      combn(
        seq_len(
          nrow(dat)
        ),
        2
      )

    for (j in seq_len(
      ncol(
        pair_columns
      )
    )) {

      a <-
        dat[
          pair_columns[1, j],
        ]

      b <-
        dat[
          pair_columns[2, j],
        ]

      output[[output_index]] <-
        tibble(
          parameter =
            "PLD",
          scenario_a =
            as.character(
              a$pld_days
            ),
          scenario_b =
            as.character(
              b$pld_days
            ),
          scenario_key_a =
            a$scenario_key,
          scenario_key_b =
            b$scenario_key,
          fixed_analysis_period =
            group_row$analysis_period,
          fixed_period_type =
            group_row$period_type,
          fixed_pld_days =
            NA_real_,
          fixed_precompetency_days =
            group_row$precompetency_days
        )

      output_index <-
        output_index +
        1L
    }
  }

  # Pre-competency pairs

  period_pld_groups <- scenarios %>%
    distinct(
      analysis_period,
      period_type,
      pld_days
    )

  for (i in seq_len(
    nrow(
      period_pld_groups
    )
  )) {

    group_row <-
      period_pld_groups[i, ]

    dat <- scenarios %>%
      filter(
        analysis_period ==
          group_row$analysis_period,
        period_type ==
          group_row$period_type,
        pld_days ==
          group_row$pld_days
      ) %>%
      arrange(
        precompetency_days
      )

    if (nrow(dat) < 2) {
      next
    }

    pair_columns <-
      combn(
        seq_len(
          nrow(dat)
        ),
        2
      )

    for (j in seq_len(
      ncol(
        pair_columns
      )
    )) {

      a <-
        dat[
          pair_columns[1, j],
        ]

      b <-
        dat[
          pair_columns[2, j],
        ]

      output[[output_index]] <-
        tibble(
          parameter =
            "precompetency",
          scenario_a =
            as.character(
              a$precompetency_days
            ),
          scenario_b =
            as.character(
              b$precompetency_days
            ),
          scenario_key_a =
            a$scenario_key,
          scenario_key_b =
            b$scenario_key,
          fixed_analysis_period =
            group_row$analysis_period,
          fixed_period_type =
            group_row$period_type,
          fixed_pld_days =
            group_row$pld_days,
          fixed_precompetency_days =
            NA_real_
        )

      output_index <-
        output_index +
        1L
    }
  }

  # Release-timing pairs
  #
  # Prefer seasonal periods for the full study. If no seasonal aggregation
  # exists, compare individual releases instead.

  timing_period_type <-
    if (
      any(
        scenarios$period_type ==
          "season"
      )
    ) {

      "season"

    } else {

      "release"
    }

  timing_scenarios <- scenarios %>%
    filter(
      period_type ==
        timing_period_type
    )

  pld_pc_groups <- timing_scenarios %>%
    distinct(
      pld_days,
      precompetency_days
    )

  for (i in seq_len(
    nrow(
      pld_pc_groups
    )
  )) {

    group_row <-
      pld_pc_groups[i, ]

    dat <- timing_scenarios %>%
      filter(
        pld_days ==
          group_row$pld_days,
        precompetency_days ==
          group_row$precompetency_days
      ) %>%
      arrange(
        analysis_period
      )

    if (nrow(dat) < 2) {
      next
    }

    pair_columns <-
      combn(
        seq_len(
          nrow(dat)
        ),
        2
      )

    for (j in seq_len(
      ncol(
        pair_columns
      )
    )) {

      a <-
        dat[
          pair_columns[1, j],
        ]

      b <-
        dat[
          pair_columns[2, j],
        ]

      output[[output_index]] <-
        tibble(
          parameter =
            "release_timing",
          scenario_a =
            a$analysis_period,
          scenario_b =
            b$analysis_period,
          scenario_key_a =
            a$scenario_key,
          scenario_key_b =
            b$scenario_key,
          fixed_analysis_period =
            NA_character_,
          fixed_period_type =
            timing_period_type,
          fixed_pld_days =
            group_row$pld_days,
          fixed_precompetency_days =
            group_row$precompetency_days
        )

      output_index <-
        output_index +
        1L
    }
  }

  bind_rows(
    output
  )
}

# Sensitivity of individual dispersal probabilities
#
# This directly implements the proposal statement that sensitivity analyses
# will examine changes in directed dispersal probabilities.
#
# All possible source-destination site-depth combinations are retained in
# connectivity, including zero values for released source nodes. This means
# the comparison can distinguish:
#
# - an edge that persists but changes strength;
# - a new non-zero edge appearing in scenario B; and
# - an existing non-zero edge disappearing in scenario B.

compare_two_connectivity_scenarios <- function(
    scenario_a,
    scenario_b) {

  joined <- scenario_a %>%
    select(
      source_node,
      source_depth_stratum,
      destination_node,
      destination_depth_stratum,
      connectivity_a =
        connectivity
    ) %>%
    inner_join(
      scenario_b %>%
        select(
          source_node,
          source_depth_stratum,
          destination_node,
          destination_depth_stratum,
          connectivity_b =
            connectivity
        ),
      by = c(
        "source_node",
        "source_depth_stratum",
        "destination_node",
        "destination_depth_stratum"
      )
    )

  if (nrow(
    joined
  ) == 0) {

    return(
      tibble(
        n_compared_edges =
          0L,
        n_nonzero_a =
          0L,
        n_nonzero_b =
          0L,
        n_edges_gained =
          0L,
        n_edges_lost =
          0L,
        spearman_rho =
          NA_real_,
        mean_absolute_change =
          NA_real_,
        max_absolute_change =
          NA_real_
      )
    )
  }

  joined %>%
    summarise(
      n_compared_edges =
        n(),
      n_nonzero_a =
        sum(
          connectivity_a > 0
        ),
      n_nonzero_b =
        sum(
          connectivity_b > 0
        ),
      n_edges_gained =
        sum(
          connectivity_a == 0 &
            connectivity_b > 0
        ),
      n_edges_lost =
        sum(
          connectivity_a > 0 &
            connectivity_b == 0
        ),
      spearman_rho =
        safe_spearman(
          connectivity_a,
          connectivity_b
        ),
      mean_absolute_change =
        mean(
          abs(
            connectivity_b -
              connectivity_a
          )
        ),
      max_absolute_change =
        max(
          abs(
            connectivity_b -
              connectivity_a
          )
        )
    )
}

calculate_connectivity_sensitivity <- function(
    connectivity) {

  scenario_table <- connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    )

  pair_index <-
    build_sensitivity_pair_index(
      scenario_table
    )

  connectivity_keyed <- connectivity %>%
    mutate(
      scenario_key =
        make_scenario_name(
          analysis_period,
          pld_days,
          precompetency_days
        )
    )

  output <- vector(
    "list",
    nrow(
      pair_index
    )
  )

  for (i in seq_len(
    nrow(
      pair_index
    )
  )) {

    pair <-
      pair_index[i, ]

    comparison <-
      compare_two_connectivity_scenarios(
        connectivity_keyed %>%
          filter(
            scenario_key ==
              pair$scenario_key_a
          ),
        connectivity_keyed %>%
          filter(
            scenario_key ==
              pair$scenario_key_b
          )
      )

    output[[i]] <-
      bind_cols(
        pair %>%
          select(
            -scenario_key_a,
            -scenario_key_b
          ),
        comparison
      )
  }

  bind_rows(
    output
  )
}

# Sensitivity of within- and cross-depth exchange
#
# The proposal distinguishes two related quantities:
#
# 1. within- and between-depth potential habitat-encounter frequency relative to all
#    particles released in a scenario; and
# 2. F_cross, the proportion of total source-normalised connectivity occurring
#    between depth strata.
#
# Keeping both makes the interpretation clear. The release-based fractions are
# descriptive event frequencies, whereas F_cross gives equal influence to each
# represented source node through the source-normalised connectivity matrix.

make_depth_exchange_scenario_summary <- function(
    within_between_depth,
    cross_depth_fraction) {

  within_summary <- within_between_depth %>%
    filter(
      depth_transition ==
        "within-depth"
    ) %>%
    select(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      within_depth_fraction =
        settlement_fraction
    )

  between_summary <- within_between_depth %>%
    filter(
      depth_transition ==
        "between-depth"
    ) %>%
    select(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      cross_depth_fraction_of_releases =
        settlement_fraction
    )

  within_summary %>%
    full_join(
      between_summary,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days"
      )
    ) %>%
    full_join(
      cross_depth_fraction %>%
        select(
          analysis_period,
          period_type,
          pld_days,
          precompetency_days,
          F_cross
        ),
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days"
      )
    )
}

calculate_depth_exchange_sensitivity <- function(
    within_between_depth,
    cross_depth_fraction) {

  scenario_summary <-
    make_depth_exchange_scenario_summary(
      within_between_depth =
        within_between_depth,
      cross_depth_fraction =
        cross_depth_fraction
    ) %>%
    mutate(
      scenario_key =
        make_scenario_name(
          analysis_period,
          pld_days,
          precompetency_days
        )
    )

  pair_index <-
    build_sensitivity_pair_index(
      scenario_summary
    )

  output <- vector(
    "list",
    nrow(
      pair_index
    )
  )

  for (i in seq_len(
    nrow(
      pair_index
    )
  )) {

    pair <-
      pair_index[i, ]

    a <- scenario_summary %>%
      filter(
        scenario_key ==
          pair$scenario_key_a
      )

    b <- scenario_summary %>%
      filter(
        scenario_key ==
          pair$scenario_key_b
      )

    output[[i]] <-
      pair %>%
      select(
        -scenario_key_a,
        -scenario_key_b
      ) %>%
      mutate(
        within_depth_fraction_a =
          a$within_depth_fraction[1],
        within_depth_fraction_b =
          b$within_depth_fraction[1],
        within_depth_fraction_difference =
          within_depth_fraction_b -
          within_depth_fraction_a,
        cross_depth_fraction_of_releases_a =
          a$cross_depth_fraction_of_releases[1],
        cross_depth_fraction_of_releases_b =
          b$cross_depth_fraction_of_releases[1],
        cross_depth_fraction_of_releases_difference =
          cross_depth_fraction_of_releases_b -
          cross_depth_fraction_of_releases_a,
        F_cross_a =
          a$F_cross[1],
        F_cross_b =
          b$F_cross[1],
        F_cross_difference =
          F_cross_b -
          F_cross_a
      )
  }

  bind_rows(
    output
  )
}

# Sensitivity of Infomap subnetwork composition
#
# Infomap module numbers cannot be compared directly between runs because
# labels are arbitrary. For example, "module 1" in one scenario is not
# necessarily the same ecological group as "module 1" in another.
#
# We therefore compare pairwise co-membership. For every pair of site-depth
# nodes represented in both scenarios, record whether the pair is grouped in
# the same Infomap module under each scenario.
#
# pairwise_membership_change_fraction = 0 means identical grouping.
# Larger values indicate progressively greater reorganisation.

compare_two_infomap_scenarios <- function(
    scenario_a,
    scenario_b) {

  joined <- scenario_a %>%
    select(
      node_id,
      depth_stratum,
      infomap_module_a =
        infomap_module
    ) %>%
    inner_join(
      scenario_b %>%
        select(
          node_id,
          depth_stratum,
          infomap_module_b =
            infomap_module
        ),
      by = c(
        "node_id",
        "depth_stratum"
      )
    )

  n_nodes <-
    nrow(
      joined
    )

  if (n_nodes < 2) {

    return(
      tibble(
        n_common_nodes =
          n_nodes,
        n_node_pairs =
          0L,
        n_changed_node_pairs =
          0L,
        pairwise_membership_change_fraction =
          NA_real_,
        pairwise_membership_agreement =
          NA_real_
      )
    )
  }

  node_pairs <-
    t(
      combn(
        seq_len(
          n_nodes
        ),
        2
      )
    )

  same_module_a <-
    joined$infomap_module_a[
      node_pairs[, 1]
    ] ==
    joined$infomap_module_a[
      node_pairs[, 2]
    ]

  same_module_b <-
    joined$infomap_module_b[
      node_pairs[, 1]
    ] ==
    joined$infomap_module_b[
      node_pairs[, 2]
    ]

  changed <-
    same_module_a !=
    same_module_b

  tibble(
    n_common_nodes =
      n_nodes,
    n_node_pairs =
      length(
        changed
      ),
    n_changed_node_pairs =
      sum(
        changed
      ),
    pairwise_membership_change_fraction =
      mean(
        changed
      ),
    pairwise_membership_agreement =
      1 -
      pairwise_membership_change_fraction
  )
}

calculate_infomap_sensitivity <- function(
    infomap_membership) {

  scenario_table <- infomap_membership %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    )

  pair_index <-
    build_sensitivity_pair_index(
      scenario_table
    )

  membership_keyed <- infomap_membership %>%
    mutate(
      scenario_key =
        make_scenario_name(
          analysis_period,
          pld_days,
          precompetency_days
        )
    )

  output <- vector(
    "list",
    nrow(
      pair_index
    )
  )

  for (i in seq_len(
    nrow(
      pair_index
    )
  )) {

    pair <-
      pair_index[i, ]

    comparison <-
      compare_two_infomap_scenarios(
        membership_keyed %>%
          filter(
            scenario_key ==
              pair$scenario_key_a
          ),
        membership_keyed %>%
          filter(
            scenario_key ==
              pair$scenario_key_b
          )
      )

    output[[i]] <-
      bind_cols(
        pair %>%
          select(
            -scenario_key_a,
            -scenario_key_b
          ),
        comparison
      )
  }

  bind_rows(
    output
  )
}

# Net vertical displacement
#
# delta_depth_m =
# destination benthic depth - source release depth
#
# Positive values indicate movement towards deeper habitat.
# Negative values indicate movement towards shallower habitat.

calculate_vertical_displacement <- function(
    settlements_period) {

  settlements_period %>%
    mutate(
      delta_depth_m =
        settlement_depth_m -
        source_depth_m,
      vertical_direction =
        case_when(
          delta_depth_m > 0 ~
            "deeper",
          delta_depth_m < 0 ~
            "shallower",
          TRUE ~
            "same exact depth"
        )
    )
}

# Cross-depth potential habitat-encounter classification
#
# Distinguish whether the qualifying potential habitat-encounter event occurs:
#
# a) within the same depth stratum and same seamount
# b) within the same depth stratum but at another seamount, or
# c) across depth strata on the same seamount
# d) across both seamount and depth stratum

classify_depth_transitions <- function(
    settlements_period) {

  settlements_period %>%
    mutate(
      cross_depth =
        source_depth_stratum !=
        destination_depth_stratum,
      same_seamount =
        source_node ==
        destination_node,
      transition_type =
        case_when(
          !cross_depth &
            same_seamount ~
            "same-seamount within-depth",
          !cross_depth &
            !same_seamount ~
            "between-seamount within-depth",
          cross_depth &
            same_seamount ~
            "same-seamount cross-depth",
          cross_depth &
            !same_seamount ~
            "between-seamount cross-depth"
        )
    )
}

# Summary of within- and cross-depth transition types
#
# The four categories distinguish:
#
# 1. same-seamount, within-depth potential habitat encounter;
# 2. between-seamount, within-depth connectivity;
# 3. same-seamount, cross-depth potential habitat encounter; and
# 4. between-seamount, cross-depth connectivity.
#
# Two fractions are retained:
#
# - fraction_of_releases: frequency relative to all released particles;
# - fraction_of_settlement_events: contribution among successful habitat encounters
#   potential habitat-encounter events only.

calculate_transition_type_summary <- function(
    settlements_period,
    releases_period) {

  classified_events <-
    classify_depth_transitions(
      settlements_period
    )

  total_releases <- releases_period %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    summarise(
      N_released =
        sum(
          N_released
        ),
      .groups =
        "drop"
    )

  total_events <- classified_events %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      name =
        "N_settlement_events"
    )

  transition_counts <- classified_events %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      transition_type,
      name =
        "N_transition"
    )

  transition_levels <- tibble(
    transition_type =
      c(
        "same-seamount within-depth",
        "between-seamount within-depth",
        "same-seamount cross-depth",
        "between-seamount cross-depth"
      )
  )

  dplyr::cross_join(
    total_releases,
    transition_levels
  ) %>%
    left_join(
      transition_counts,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "transition_type"
      )
    ) %>%
    left_join(
      total_events,
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days"
      )
    ) %>%
    mutate(
      N_transition =
        coalesce(
          N_transition,
          0L
        ),
      N_settlement_events =
        coalesce(
          N_settlement_events,
          0L
        ),
      fraction_of_releases =
        N_transition /
        N_released,
      fraction_of_settlement_events =
        if_else(
          N_settlement_events > 0,
          N_transition /
            N_settlement_events,
          NA_real_
        )
    )
}

# Overall cross-depth fraction
#
# F_cross is calculated from the source-normalised connectivity matrix:
#
# F_cross = sum(C for alpha != beta) / sum(C for all alpha,beta).
#
# This avoids weighting source site-depth habitats according to the number of
# particles simulated from them. Same-seamount and between-seamount cross-depth
# contributions are also retained separately for interpretation.

calculate_cross_depth_fraction <- function(
    connectivity) {

  connectivity %>%
    mutate(
      cross_depth =
        source_depth_stratum != destination_depth_stratum,
      same_seamount =
        source_node == destination_node
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    summarise(
      total_connectivity =
        sum(
          connectivity
        ),
      cross_depth_connectivity =
        sum(
          connectivity[cross_depth]
        ),
      same_seamount_cross_depth_connectivity =
        sum(
          connectivity[cross_depth & same_seamount]
        ),
      between_seamount_cross_depth_connectivity =
        sum(
          connectivity[cross_depth & !same_seamount]
        ),
      F_cross =
        if_else(
          total_connectivity > 0,
          cross_depth_connectivity / total_connectivity,
          NA_real_
        ),
      .groups = "drop"
    )
}

# Infomap connectivity subnetworks
#
# Infomap identifies groups of site-depth nodes linked by relatively
# directed network flow.
#
# Edge weights use connectivity probabilities directly:
# stronger dispersal probabilities represent stronger flow.

calculate_graph_infomap <- function(
    graph,
    nb_trials = infomap_trials,
    seed = infomap_seed) {

  n_nodes <-
    vcount(
      graph
    )

  # If there are no edges, each node is treated as its own module.

  if (ecount(graph) == 0) {

    membership <- seq_len(
      n_nodes
    )

    module_sizes <- rep(
      1L,
      n_nodes
    )

    return(
      list(
        summary =
          tibble(
            n_modules =
              n_nodes,
            largest_module_size =
              if_else(
                n_nodes > 0,
                1L,
                0L
              ),
            largest_module_proportion =
              if_else(
                n_nodes > 0,
                1 / n_nodes,
                NA_real_
              ),
            codelength =
              NA_real_,
            infomap_trials =
              nb_trials,
            infomap_seed =
              seed
          ),
        membership =
          tibble(
            graph_node =
              V(graph)$name,
            infomap_module =
              membership,
            module_size =
              module_sizes
          )
      )
    )
  }

  set.seed(
    seed
  )

  infomap <- igraph::cluster_infomap(
    graph,
    e.weights =
      E(graph)$connectivity,
    nb.trials =
      nb_trials
  )

  membership <-
    as.integer(
      igraph::membership(
        infomap
      )
    )

  module_sizes <- table(
    membership
  )

  largest_module_size <-
    max(
      as.integer(
        module_sizes
      )
    )

  summary_table <- tibble(
    n_modules =
      length(
        module_sizes
      ),
    largest_module_size =
      largest_module_size,
    largest_module_proportion =
      largest_module_size /
      n_nodes,
    codelength =
      infomap$codelength,
    infomap_trials =
      nb_trials,
    infomap_seed =
      seed
  )

  membership_table <- tibble(
    graph_node =
      V(graph)$name,
    infomap_module =
      membership
  ) %>%
    mutate(
      module_size =
        as.integer(
          module_sizes[
            as.character(
              infomap_module
            )
          ]
        )
    )

  list(
    summary =
      summary_table,
    membership =
      membership_table
  )
}

# Calculate Infomap subnetworks for every scenario

calculate_all_infomap <- function(
    connectivity,
    network_nodes,
    nb_trials = infomap_trials,
    seed = infomap_seed) {

  scenarios <- connectivity %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days
    )

  summary_output <- vector(
    "list",
    nrow(
      scenarios
    )
  )

  membership_output <- vector(
    "list",
    nrow(
      scenarios
    )
  )

  for (i in seq_len(
    nrow(
      scenarios
    )
  )) {

    scenario <-
      scenarios[i, ]

    graph <- build_multilayer_graph(
      connectivity =
        connectivity,
      network_nodes =
        network_nodes,
      analysis_period_value =
        scenario$analysis_period,
      pld_value =
        scenario$pld_days,
      precompetency_value =
        scenario$precompetency_days
    )

    infomap_result <-
      calculate_graph_infomap(
        graph = graph,
        nb_trials = nb_trials,
        seed = seed
      )

    summary_output[[i]] <-
      infomap_result$summary %>%
      mutate(
        analysis_period =
          scenario$analysis_period,
        period_type =
          scenario$period_type,
        pld_days =
          scenario$pld_days,
        precompetency_days =
          scenario$precompetency_days,
        .before =
          1
      )

    vertices <- igraph::as_data_frame(
      graph,
      what = "vertices"
    )

    membership_output[[i]] <-
      infomap_result$membership %>%
      left_join(
        vertices,
        by = c(
          "graph_node" =
            "name"
        )
      ) %>%
      transmute(
        analysis_period =
          scenario$analysis_period,
        period_type =
          scenario$period_type,
        pld_days =
          scenario$pld_days,
        precompetency_days =
          scenario$precompetency_days,
        node_id,
        name =
          seamount_name,
        chain,
        depth_stratum,
        infomap_module,
        module_size
      )
  }

  list(
    summary =
      bind_rows(
        summary_output
      ),
    membership =
      bind_rows(
        membership_output
      )
  )
}

# Cross-depth edge-removal comparison
#
# Construct a modified depth-resolved network in which all site-depth
# nodes are retained but connections between different depth strata
# are removed. Within-depth connectivity is unchanged.
#
# The comparison is deliberately restricted to the two analyses described
# in the Methods: betweenness centrality and Infomap subnetwork structure.
# Directed reachability and strongly connected components are not calculated.

remove_cross_depth_edges <- function(
    connectivity) {

  connectivity %>%
    mutate(
      N_settled =
        if_else(
          source_depth_stratum ==
            destination_depth_stratum,
          N_settled,
          0L
        ),
      connectivity =
        if_else(
          source_depth_stratum ==
            destination_depth_stratum,
          connectivity,
          0
        )
    )
}

# Compare node betweenness before and after removal of cross-depth edges.
#
# Positive differences indicate greater betweenness in the full network.
# Negative differences indicate that removing cross-depth edges causes
# alternative paths to increase the node's betweenness.

compare_cross_depth_betweenness <- function(
    full_betweenness,
    no_cross_depth_betweenness) {

  full_betweenness %>%
    rename(
      betweenness_full =
        betweenness,
      betweenness_normalised_full =
        betweenness_normalised
    ) %>%
    left_join(
      no_cross_depth_betweenness %>%
        select(
          analysis_period,
          period_type,
          pld_days,
          precompetency_days,
          node_id,
          depth_stratum,
          betweenness,
          betweenness_normalised
        ) %>%
        rename(
          betweenness_no_cross_depth =
            betweenness,
          betweenness_normalised_no_cross_depth =
            betweenness_normalised
        ),
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    ) %>%
    mutate(
      betweenness_rank_full =
        rank(
          -betweenness_normalised_full,
          ties.method =
            "average",
          na.last =
            "keep"
        ),
      betweenness_rank_no_cross_depth =
        rank(
          -betweenness_normalised_no_cross_depth,
          ties.method =
            "average",
          na.last =
            "keep"
        ),
      betweenness_difference =
        betweenness_full -
        betweenness_no_cross_depth,
      betweenness_normalised_difference =
        betweenness_normalised_full -
        betweenness_normalised_no_cross_depth,
      # Positive values mean the node has a numerically larger rank after
      # cross-depth edges are removed.
      betweenness_rank_difference =
        betweenness_rank_no_cross_depth -
        betweenness_rank_full
    ) %>%
    ungroup()
}

# Compare Infomap network-level summaries before and after
# cross-depth edge removal.
#
# Module identifiers themselves are arbitrary between separate Infomap
# runs, so module numbers are not compared directly.

compare_cross_depth_infomap <- function(
    full_infomap,
    no_cross_depth_infomap) {

  full_infomap %>%
    rename(
      n_modules_full =
        n_modules,
      largest_module_size_full =
        largest_module_size,
      largest_module_proportion_full =
        largest_module_proportion,
      codelength_full =
        codelength
    ) %>%
    left_join(
      no_cross_depth_infomap %>%
        rename(
          n_modules_no_cross_depth =
            n_modules,
          largest_module_size_no_cross_depth =
            largest_module_size,
          largest_module_proportion_no_cross_depth =
            largest_module_proportion,
          codelength_no_cross_depth =
            codelength
        ),
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days"
      )
    ) %>%
    mutate(
      module_number_difference =
        n_modules_no_cross_depth -
        n_modules_full,
      largest_module_size_difference =
        largest_module_size_full -
        largest_module_size_no_cross_depth,
      largest_module_proportion_difference =
        largest_module_proportion_full -
        largest_module_proportion_no_cross_depth,
      codelength_difference =
        codelength_full -
        codelength_no_cross_depth
    )
}

# Compare Infomap membership before and after cross-depth edge removal
#
# Infomap module labels are arbitrary between separate runs, so raw module
# numbers cannot be compared directly. Instead, comparison is based on
# pairwise co-membership:
#
# - for each pair of site-depth nodes, ask whether the two nodes belong to
#   the same module in the full network;
# - ask the same question after cross-depth edges are removed;
# - quantify the proportion of node pairs for which that relationship changes.
#
# This provides a label-invariant measure of change in subnetwork membership.

compare_cross_depth_infomap_membership <- function(
    full_membership,
    no_cross_depth_membership) {

  joined <- full_membership %>%
    select(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      node_id,
      depth_stratum,
      infomap_module_full =
        infomap_module
    ) %>%
    inner_join(
      no_cross_depth_membership %>%
        select(
          analysis_period,
          period_type,
          pld_days,
          precompetency_days,
          node_id,
          depth_stratum,
          infomap_module_no_cross_depth =
            infomap_module
        ),
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "node_id",
        "depth_stratum"
      )
    )

  scenarios <- joined %>%
    distinct(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days
    )

  output <- vector(
    "list",
    nrow(
      scenarios
    )
  )

  for (i in seq_len(
    nrow(
      scenarios
    )
  )) {

    scenario <-
      scenarios[i, ]

    dat <- joined %>%
      filter(
        analysis_period ==
          scenario$analysis_period,
        pld_days ==
          scenario$pld_days,
        precompetency_days ==
          scenario$precompetency_days
      )

    n_nodes <-
      nrow(
        dat
      )

    if (n_nodes < 2) {

      n_pairs <-
        0L

      n_changed_pairs <-
        0L

      pairwise_membership_change_fraction <-
        NA_real_

      pairwise_membership_agreement <-
        NA_real_

    } else {

      node_pairs <-
        t(
          combn(
            seq_len(
              n_nodes
            ),
            2
          )
        )

      same_module_full <-
        dat$infomap_module_full[
          node_pairs[, 1]
        ] ==
        dat$infomap_module_full[
          node_pairs[, 2]
        ]

      same_module_no_cross_depth <-
        dat$infomap_module_no_cross_depth[
          node_pairs[, 1]
        ] ==
        dat$infomap_module_no_cross_depth[
          node_pairs[, 2]
        ]

      changed <-
        same_module_full !=
        same_module_no_cross_depth

      n_pairs <-
        length(
          changed
        )

      n_changed_pairs <-
        sum(
          changed
        )

      pairwise_membership_change_fraction <-
        mean(
          changed
        )

      pairwise_membership_agreement <-
        1 -
        pairwise_membership_change_fraction
    }

    output[[i]] <- tibble(
      analysis_period =
        scenario$analysis_period,
      period_type =
        scenario$period_type,
      pld_days =
        scenario$pld_days,
      precompetency_days =
        scenario$precompetency_days,
      n_nodes =
        n_nodes,
      n_node_pairs =
        n_pairs,
      n_changed_node_pairs =
        n_changed_pairs,
      pairwise_membership_change_fraction =
        pairwise_membership_change_fraction,
      pairwise_membership_agreement =
        pairwise_membership_agreement
    )
  }

  bind_rows(
    output
  )
}

# Connectivity-threshold sensitivity
#
# The primary network retains all C > 0 edges. For sensitivity analysis, edges
# with C < tau are set to zero and node-level network roles and Infomap
# subnetworks are recalculated. Threshold values are configured near the top of
# the script rather than chosen automatically from the pilot.

apply_connectivity_threshold <- function(
    connectivity,
    tau) {

  connectivity %>%
    mutate(
      N_settled =
        if_else(
          connectivity >= tau,
          N_settled,
          0L
        ),
      connectivity =
        if_else(
          connectivity >= tau,
          connectivity,
          0
        )
    )
}

calculate_threshold_sensitivity <- function(
    connectivity,
    releases_period,
    network_nodes,
    thresholds = connectivity_thresholds) {

  thresholds <- sort(
    unique(
      thresholds[
        is.finite(thresholds) & thresholds > 0
      ]
    )
  )

  if (length(thresholds) == 0) {
    return(
      list(
        summary = tibble(),
        network_roles = tibble(),
        infomap_summary = tibble(),
        infomap_membership = tibble()
      )
    )
  }

  summary_output <- vector("list", length(thresholds))
  role_output <- vector("list", length(thresholds))
  infomap_summary_output <- vector("list", length(thresholds))
  infomap_membership_output <- vector("list", length(thresholds))

  for (i in seq_along(thresholds)) {
    tau <- thresholds[i]

    thresholded_connectivity <-
      apply_connectivity_threshold(
        connectivity = connectivity,
        tau = tau
      )

    summary_output[[i]] <- thresholded_connectivity %>%
      group_by(
        analysis_period,
        period_type,
        pld_days,
        precompetency_days
      ) %>%
      summarise(
        n_edges_retained = sum(connectivity > 0),
        total_connectivity_retained = sum(connectivity),
        .groups = "drop"
      ) %>%
      mutate(
        connectivity_threshold = tau,
        .before = 1
      )

    role_output[[i]] <-
      calculate_network_roles(
        connectivity = thresholded_connectivity,
        releases_period = releases_period,
        network_nodes = network_nodes
      ) %>%
      mutate(
        connectivity_threshold = tau,
        .before = 1
      )

    infomap_result <-
      calculate_all_infomap(
        connectivity = thresholded_connectivity,
        network_nodes = network_nodes,
        nb_trials = infomap_trials,
        seed = infomap_seed
      )

    infomap_summary_output[[i]] <-
      infomap_result$summary %>%
      mutate(
        connectivity_threshold = tau,
        .before = 1
      )

    infomap_membership_output[[i]] <-
      infomap_result$membership %>%
      mutate(
        connectivity_threshold = tau,
        .before = 1
      )
  }

  list(
    summary = bind_rows(summary_output),
    network_roles = bind_rows(role_output),
    infomap_summary = bind_rows(infomap_summary_output),
    infomap_membership = bind_rows(infomap_membership_output)
  )
}

# Run analysis

period_inputs <- prepare_analysis_periods(
  settlement_events =
    settlement_events,
  release_counts =
    release_counts,
  release_calendar =
    release_calendar
)

vertical_displacement_events <-
  calculate_vertical_displacement(
    period_inputs$settlements
  )

depth_transition_events <-
  classify_depth_transitions(
    period_inputs$settlements
  )

transition_type_summary <-
  calculate_transition_type_summary(
    settlements_period =
      period_inputs$settlements,
    releases_period =
      period_inputs$releases
  )

connectivity <- calculate_connectivity(
  settlements_period =
    period_inputs$settlements,
  releases_period =
    period_inputs$releases,
  network_nodes =
    network_nodes
)

cross_depth_fraction <-
  calculate_cross_depth_fraction(
    connectivity
  )

# Explicit sensitivity of directed dispersal probabilities
#
# This compares edge strengths while varying PLD, pre-competency or release
# timing one at a time. It also records edges that appear or disappear.

connectivity_sensitivity <-
  calculate_connectivity_sensitivity(
    connectivity
  )

infomap_results <-
  calculate_all_infomap(
    connectivity =
      connectivity,
    network_nodes =
      network_nodes,
    nb_trials =
      infomap_trials,
    seed =
      infomap_seed
  )

infomap_summary <-
  infomap_results$summary

infomap_membership <-
  infomap_results$membership

# Sensitivity of Infomap subnetwork composition across biological and
# temporal scenarios. Module labels are compared using pairwise co-membership
# rather than raw module numbers.

infomap_sensitivity <-
  calculate_infomap_sensitivity(
    infomap_membership
  )

# Network with cross-depth edges removed

connectivity_no_cross_depth <-
  remove_cross_depth_edges(
    connectivity
  )

# Betweenness comparison

betweenness_full <-
  calculate_all_betweenness(
    connectivity =
      connectivity,
    network_nodes =
      network_nodes
  )

betweenness_no_cross_depth <-
  calculate_all_betweenness(
    connectivity =
      connectivity_no_cross_depth,
    network_nodes =
      network_nodes
  )

cross_depth_betweenness_comparison <-
  compare_cross_depth_betweenness(
    full_betweenness =
      betweenness_full,
    no_cross_depth_betweenness =
      betweenness_no_cross_depth
  )

# Infomap comparison

infomap_no_cross_depth <-
  calculate_all_infomap(
    connectivity =
      connectivity_no_cross_depth,
    network_nodes =
      network_nodes
  )

infomap_summary_no_cross_depth <-
  infomap_no_cross_depth$summary

infomap_membership_no_cross_depth <-
  infomap_no_cross_depth$membership

cross_depth_infomap_comparison <-
  compare_cross_depth_infomap(
    full_infomap =
      infomap_summary,
    no_cross_depth_infomap =
      infomap_summary_no_cross_depth
  )

cross_depth_infomap_membership_comparison <-
  compare_cross_depth_infomap_membership(
    full_membership =
      infomap_membership,
    no_cross_depth_membership =
      infomap_membership_no_cross_depth
  )

settlement_summary <- calculate_settlement_summary(
  connectivity
)

directionality <- calculate_directionality(
  connectivity
)

network_roles <- calculate_network_roles(
  connectivity =
    connectivity,
  releases_period =
    period_inputs$releases,
  network_nodes =
    network_nodes
)

role_sensitivity <-
  calculate_role_sensitivity(
    network_roles
  )

threshold_sensitivity <-
  calculate_threshold_sensitivity(
    connectivity = connectivity,
    releases_period = period_inputs$releases,
    network_nodes = network_nodes,
    thresholds = connectivity_thresholds
  )

matrices <- build_all_matrices(
  connectivity =
    connectivity,
  network_nodes =
    network_nodes,
  candidate_nodes =
    candidate_nodes
)

depth_connectivity <- calculate_depth_connectivity(
  connectivity
)

within_between_depth <- calculate_within_between_depth(
  settlements_period =
    period_inputs$settlements,
  releases_period =
    period_inputs$releases
)

# Explicit sensitivity of within- and cross-depth exchange
#
# Compares both:
# - habitat-encounter fractions relative to all released particles; and
# - F_cross from source-normalised connectivity probabilities.

depth_exchange_sensitivity <-
  calculate_depth_exchange_sensitivity(
    within_between_depth =
      within_between_depth,
    cross_depth_fraction =
      cross_depth_fraction
  )

chain_connectivity <- calculate_chain_connectivity(
  connectivity =
    connectivity,
  candidate_nodes =
    candidate_nodes
)

chain_directionality <- calculate_chain_directionality(
  chain_connectivity
)

collapsed_connectivity <- calculate_collapsed_connectivity(
  connectivity =
    connectivity,
  releases_period =
    period_inputs$releases,
  candidate_nodes =
    candidate_nodes
)

collapsed_directionality <- calculate_collapsed_directionality(
  collapsed_connectivity
)

collapsed_roles <- calculate_collapsed_roles(
  collapsed_connectivity =
    collapsed_connectivity,
  candidate_nodes =
    candidate_nodes
)

depth_integrated_role_comparison <-
  compare_depth_resolved_collapsed_roles(
    network_roles =
      network_roles,
    collapsed_roles =
      collapsed_roles
  )

collapsed_matrices <- make_collapsed_matrices(
  collapsed_connectivity =
    collapsed_connectivity,
  candidate_nodes =
    candidate_nodes
)

# Collect analysis outputs
#
# Keep everything from 04 together in one object, including release-specific,
# seasonal (when available), and pooled results.
connectivity_analysis <- list(
  depth_levels =
    depth_levels,
  # If seasonal pooling is used in the full study, save the verified
  # release_id-to-month lookup with the analysis for reproducibility.
  release_calendar =
    release_calendar,
  analysis_settings =
    list(
      infomap_trials = infomap_trials,
      infomap_seed = infomap_seed,
      connectivity_thresholds = connectivity_thresholds,
      igraph_version = as.character(packageVersion("igraph")),
      R_version = R.version.string
    ),
  scenario_index =
    matrices$scenario_index,
  connectivity =
    connectivity,
  settlement_summary =
    settlement_summary,
  directionality =
    directionality,
  network_roles =
    network_roles,
  role_sensitivity =
    role_sensitivity,
  connectivity_sensitivity =
    connectivity_sensitivity,
  depth_exchange_sensitivity =
    depth_exchange_sensitivity,
  infomap_sensitivity =
    infomap_sensitivity,
  threshold_sensitivity_summary =
    threshold_sensitivity$summary,
  threshold_network_roles =
    threshold_sensitivity$network_roles,
  threshold_infomap_summary =
    threshold_sensitivity$infomap_summary,
  threshold_infomap_membership =
    threshold_sensitivity$infomap_membership,
  # Depth-dependent connectivity
  vertical_displacement_events =
    vertical_displacement_events,
  depth_transition_events =
    depth_transition_events,
  transition_type_summary =
    transition_type_summary,
  depth_connectivity =
    depth_connectivity,
  within_between_depth =
    within_between_depth,
  cross_depth_fraction =
    cross_depth_fraction,
  # Connectivity subnetworks
  infomap_summary =
    infomap_summary,
  infomap_membership =
    infomap_membership,
  # Cross-depth edge-removal comparison
  connectivity_no_cross_depth =
    connectivity_no_cross_depth,
  cross_depth_betweenness_comparison =
    cross_depth_betweenness_comparison,
  infomap_summary_no_cross_depth =
    infomap_summary_no_cross_depth,
  infomap_membership_no_cross_depth =
    infomap_membership_no_cross_depth,
  cross_depth_infomap_comparison =
    cross_depth_infomap_comparison,
  cross_depth_infomap_membership_comparison =
    cross_depth_infomap_membership_comparison,
  # Chain-level analysis
  chain_connectivity =
    chain_connectivity,
  chain_directionality =
    chain_directionality,
  # Matrices
  multilayer_matrices =
    matrices$multilayer,
  depth_block_matrices =
    matrices$depth_blocks,
  # Depth-collapsed comparison
  collapsed_connectivity =
    collapsed_connectivity,
  collapsed_directionality =
    collapsed_directionality,
  collapsed_roles =
    collapsed_roles,
  depth_integrated_role_comparison =
    depth_integrated_role_comparison,
  collapsed_matrices =
    collapsed_matrices
)

# Save the complete 04 analysis

saveRDS(
  connectivity_analysis,
  file.path(
    "Outputs",
    "connectivity_analysis.rds"
  )
)

# Non-zero depth-resolved connections

pooled_connectivity <- connectivity %>%
  filter(
    analysis_period == "pooled",
    connectivity > 0
  ) %>%
  mutate(
    source_depth_order =
      match(
        source_depth_stratum,
        depth_levels
      ),
    destination_depth_order =
      match(
        destination_depth_stratum,
        depth_levels
      )
  ) %>%
  arrange(
    pld_days,
    source_node,
    source_depth_order,
    destination_node,
    destination_depth_order
  ) %>%
  select(
    -source_depth_order,
    -destination_depth_order
  )

# Potential habitat-encounter probability by source-depth node

pooled_settlement_summary <- settlement_summary %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  mutate(
    depth_order =
      match(
        source_depth_stratum,
        depth_levels
      )
  ) %>%
  arrange(
    pld_days,
    source_node,
    depth_order
  ) %>%
  select(
    -depth_order
  )

# Forward-versus-reverse depth-resolved connections

pooled_directionality <- directionality %>%
  filter(
    analysis_period == "pooled",
    connectivity > 0 |
      reverse_connectivity > 0
  ) %>%
  arrange(
    pld_days,
    source_node,
    source_depth_stratum,
    destination_node,
    destination_depth_stratum
  )

# Multilayer network roles

pooled_network_roles <- network_roles %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  mutate(
    depth_order =
      match(
        depth_stratum,
        depth_levels
      )
  ) %>%
  arrange(
    pld_days,
    node_id,
    depth_order
  ) %>%
  select(
    -depth_order
  )

# Connectivity among depth strata

pooled_depth_connectivity <- depth_connectivity %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  arrange(
    pld_days,
    source_depth_stratum,
    destination_depth_stratum
  )

# Within-depth versus between-depth potential habitat encounter

pooled_within_between_depth <- within_between_depth %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  arrange(
    pld_days,
    depth_transition
  )

# Within- and between-chain connectivity

pooled_chain_connectivity <- chain_connectivity %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  arrange(
    pld_days,
    source_chain,
    destination_chain
  )

# Forward-versus-reverse chain comparison

pooled_chain_directionality <- chain_directionality %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  arrange(
    pld_days,
    source_chain,
    destination_chain
  )

# Between-chain connections only

pooled_between_chain <- chain_connectivity %>%
  filter(
    analysis_period == "pooled",
    chain_transition == "between-chain"
  ) %>%
  arrange(
    pld_days,
    source_chain,
    destination_chain
  )

# Depth-collapsed connections

pooled_collapsed_connectivity <- collapsed_connectivity %>%
  filter(
    analysis_period == "pooled",
    connectivity > 0
  ) %>%
  arrange(
    pld_days,
    source_node,
    destination_node
  )

# Depth-collapsed directionality

pooled_collapsed_directionality <- collapsed_directionality %>%
  filter(
    analysis_period == "pooled",
    connectivity > 0 |
      reverse_connectivity > 0
  ) %>%
  arrange(
    pld_days,
    source_node,
    destination_node
  )

# Depth-collapsed network roles

pooled_collapsed_roles <- collapsed_roles %>%
  filter(
    analysis_period == "pooled"
  ) %>%
  arrange(
    pld_days,
    node_id
  )

