# Taupo Bank Analysis pt. 5 :)
# Depth-resolved connectivity and network analysis
#
# Builds connectivity from the first-settlement events calculated in 03.
# Depth is retained throughout so within- and between-depth connections
# can be analysed in the same multilayer network.
#
# Analysis outputs are saved together in connectivity_analysis.rds.
# Figures are made separately in 05.


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


# Inputs from settlement analysis

settlement_events <- readRDS(
  "settlement_events.rds"
)

release_counts <- readRDS(
  "release_counts.rds"
)

network_nodes <- readRDS(
  "network_nodes.rds"
)


# Site information

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


make_scenario_name <- function(
    analysis_period,
    pld_days,
    precompetency_days) {
  
  paste0(
    analysis_period,
    "__PLD",
    pld_days,
    "__PC",
    precompetency_days
  )
}


# Release periods and pooled data
#
# Keep individual releases for temporal comparisons, then pool the
# underlying particle counts for the overall result.
# Pooled connectivity is recalculated from counts rather than averaging
# release-specific probabilities.

prepare_analysis_periods <- function(
    settlement_events,
    release_counts) {
  
  release_settlements <- settlement_events %>%
    mutate(
      analysis_period =
        make_period_label(
          release_id
        ),
      
      period_type =
        "release"
    )
  
  
  pooled_settlements <- settlement_events %>%
    mutate(
      analysis_period =
        "pooled",
      
      period_type =
        "pooled"
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
  
  
  list(
    settlements =
      bind_rows(
        release_settlements,
        pooled_settlements
      ),
    
    releases =
      bind_rows(
        release_releases,
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
# seamount-depth node (i,alpha) that first settle at destination
# seamount-depth node (j,beta).
#
# Keep zero-connectivity destinations so a zero means the source was
# simulated but no particles settled there.

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


# Total settlement probability
#
# Because 03 keeps only the first valid settlement event per particle,
# summing outgoing connectivity gives the proportion settling somewhere
# in the network.

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
# Keep reverse connectivity as NA when the reverse source was not
# simulated rather than treating it as zero.

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
# Same-seamount settlement is kept separate as retention.

calculate_strength_retention_isolation <- function(
    connectivity,
    releases_period,
    network_nodes) {
  
  external_edges <- connectivity %>%
    filter(
      source_node !=
        destination_node
    )
  
  
  # External out-strength: export from this seamount-depth node
  # to other seamounts.
  
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
  
  
  # External in-strength: settlement into this seamount-depth node
  # from other seamounts.
  
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
  
  
  # Retention
  #
  # R_(i,alpha) = sum_beta P_(i,alpha -> i,beta)
  #
  # Includes settlement back onto the same seamount at either the same
  # or a different depth stratum.
  
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
  
  
  # Roles are reported for nodes from which particles were released.
  
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
      
      
      # Source-sink index
      #
      # SSI =
      # (s_OUT - s_IN) / (s_OUT + s_IN)
      #
      # Positive values indicate relatively greater export.
      # Negative values indicate relatively greater import.
      # Values near zero indicate more balanced external exchange.
      
      source_sink_index =
        if_else(
          s_out + s_in > 0,
          
          (s_out - s_in) /
            (s_out + s_in),
          
          NA_real_
        ),
      
      
      # External strength
      #
      # Lower external strength means weaker exchange with the wider
      # network and therefore greater relative isolation.
      
      s_external =
        s_out +
        s_in
    )
}


# Weighted shortest-path cost
#
# Network edges are probabilities, where larger values mean stronger
# connectivity. Convert them to costs for shortest-path calculations:
#
# cost = -log(P)
#
# This means the minimum-cost path is the path with the largest product
# of connectivity probabilities.

connectivity_to_cost <- function(
    connectivity) {
  
  -log(
    connectivity
  )
}


# Betweenness centrality
#
# Higher betweenness means a node lies on more minimum-cost directed
# paths, so it has greater potential importance as a stepping stone.

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
# Keep both within-depth and cross-depth edges.
#
# Exact self-loops are removed for betweenness, but same-seamount
# cross-depth links remain because they connect different depth nodes.

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


# Betweenness for every release-specific and pooled scenario.

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


# Combine strengths, retention, SSI, isolation and betweenness.

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
# 0 = source simulated but no settlement
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
# alpha = beta is within-depth; alpha != beta is cross-depth.

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


# Build multilayer matrices and all 16 depth blocks for every
# analysis scenario.

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
# P_(alpha -> beta) = N_(alpha -> beta) / N_alpha
#
# Collapse seamount identity but retain source and destination depth.

calculate_depth_connectivity <- function(
    settlements_period,
    releases_period) {
  
  releases_by_depth <- releases_period %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_depth_stratum
    ) %>%
    summarise(
      N_released =
        sum(
          N_released
        ),
      
      .groups =
        "drop"
    )
  
  
  settlements_by_depth <- settlements_period %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_depth_stratum,
      destination_depth_stratum,
      name = "N_settled"
    )
  
  
  dplyr::cross_join(
    releases_by_depth,
    
    tibble(
      destination_depth_stratum =
        depth_levels
    )
  ) %>%
    left_join(
      settlements_by_depth,
      
      by = c(
        "analysis_period",
        "period_type",
        "pld_days",
        "precompetency_days",
        "source_depth_stratum",
        "destination_depth_stratum"
      )
    ) %>%
    mutate(
      N_settled =
        coalesce(
          N_settled,
          0L
        ),
      
      depth_connectivity =
        N_settled /
        N_released,
      
      depth_transition =
        if_else(
          source_depth_stratum ==
            destination_depth_stratum,
          
          "within-depth",
          "between-depth"
        ),
      
      source_depth_stratum =
        factor(
          source_depth_stratum,
          levels =
            depth_levels
        ),
      
      destination_depth_stratum =
        factor(
          destination_depth_stratum,
          levels =
            depth_levels
        )
    ) %>%
    arrange(
      period_type,
      analysis_period,
      pld_days,
      source_depth_stratum,
      destination_depth_stratum
    )
}


# Overall within-depth versus between-depth settlement
#
# Both are calculated relative to all particles released in the
# scenario.

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
# Collapse seamount identity and depth within each chain, while keeping
# the direction of exchange between chains.

calculate_chain_connectivity <- function(
    settlements_period,
    releases_period,
    candidate_nodes) {
  
  chain_lookup <- candidate_nodes %>%
    select(
      node_id,
      chain
    ) %>%
    distinct()
  
  
  releases_by_chain <- releases_period %>%
    left_join(
      chain_lookup,
      
      by = c(
        "source_node" =
          "node_id"
      )
    ) %>%
    group_by(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      
      source_chain =
        chain
    ) %>%
    summarise(
      N_released =
        sum(
          N_released
        ),
      
      .groups =
        "drop"
    )
  
  
  settlements_by_chain <- settlements_period %>%
    left_join(
      chain_lookup %>%
        rename(
          source_chain =
            chain
        ),
      
      by = c(
        "source_node" =
          "node_id"
      )
    ) %>%
    left_join(
      chain_lookup %>%
        rename(
          destination_chain =
            chain
        ),
      
      by = c(
        "destination_node" =
          "node_id"
      )
    ) %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_chain,
      destination_chain,
      name = "N_settled"
    )
  
  
  destination_chains <- tibble(
    destination_chain =
      sort(
        unique(
          chain_lookup$chain
        )
      )
  )
  
  
  dplyr::cross_join(
    releases_by_chain,
    destination_chains
  ) %>%
    left_join(
      settlements_by_chain,
      
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
      N_settled =
        coalesce(
          N_settled,
          0L
        ),
      
      chain_connectivity =
        N_settled /
        N_released,
      
      chain_transition =
        if_else(
          source_chain ==
            destination_chain,
          
          "within-chain",
          "between-chain"
        )
    )
}


# Chain directionality
#
# Compares P_(A -> B) with P_(B -> A).
#
# The reverse value stays NA when the reverse source chain was not
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


# Depth-collapsed connectivity
#
# Collapse source and destination depth to give one connectivity value
# for each seamount pair.
#
# Recalculate this directly from particle counts rather than summing or
# averaging depth-specific probabilities.

calculate_collapsed_connectivity <- function(
    settlements_period,
    releases_period,
    candidate_nodes) {
  
  collapsed_releases <- releases_period %>%
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
      
      .groups =
        "drop"
    )
  
  
  collapsed_settlements <- settlements_period %>%
    count(
      analysis_period,
      period_type,
      pld_days,
      precompetency_days,
      source_node,
      destination_node,
      name = "N_settled"
    )
  
  
  destination_nodes <- candidate_nodes %>%
    distinct(
      destination_node =
        node_id
    )
  
  
  dplyr::cross_join(
    collapsed_releases,
    destination_nodes
  ) %>%
    left_join(
      collapsed_settlements,
      
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
      N_settled =
        coalesce(
          N_settled,
          0L
        ),
      
      connectivity =
        N_settled /
        N_released
    )
}


# Directionality in the depth-collapsed network.

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
# - betweenness uses -log(P) edge costs.

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
# Zero indicates a released source with no settlement at that
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


# Run analysis

period_inputs <- prepare_analysis_periods(
  settlement_events =
    settlement_events,
  
  release_counts =
    release_counts
)


connectivity <- calculate_connectivity(
  settlements_period =
    period_inputs$settlements,
  
  releases_period =
    period_inputs$releases,
  
  network_nodes =
    network_nodes
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


matrices <- build_all_matrices(
  connectivity =
    connectivity,
  
  network_nodes =
    network_nodes,
  
  candidate_nodes =
    candidate_nodes
)


depth_connectivity <- calculate_depth_connectivity(
  settlements_period =
    period_inputs$settlements,
  
  releases_period =
    period_inputs$releases
)


within_between_depth <- calculate_within_between_depth(
  settlements_period =
    period_inputs$settlements,
  
  releases_period =
    period_inputs$releases
)


chain_connectivity <- calculate_chain_connectivity(
  settlements_period =
    period_inputs$settlements,
  
  releases_period =
    period_inputs$releases,
  
  candidate_nodes =
    candidate_nodes
)


chain_directionality <- calculate_chain_directionality(
  chain_connectivity
)


collapsed_connectivity <- calculate_collapsed_connectivity(
  settlements_period =
    period_inputs$settlements,
  
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


collapsed_matrices <- make_collapsed_matrices(
  collapsed_connectivity =
    collapsed_connectivity,
  
  candidate_nodes =
    candidate_nodes
)


# Collect analysis outputs
#
# Keep everything from 04 together in one object, including both
# release-specific and pooled results.

connectivity_analysis <- list(
  depth_levels =
    depth_levels,
  
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
  
  depth_connectivity =
    depth_connectivity,
  
  within_between_depth =
    within_between_depth,
  
  chain_connectivity =
    chain_connectivity,
  
  chain_directionality =
    chain_directionality,
  
  multilayer_matrices =
    matrices$multilayer,
  
  depth_block_matrices =
    matrices$depth_blocks,
  
  collapsed_connectivity =
    collapsed_connectivity,
  
  collapsed_directionality =
    collapsed_directionality,
  
  collapsed_roles =
    collapsed_roles,
  
  collapsed_matrices =
    collapsed_matrices
)


# Save the complete 04 analysis

saveRDS(
  connectivity_analysis,
  "connectivity_analysis.rds"
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


# Settlement probability by source-depth node

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


# Within-depth versus between-depth settlement

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

# Check outputs

names(connectivity_analysis)

taupo_retention <- pooled_network_roles %>%
  filter(
    node_id == "T22"
  ) %>%
  select(
    pld_days,
    depth_stratum,
    N_released,
    s_out,
    s_in,
    retention,
    s_external,
    source_sink_index
  ) %>%
  arrange(
    pld_days,
    depth_stratum
  )

taupo_retention

pooled_connectivity

pooled_depth_connectivity

pooled_within_between_depth

pooled_chain_connectivity

pooled_collapsed_connectivity

View(taupo_retention)

View(pooled_connectivity)

View(pooled_depth_connectivity)

View(pooled_chain_connectivity)
