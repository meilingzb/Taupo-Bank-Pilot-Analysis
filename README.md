# Taupo Bank Pilot Analysis

R code supporting the Taupo Bank pilot analysis for an Honours research project investigating three-dimensional, depth-resolved potential coral larval connectivity across the Lord Howe and Tasmantid seamount chains in the Tasman Sea.

The pilot uses three-dimensional Lagrangian particle trajectories generated from the NASA ECCO LLC4320v2 circulation model. The analysis identifies potential encounters between competent particles and model-resolved benthic habitat, constructs directed depth-resolved connectivity networks, and evaluates how connectivity varies among depth strata and biological scenarios.

## Repository scope

This repository contains the analysis and figure-generation code used for the Taupo Bank pilot study presented in the Honours research proposal.

It does not contain the full LLC4320v2 model output, particle trajectory NetCDF files, GEBCO bathymetry, or other large source datasets. These files are excluded because of their size and, where applicable, external data-source requirements.

The pilot represents releases from Taupo Bank (T22) only. Network-role metrics calculated from the pilot should therefore be interpreted primarily as tests of the analytical workflow rather than as complete ecological estimates for the wider seamount network.

## Analysis workflow

The scripts are intended to be run in numerical order.

### `00_setup_and_inputs.R`

Loads and cleans the study-site inventory, standardises site information, loads GEBCO bathymetry and prepares common project inputs used by subsequent scripts.

### `01_data_exploration.R`

Performs exploratory checks of the study-site inventory and Taupo Bank pilot particle trajectories, including release depths, particle retention and depth-dependent displacement.

### `02_NETVOLC_delineation.R`

Implements an adapted form of the NETVOLC method of Euillades et al. (2013) to delineate the basal footprint of each seamount from GEBCO bathymetry. Candidate boundaries are checked for convergence across alternative spatial crop sizes and visually inspected against bathymetry.

### `03_settlement_criteria.R`

Constructs model-resolved benthic surfaces and identifies the first qualifying potential habitat encounter for each competent particle under alternative pelagic larval duration and pre-competency scenarios.

Potential encounters require the particle to:

- be competent and remain within its specified pelagic larval duration;
- occur within one area-equivalent local LLC4320v2 horizontal model-cell width of mapped benthic habitat;
- remain above the seabed; and
- occur within one local LLC4320v2 vertical model-layer thickness of the benthic surface.

Only the first qualifying encounter is retained for each particle and biological scenario.

### `04_network_analysis.R`

Converts potential habitat encounters into source-normalised directed connectivity probabilities and constructs depth-resolved connectivity networks.

Analyses include:

- external in- and out-degree;
- external in- and out-strength;
- source–sink index;
- same-seamount retention;
- relative isolation;
- weighted directed betweenness centrality;
- Infomap network structure;
- within- and cross-depth connectivity;
- depth-integrated connectivity;
- chain-level aggregation; and
- sensitivity analyses.

### Figure scripts

`05_study_site_figure.R`  
Generates the study-region and seamount-site map.

`06_depth_dispersal_figure.R`  
Generates the depth-dependent particle-displacement figure.

`07_geographic_connectivity_figure.R`  
Generates the geographic representation of connectivity for the reference pilot scenario.

`08_network_representation_figure.R`  
Compares multilayer, depth-independent and depth-integrated representations of the pilot network.

`09_biological_sensitivity_figure.R`  
Summarises habitat-encounter probability and the fraction of source-normalised connectivity occurring across depth strata under alternative PLD and pre-competency scenarios.

`Fig_S1_particle_dispersal.R`  
Generates the supplementary particle-dispersal maps.

`Fig_S2_connectivity_matrices.R`  
Generates connectivity matrices for all nine PLD × pre-competency scenarios.

## Biological scenarios

The pilot evaluates:

- pelagic larval durations of 30, 60 and 90 days; and
- pre-competency periods of 14, 21 and 28 days.

Particles are treated as neutrally buoyant and transported passively by the modelled circulation.

The analysis is species-agnostic and therefore estimates potential connectivity under defined biological assumptions rather than realised connectivity for a particular coral species.

## Depth representation

Habitat shallower than 200 m is excluded.

Connectivity is analysed within four depth strata:

- 200–700 m
- 700–1500 m
- 1500–2500 m
- ≥2500 m

Each network node represents a seamount–depth combination.

## Data requirements

Large source datasets are not stored in this repository.

The workflow requires, where relevant:

- GEBCO bathymetry;
- LLC4320v2 particle trajectories;
- LLC4320v2 horizontal-grid information;
- LLC4320v2 vertical-grid information; and
- study-site inventory files.

Local file paths used by the analysis assume a project structure containing `Data/` and `Outputs/` directories.

## R packages

The workflow uses packages including:

`dplyr`, `tidyr`, `tibble`, `sf`, `terra`, `ncdf4`, `igraph`, `ggplot2`, `ggraph`, `tidygraph`, `patchwork`, and `ggrepel`.

Additional package dependencies are identified at the beginning of the relevant scripts.

## Interpretation

Connectivity values represent modelled dispersal probabilities under the assumptions of the particle simulation and habitat-encounter criteria.

A qualifying event represents a potential opportunity for interaction with suitable benthic habitat. It should not be interpreted as confirmed larval settlement, recruitment or demographic connectivity.

Because the Taupo Bank pilot contains only a single source seamount, complete source–sink relationships and system-wide network roles cannot yet be inferred. These analyses are included to validate and demonstrate the workflow that will subsequently be applied across the full study region.

# AI declaration
Generative AI assistance (Claude, Anthropic) was used during development of this analysis workflow for code review against proposal methods, debugging, development and refinement of specific functions, and code organisation and documentation. Full details of AI use are provided in the Honours proposal.

## Author

Meiling Binney  
Bachelor of Advanced Science (Honours)  
The University of Queensland

Honours project: *Modelling three-dimensional deep-sea coral larval connectivity: Lord Howe and Tasmantid Seamount Chains*
