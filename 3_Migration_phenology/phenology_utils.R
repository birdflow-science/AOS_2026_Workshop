library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)
library(scales)
library(tibble)
library(RColorBrewer)

birdflow_palette <- RColorBrewer::brewer.pal(3, "Set1")

# Identify the BirdFlow dates corresponding to a selected season.
get_season_dates <- function(bf, season = "prebreeding", season_buffer = 1) {
  if (!(season %in% c("prebreeding", "postbreeding", "breeding", "nonbreeding"))) {
    stop("season must be one of: prebreeding, postbreeding, breeding, nonbreeding")
  }
  
  season_timesteps <- BirdFlowR::lookup_season_timesteps(
    bf,
    season,
    season_buffer = season_buffer
  )
  
  bf$dates |>
    dplyr::filter(.data$timestep %in% season_timesteps) |>
    dplyr::pull(.data$date) |>
    as.Date()
}

# Simulate BirdFlow migration routes for a selected season.
# Optionally remove unrealistic tracks and replace them with new simulations.
simulate_birdflow_routes <- function(
    bf,
    season = "prebreeding",
    n_routes = 2000,
    season_buffer = 1,
    from_marginals = FALSE,
    filter_tracks = FALSE
) {
  dates <- get_season_dates(
    bf = bf,
    season = season,
    season_buffer = season_buffer
  )
  
  sim_routes <- BirdFlowR::route(
    bf,
    n = n_routes,
    from_marginals = from_marginals,
    season = season
  )
  
  if (filter_tracks) {
    if (!exists("filter_unrealistic_turns", mode = "function")) {
      stop("filter_tracks = TRUE requires filter_unrealistic_turns() to be loaded.")
    }
    
    sim_routes$data <- filter_unrealistic_turns(sim_routes$data)
    
    while (dplyr::n_distinct(sim_routes$data$route_id) < n_routes) {
      n_missing <- n_routes - dplyr::n_distinct(sim_routes$data$route_id)
      
      extra_routes <- BirdFlowR::route(
        bf,
        n = n_missing,
        from_marginals = from_marginals,
        season = season
      )
      
      extra_routes$data <- filter_unrealistic_turns(extra_routes$data)
      
      max_route_id <- max(sim_routes$data$route_id, na.rm = TRUE)
      
      extra_routes$data <- extra_routes$data |>
        dplyr::mutate(route_id = .data$route_id + max_route_id)
      
      sim_routes$data <- dplyr::bind_rows(
        sim_routes$data,
        extra_routes$data
      )
    }
  }
  
  sim_routes$data |>
    dplyr::mutate(date = as.Date(.data$date)) |>
    dplyr::filter(.data$date %in% dates)
}


# Summarize BirdFlow modeled abundance within a geographic box.
calculate_regional_abundance <- function(
    bf,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    season = "prebreeding",
    season_buffer = 1,
    from_marginals = FALSE
) {
  
  dates <- get_season_dates(
    bf = bf,
    season = season,
    season_buffer = season_buffer
  )
  
  # BirdFlow state distributions for each date
  st_dist <- BirdFlowR::get_distr(
    bf,
    which = dates,
    from_marginals = from_marginals
  )
  
  # Coordinates of BirdFlow spatial cells
  xy <- BirdFlowR::i_to_xy(
    seq_len(nrow(st_dist)),
    bf
  )
  
  latlon <- BirdFlowR::xy_to_latlon(
    bf = bf,
    x = xy$x,
    y = xy$y
  )
  
  latlon$cell_index <- seq_len(nrow(latlon))
  
  # Identify cells whose centers fall inside the selected region
  cells_inside <- latlon |>
    dplyr::filter(
      .data$lat >= lat_min,
      .data$lat <= lat_max,
      .data$lon >= lon_min,
      .data$lon <= lon_max
    )
  
  if (nrow(cells_inside) == 0) {
    stop("No BirdFlow cells fall inside the selected region.")
  }
  
  # Sum modeled probability mass within the region
  regional_dist <- st_dist[
    cells_inside$cell_index,
    ,
    drop = FALSE
  ]
  
  tibble::tibble(
    date = as.Date(dates),
    abundance = colSums(
      regional_dist,
      na.rm = TRUE
    ) / colSums(
      st_dist,
      na.rm = TRUE
    )
  )
}


# Plot modeled BirdFlow abundance within a region through time.
plot_regional_abundance <- function(
    regional_abundance,
    title = NULL,
    palette = birdflow_palette
) {
  
  ggplot2::ggplot(
    regional_abundance,
    ggplot2::aes(
      x = .data$date,
      y = .data$abundance
    )
  ) +
    ggplot2::geom_point(
      color = palette[3],
      size = 2.5,
      alpha = 0.7
    ) +
    ggplot2::geom_line(
      color = palette[3],
      linewidth = 1.2
    ) +
    ggplot2::scale_x_date(
      date_labels = "%b %d"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1)
    ) +
    ggplot2::labs(
      title = title,
      x = "Date",
      y = "Proportion of modeled population"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}


# Identify immigration and emigration events for simulated routes
# crossing the boundary of a selected geographic region.
calculate_weekly_movement <- function(
    sim_data,
    lon_min,
    lon_max,
    lat_min,
    lat_max
) {
  
  route_days <- sim_data |>
    dplyr::mutate(
      date = as.Date(.data$date),
      inside_flag =
        .data$lat >= lat_min &
        .data$lat <= lat_max &
        .data$lon >= lon_min &
        .data$lon <= lon_max
    ) |>
    dplyr::group_by(.data$route_id, .data$date) |>
    dplyr::summarise(
      inside = any(.data$inside_flag),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$route_id, .data$date) |>
    dplyr::group_by(.data$route_id) |>
    dplyr::mutate(
      prev_inside = dplyr::lag(
        .data$inside,
        default = dplyr::first(.data$inside)
      ),
      enters = !.data$prev_inside & .data$inside,
      leaves = .data$prev_inside & !.data$inside
    ) |>
    dplyr::ungroup()
  
  n_routes <- dplyr::n_distinct(route_days$route_id)
  
  route_days |>
    dplyr::group_by(.data$date) |>
    dplyr::summarise(
      n_enters = sum(.data$enters),
      n_leaves = sum(.data$leaves),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      date_mid = as.POSIXct(.data$date),
      enters_scaled = .data$n_enters / n_routes,
      leaves_scaled = .data$n_leaves / n_routes
    )
}

# Plot smoothed immigration and emigration phenology curves.
plot_immigration_emigration <- function(
    weekly_route_changes,
    title = NULL,
    palette = birdflow_palette
) {
  plot_data <- weekly_route_changes |>
    dplyr::filter(
      is.finite(.data$enters_scaled),
      is.finite(.data$leaves_scaled)
    )
  
  if (nrow(plot_data) < 4) {
    stop("Not enough route records to smooth immigration and emigration curves.")
  }
  
  ss_enters <- smooth.spline(
    x = as.numeric(plot_data$date_mid),
    y = plot_data$enters_scaled,
    spar = 0.15
  )
  
  ss_leaves <- smooth.spline(
    x = as.numeric(plot_data$date_mid),
    y = plot_data$leaves_scaled,
    spar = 0.15
  )
  
  grid_numeric <- seq(
    min(as.numeric(plot_data$date_mid)),
    max(as.numeric(plot_data$date_mid)),
    length.out = 1000
  )
  
  grid_df <- tibble::tibble(
    date_mid = as.POSIXct(grid_numeric, origin = "1970-01-01")
  ) |>
    dplyr::mutate(
      enters_smooth = predict(ss_enters, x = as.numeric(.data$date_mid))$y,
      leaves_smooth = predict(ss_leaves, x = as.numeric(.data$date_mid))$y,
      ymin = pmin(.data$enters_smooth, .data$leaves_smooth),
      ymax = pmax(.data$enters_smooth, .data$leaves_smooth),
      movement_type = factor(
        dplyr::if_else(.data$enters_smooth >= .data$leaves_smooth, "Immigration", "Emigration"),
        levels = c("Immigration", "Emigration")
      ),
      difference = abs(.data$enters_smooth - .data$leaves_smooth)
    ) |>
    dplyr::filter(.data$difference > 1e-6) |>
    dplyr::arrange(.data$date_mid) |>
    dplyr::mutate(
      changed = .data$movement_type != dplyr::lag(.data$movement_type, default = dplyr::first(.data$movement_type)),
      segment_id = cumsum(.data$changed)
    )
  
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = grid_df,
      ggplot2::aes(
        x = .data$date_mid,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$movement_type,
        group = .data$segment_id
      ),
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = grid_df,
      ggplot2::aes(x = .data$date_mid, y = .data$enters_smooth),
      color = palette[1],
      linewidth = 1.2
    ) +
    ggplot2::geom_line(
      data = grid_df,
      ggplot2::aes(x = .data$date_mid, y = .data$leaves_smooth),
      color = palette[2],
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = .data$date_mid, y = .data$enters_scaled),
      color = palette[1],
      size = 2.5,
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = .data$date_mid, y = .data$leaves_scaled),
      color = palette[2],
      size = 2.5,
      alpha = 0.7
    ) +
    ggplot2::scale_fill_manual(
      name = "BirdFlow movement",
      values = c(
        "Immigration" = palette[1],
        "Emigration" = palette[2]
      )
    ) +
    ggplot2::scale_x_datetime(date_labels = "%b %d") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(
      title = title,
      x = "Date",
      y = "Proportion of simulated routes"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}


# Run the complete regional migration phenology workflow.
run_phenology_case_study <- function(
    bf,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    season,
    n_routes = 1000,
    season_buffer = 1,
    filter_tracks = FALSE
) {
  
  sim_data <- simulate_birdflow_routes(
    bf = bf,
    season = season,
    n_routes = n_routes,
    season_buffer = season_buffer,
    filter_tracks = filter_tracks
  )
  
  weekly_movement <- calculate_weekly_movement(
    sim_data = sim_data,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max
  )
  
  list(
    sim_data = sim_data,
    weekly_movement = weekly_movement
  )
}
