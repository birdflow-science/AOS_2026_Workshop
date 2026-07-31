#' Calculate migratory connectivity from BirdFlowR routes.
#'
#' This function calculates migratory connectivity based on sample of routes
#' objects sampled either from a BirdFlowR model or from observational track data.
#'
#' It calculates the migratory connectivity metric MC as defined in Cohen 2018.
#' MC represents an abundance-weighted correlation that is calculated between the origin
#' locations and target locations, taking into account all grid transitions for the
#' specified time period.
#'
#' The route implementation does not correct for spatial sampling inbalances of the
#' provided routes, i.e. the transition matrix is calculated directly from the route
#' transitions.
#'
#' When sampling a high number of routes from a BirdFlow model, the output of
#' [calc_route_mc()] will become asymptotically identical to the output of
#' [calc_birdflow_mc()].
#'
#' @param rts Output from [BirdFlowR::route()]
#' @param bf The BirdFlow model used to make `rts`
#' @param exact logical. Whether to match route time steps exactly to period
#' requested with [BirdFlowR::lookup_timestep_sequence()] (TRUE)
#' or use the route's closest available time steps (FALSE).
#' @param delta_steps The number of time steps (weeks) you allow for the nearest
#' timestep search when exact = FALSE. Default value is 2.
#' @inheritDotParams BirdFlowR::lookup_timestep_sequence -x
#' @return migratory connectivity estimated from the `rts` object.
#' @export
#' @seealso [calc_birdflow_mc()]
#' @examples
#' bf <- BirdFlowModels::amewoo
#' # generate 100 synthetic routes
#' rts <- route(bf, 100, season = "prebreeding")
#' # calculate MC across prebreeding season
#' calc_route_mc(rts, bf, season="prebreeding")
#' # calculate MC across a subset of weeks:
#' calc_route_mc(rts, bf, start=10, end=20)
#' # set exact to false to not enforce exact matches of route timestamps and
#' # requested start and end weeks (in this example end week 30 is after the last
#' # timestamp of the input routes):
#' calc_route_mc(rts, bf, start=10, end=30, exact=FALSE)
#' @references
#' Cohen EB, Hostetler JA, Hallworth MT, Rushing CS, Sillett TS, Marra PP.
#' Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution. 2018 Mar;9(3):513-24.
#' \doi{10.1111/2041-210X.12916}
calc_route_mc <- function(rts, bf, exact=TRUE, delta_steps = 2, ...) {
  # Using MigConnectivity::estPSI treating the routes as tracking data
  stopifnot(is.logical(exact))
  
  ts <- lookup_timestep_sequence(bf, ...)
  origin_t <- ts[1]
  target_t <- ts[length(ts)]

  if(exact){
    origin <- rts$data[rts$data$timestep %in% origin_t, ]
    target <- rts$data[rts$data$timestep %in% target_t,]
  }
  else{
    origin <- time_filter(rts, season = "postbreeding", delta_steps)$origin
    target <- time_filter(rts, season = "postbreeding", delta_steps)$target
  }

  # Initialize origin and target distributions
  origin_abun <- target_abun <- rep(0, n_active(bf))
  # Tally the origin and target cells
  target_table <- table(target$i)
  origin_table <- table(origin$i)
  # Populate origin and target distributions
  target_abun[as.numeric(names(target_table))]=target_table/sum(target_table)
  origin_abun[as.numeric(names(origin_table))]=origin_table/sum(origin_table)

  # Calculate distance matrices
  dist <- great_circle_distances(bf)  # all active cells

  # matrix product to have all combinations for two samples for origin_abun
  # note that this outer product remains normalized, because origin_abun is normalized
  origin_abun_prod <- origin_abun %*% t(origin_abun)
  # mu_D effectively equals weighted.mean(origin_dist, origin_abun_prod)
  mu_D <- sum(dist * origin_abun_prod)
  # SD in origin distance between any two given pixels
  sd_D <- sqrt(sum((dist - mu_D)^2 * origin_abun_prod))

  # same calculation, now for target_abun
  target_abun_prod <- target_abun %*% t(target_abun)
  mu_V <- sum(dist * target_abun_prod)
  sd_V <- sqrt(sum((dist - mu_V)^2 * target_abun_prod))

  # construct transition matrix
  psi <- as.matrix(table(origin$i,target$i))
  # normalize the matrix
  psi <- psi/rowSums(psi)

  origin_abun <- origin_abun[as.numeric(rownames(psi))]
  target_abun <- target_abun[as.numeric(colnames(psi))]

  # multiply transition matrix and relative abundance
  psi_abun <- psi * origin_abun

  # standardizing origin distance matrix
  origin_std <- (dist[as.numeric(rownames(psi)), as.numeric(rownames(psi))] - mu_D) / sd_D
  target_std <- (dist[as.numeric(colnames(psi)), as.numeric(colnames(psi))] - mu_V) / sd_V

  # reduce the distance matrix and abundance vectors
  dist <- dist[as.numeric(colnames(psi)), as.numeric(rownames(psi))]

  dim(psi_abun)
  dim(origin_std)
  dim(target_std)

  # calculate MC
  #MC=sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)
  MC=sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)

  return(MC)
}

get_circular_seq <- function(from, to, max_timestep = 52) {
  if (from <= to) {
    seq(from, to)
  } else {
    c(seq(from, max_timestep), seq(1, to))
  }
}

circular_week_distance <- function(x, target, max_week = 52) {
  pmin(abs(x - target), max_week - abs(x - target))
}

# Filter the routes that cover the "long jump"
time_filter = function(birdflow_routes_obj, bf, season){
  
  ts <- lookup_timestep_sequence(bf, season)
  origin_t <- ts[1]
  target_t <- ts[length(ts)]
  breeding_ts = lookup_timestep_sequence(bf, season = "breeding")
  nonbreeding_ts = lookup_timestep_sequence(bf, season = "nonbreeding")
  ids <- unique(birdflow_routes_obj$data$route_id)
  
  # allow a track to start as early as breeding season, and end as late as nonbreeding season
  # but find the nearest to both ends of postbreeding season if the track is long 
  
  if (season == "postbreeding"){
    start_range = breeding_ts
    end_range = nonbreeding_ts
  }
  if (season == "prebreeding"){
    start_range = nonbreeding_ts
    end_range = breeding_ts
  }
  
  table1 = birdflow_routes_obj$data |>
    distinct(route_id, route_type) |>
    count(route_type, name = "n_routes")
  
  # Subset each route so that:
  # 1. It has at least one detection in both start_range and end_range, in chronological order; When there are more than one chucks of start dates (spanning over a year), find the first chuck of start dates
  # 2. trim each track so that it begins with the last detection within start_range, and ends with the first detection within end_range.
  # 3. The resulting track only includes detections between those two points (inclusive).

  
  gap_tol <- 30  # tolerance in days between consecutive start dates (controls how big of a gap you allow between consecutive “start” detections before treating them as a new “chunk.”)
  
  birdflow_routes_filtered <- birdflow_routes_obj$data |>
    arrange(route_id, date) |>
    group_by(route_id) |>
    mutate(
      in_start_range = timestep %in% start_range,
      in_end_range   = timestep %in% end_range,
      length = n(),
      route_id = as.character(route_id)
    ) |>
    group_modify(~{
      df <- .x
      
      # If no start, return NA for this route
      if (!any(df$in_start_range)) {
        df$has_start <- FALSE
        df$has_end <- FALSE
        df$last_start_date <- as.Date(NA)
        df$first_end_date <- as.Date(NA)
        df$in_valid_segment <- FALSE
        return(df)
      }
      
      # Extract all start dates
      # Wrap max/min calls in suppressWarnings() to silence empty vector warnings
      suppressWarnings({
        start_dates <- sort(unique(df$date[df$in_start_range]))
      })
      
      # Identify chunks of consecutive start dates separated by > gap_tol days
      start_groups <- cumsum(c(1, diff(start_dates) > gap_tol))
      
      # Pick the *first chunk* of start dates
      first_chunk_dates <- start_dates[start_groups == min(start_groups)]
      
      # Use the max date from that chunk as the "last_start_date"
      last_start_date <- max(first_chunk_dates)
      
      # Find the first end date that occurs after the last start date
      valid_end_dates <- df$date[df$in_end_range & df$date > last_start_date]
      first_end_date <- if (length(valid_end_dates)) min(valid_end_dates) else as.Date(NA)
      
      # Mark segment membership
      df$has_start <- TRUE
      df$last_start_date <- last_start_date
      df$first_end_date <- first_end_date
      df$has_end <- !is.na(first_end_date)
      df$in_valid_segment <- df$date >= last_start_date & df$date <= first_end_date
      df
    }) |>
    filter(has_start & has_end) |>
    filter(in_valid_segment) |>
    ungroup()
  

origin <- birdflow_routes_filtered |>
  group_by(route_id) |>
  slice_min(date, n = 1, with_ties = FALSE) |>
  ungroup()

target <- birdflow_routes_filtered |>
  group_by(route_id) |>
  slice_max(date, n = 1, with_ties = FALSE) |>
  ungroup()

origin_target = merge(origin, target, by = c("route_id", "route_type"), suffixes = c("1","2")) |> 
  filter(date2>date1) |>  # this ensures the direction aligns with defined season
  filter(i1!=i2) # filter out non-moving individuals

print(paste("Original number of tracks:", length(ids), "   ","Tracks remaining:", nrow(origin_target)))

table2 = origin_target |>
  distinct(route_id, route_type) |>
  count(route_type, name = "n_routes")

return(list(
  origin = origin,
  target = target,
  origin_target = origin_target,
  routes_meta_before = table1,
  routes_meta_after = table2
))

}

make_origin_target_df <- function(sim_rts) {
  sim_origin = sim_rts$data |>
    group_by(route_id) |>
    arrange(date) |>
    filter(date == min(date)) |>
    ungroup()
  sim_target = sim_rts$data |>
    group_by(route_id) |>
    arrange(date) |>
    filter(date == max(date)) |>
    ungroup()
  sim_origin_target = merge(sim_origin, sim_target, by = c("route_id", "route_type"), suffixes = c("1","2"))
  return(sim_origin_target)
}
