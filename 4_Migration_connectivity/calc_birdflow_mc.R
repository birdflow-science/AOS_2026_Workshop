#' Calculate migratory connectivity (MC) for a BirdFlow model.
#'
#' This function calculates migratory connectivity based on the transitions
#' in a BirdFlowR model using the raster cells as regions.
#' It calculates the migratory connectivity metric MC as defined in Cohen 2018.
#' MC represents an abundance-weighted correlation that is calculated between the origin
#' locations and target locations, taking into account all grid transitions for the
#' specified time period.
#'
#' @param bf A BirdFlow model
#' @param from_marginals Use TRUE (the default) to sample from distributions derived from
#' the fitted model parameters stored in the marginals. Use FALSE to use distributions
#' derived directly from eBird Status and Trends when sampling starting locations.
#' @inheritDotParams BirdFlowR::lookup_timestep_sequence -x
#' @return The migratory connectivity of the BirdFlow model over the time period
#' indicated by `...`
#' @export
#' @import BirdFlowR
#' @seealso [calc_route_mc()]
#' @examples
#' # load a BirdFlow model
#' bf <- BirdFlowModels::amewoo
#' # calculate migratory connectivity across the spring migration season:
#' mc <- calc_birdflow_mc(bf, season = "prebreeding")
#' print(mc)
#' # calculate mc from week 10 to week 20:
#' calc_birdflow_mc(bf, start=10, end=20)
#' @references
#' Cohen EB, Hostetler JA, Hallworth MT, Rushing CS, Sillett TS, Marra PP.
#' Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution. 2018 Mar;9(3):513-24.
#' \doi{10.1111/2041-210X.12916}
#'
calc_birdflow_mc <- function(bf, from_marginals = TRUE, origin_t, target_t, ...) {
  # Figure out time
  ts <- lookup_timestep_sequence(bf, ...)
  origin_t <- ts[1]
  target_t <- ts[length(ts)]

  # Dynamic masks
  origin_dm <- get_dynamic_mask(bf, origin_t) # origin dynamic cells
  target_dm <- get_dynamic_mask(bf, target_t) # target dynamic cells

  # distance matrices
  dist <- great_circle_distances(bf)  # all active cells
  origin_dist <- dist[origin_dm, origin_dm] # origin dynamic
  target_dist <- dist[target_dm, target_dm] # target dynamic

  # Transition probabilities
  psi <- t(combine_transitions(bf, ...))

  # Origin and target relative abundance from marginals
  origin_abun <- get_distr(bf, origin_t, from_marginals = from_marginals)[origin_dm]
  if(from_marginals){
    target_abun <- get_distr(bf, target_t, from_marginals = from_marginals)[target_dm]
  } else{
    target_abun <- (origin_abun %*% psi)[1,]
  }

  # Double check dimensions
  stopifnot(isTRUE(all.equal(nrow(psi), sum(origin_dm))))
  stopifnot(isTRUE(all.equal(ncol(psi), sum(target_dm))))

  # matrix product to have all compbinations for two samples for origin_abun
  # note that this outer product remains normalized, because origin_abun is normalized
  origin_abun_prod <- origin_abun %*% t(origin_abun)
  # mu_D effectively equals weighted.mean(origin_dist, origin_abun_prod)
  mu_D <- sum(origin_dist * origin_abun_prod)
  # SD in origin distance between any two given pixels
  sd_D <- sqrt(sum((origin_dist - mu_D)^2 * origin_abun_prod))

  # same calculation, now for target_abun
  target_abun_prod <- target_abun %*% t(target_abun)
  mu_V <- sum(target_dist * target_abun_prod)
  sd_V <- sqrt(sum((target_dist - mu_V)^2 * target_abun_prod))

  # multiply transition matrix and relative abundance
  psi_abun <- psi * origin_abun

  # standardizing origin distance matrix
  origin_std <- (origin_dist - mu_D) / sd_D
  target_std <- (target_dist - mu_V) / sd_V

  # calculate MC
  MC=sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)

  return(MC)
}
