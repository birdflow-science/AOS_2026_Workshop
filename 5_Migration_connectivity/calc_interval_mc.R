calc_interval_mc <- function(interval, bf) {

    # Initialize origin and target distributions
    origin_abun <- target_abun <- rep(0, n_active(bf))
    # Tally the origin and target cells
    target_table <- table(interval$i2)
    origin_table <- table(interval$i1)
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
      psi <- as.matrix(table(interval$i1,interval$i2))

  # normalize the matrix
  psi <- psi/rowSums(psi)

  # Align vectors to psi
  origin_abun <- origin_abun[as.numeric(rownames(psi))]
  target_abun <- target_abun[as.numeric(colnames(psi))]

  # multiply transition matrix and relative abundance (element-wise multiplication)
  # get a matrix where each cell (i, j) represents the Joint Probability: The probability that a bird exists and it starts at i AND it moves to j
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
  # origin_std %*% psi_abun: This takes the "starting neighborhood" information and pushes it through the transition links. It creates a matrix of "Where we came from, weighted by how many of us moved."
  # t(psi_abun) %*% ...: This completes the bridge. It results in a matrix that represents the covariance of the starting positions for birds that end up in various target positions.
  # * target_std: Finally, it checks if that "starting relationship" matches the "ending relationship" (target_std).
  MC = sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)

  return(MC)
}

# If you have the "transition matrix" with starting cells as rows, and ending cells as cols
matrix_to_interval <- function(matrix) {
  # Convert matrix to long format with i1, i2, and counts
  df <- as.data.frame(as.table(matrix))
  colnames(df) <- c("i1", "i2", "count")
  
  # Keep only rows with count > 0
  df <- df[df$count > 0, ]
  
  # Repeat each i1–i2 pair by the number of individuals
  interval_df <- df[rep(seq_len(nrow(df)), df$count), c("i1", "i2")]
  
  # Reset row names
  rownames(interval_df) <- NULL
  
  return(interval_df)
}

# calc_interval_mc <- function(interval, bf, weight = FALSE) {
#   
#   # --- 1. PREPARE PROBABILITIES ---
#   if (weight == TRUE) {
#     # If weights (effort) are provided, we treat the 'observed' data as the truth.
#     # We calculate the Joint Probability Matrix: P(start = i AND end = j)
#     w <- interval[["obs_weight"]]
#     
#     # Create a weighted joint distribution matrix
#     # This naturally handles zero-effort by giving those transitions 0 weight
#     joint_mat <- tapply(w, list(interval$i1, interval$i2), sum, default = 0)
#     
#     # Normalize the entire matrix so it sums to 1
#     total_w <- sum(joint_mat)
#     if (total_w == 0) return(NA) # Handle case where all effort/weights are zero
#     joint_mat <- joint_mat / total_w
#     
#     # Marginal distributions (the 'weighted' abundance at start and end)
#     origin_abun_v <- rowSums(joint_mat)
#     target_abun_v <- colSums(joint_mat)
#     
#     # The Transition Matrix (psi) is the conditional probability: P(end | start)
#     # Row-normalize the joint matrix
#     psi <- joint_mat / rowSums(joint_mat)
#     # Replace NaNs (from rowSums = 0) with 0
#     psi[is.na(psi)] <- 0
#     
#   } else {
#     # Standard unweighted version
#     joint_mat <- as.matrix(table(interval$i1, interval$i2))
#     joint_mat <- joint_mat / sum(joint_mat)
#     
#     origin_abun_v <- rowSums(joint_mat)
#     target_abun_v <- colSums(joint_mat)
#     
#     psi <- joint_mat / rowSums(joint_mat)
#     psi[is.na(psi)] <- 0
#   }
#   
#     # multiply transition matrix and relative abundance (element-wise multiplication)
#     # get a matrix where each cell (i, j) represents the Joint Probability: The probability that a bird exists and it starts at i AND it moves to j
#     psi_abun <- psi * origin_abun_v
#   
#   # --- 2. DISTANCE STANDARDIZATION ---
#   # Get distance matrix for all active cells
#   full_dist <- great_circle_distances(bf)
#   
#   # Helper to get standardized distance matrix based on an abundance vector
#   standardize_dist <- function(d, abun) {
#     # map the abundance vector back to all active cells for the product
#     # but the input 'abun' is already sized to the active cells of the BirdFlow object
#     # from the tapply/table logic above.
#     
#     # Map the marginal abundance to the matrix indices
#     idx <- as.numeric(names(abun))
#     
#     # mu (mean distance): sum of (distances * product of abundances)
#     # Using outer product to weight the distances between all pairs of pixels
#     mu <- sum(d[idx, idx] * (abun %*% t(abun)))
#     
#     # variance
#     var_d <- sum(((d[idx, idx] - mu)^2) * (abun %*% t(abun)))
#     sd_d <- sqrt(var_d)
#     
#     if (sd_d == 0) return(matrix(0, length(idx), length(idx)))
#     
#     return((d[idx, idx] - mu) / sd_d)
#   }
#   
#   # Standardize using the (weighted) marginals
#   # Rownames of psi = origin indices; Colnames = target indices
#   orig_idx <- as.numeric(rownames(psi))
#   targ_idx <- as.numeric(colnames(psi))
#   
#   # Extract relevant sub-vectors for standardization
#   origin_abun_sub <- origin_abun_v
#   target_abun_sub <- target_abun_v
#   
#   origin_std <- standardize_dist(full_dist, origin_abun_sub)
#   target_std <- standardize_dist(full_dist, target_abun_sub)
#   
#   # --- 3. FINAL MC CALCULATION ---
#   # Optimized Matrix Version:
#   # origin_std %*% psi_abun: This takes the "starting neighborhood" information and pushes it through the transition links. It creates a matrix of "Where we came from, weighted by how many of us moved."
#   # t(psi_abun) %*% ...: This completes the bridge. It results in a matrix that represents the covariance of the starting positions for birds that end up in various target positions.
#   # * target_std: Finally, it checks if that "starting relationship" matches the "ending relationship" (target_std).
#   MC <- sum(t(psi_abun) %*% origin_std %*% psi_abun * target_std)
#   
#   return(MC)
# }

