# ekf_spring.R
# Extended Kalman Filter for a spring/mass/damper system
# Used to demonstrate the Dirichlet calibration framework
#
# TRUE SYSTEM: Duffing-type nonlinear spring/mass/damper
#   q'' = -(k/m)*q - (c/m)*q' - epsilon*q^3 + w_t
#   y_t = q_t + v_t
#
# MISSPECIFIED FORWARD MODEL (used in EKF):
#   q'' = -(k/m)*q - (c/m)*q'  (linear, ignores epsilon*q^3 term)
#
# The degree of misspecification is controlled by epsilon.
# When epsilon = 0 the EKF forward model is correct.
#
# State vector: x = c(q, q')  (position, velocity)
# Discretized via Euler-Maruyama with step dt.
#
# The EKF Jacobian F = d(f)/d(x) is computed analytically for the
# linear forward model, which is exact since the forward model is linear.
# If you switch to the nonlinear forward model in the EKF the Jacobian
# picks up the -3*epsilon*q^2 term in the (1,1) position.

# ==============================================================================
# 1. NUMERICAL INTEGRATION OF TRUE NONLINEAR SYSTEM
# ==============================================================================

#' Simulate the true nonlinear spring/mass/damper system
#'
#' Uses Euler-Maruyama discretization of the SDE:
#'   dq  = q' dt
#'   dq' = (-(k/m)*q - (c/m)*q' - epsilon*q^3) dt + sigma_w dW
#'
#' @param n_steps Integer. Number of time steps.
#' @param dt Numeric. Time step size.
#' @param k Numeric. Spring constant.
#' @param m Numeric. Mass.
#' @param c Numeric. Damping coefficient.
#' @param epsilon Numeric. Nonlinearity parameter. 0 = linear system.
#' @param sigma_w Numeric. Process noise standard deviation.
#' @param sigma_v Numeric. Observation noise standard deviation.
#' @param q0 Numeric. Initial position.
#' @param qdot0 Numeric. Initial velocity.
#' @param seed Optional integer seed.
#' @return List with:
#'   $states: (n_steps+1) x 2 matrix of true states (q, q')
#'   $obs: numeric vector of noisy observations
#'   $times: numeric vector of time points
simulate_spring <- function(n_steps  = 200,
                            dt       = 0.05,
                            k        = 2.0,
                            m        = 1.0,
                            c        = 0.3,
                            epsilon  = 0.0,
                            sigma_w  = 0.1,
                            sigma_v  = 0.2,
                            q0       = 1.0,
                            qdot0    = 0.0,
                            seed     = NULL) {
  if (!is.null(seed)) set.seed(seed)

  states       <- matrix(0, nrow = n_steps + 1, ncol = 2)
  states[1, ]  <- c(q0, qdot0)
  colnames(states) <- c("q", "qdot")

  # Process noise: enters through acceleration (second state component)
  w <- rnorm(n_steps, mean = 0, sd = sigma_w * sqrt(dt))

  for (i in seq_len(n_steps)) {
    q    <- states[i, 1]
    qdot <- states[i, 2]

    # Euler-Maruyama step
    q_new    <- q + qdot * dt
    qdot_new <- qdot + (-(k/m)*q - (c/m)*qdot - epsilon*q^3) * dt + w[i]

    states[i + 1, ] <- c(q_new, qdot_new)
  }

  # Observations: noisy position measurements
  obs <- states[, 1] + rnorm(n_steps + 1, mean = 0, sd = sigma_v)

  list(
    states = states,
    obs    = obs,
    times  = seq(0, n_steps * dt, by = dt)
  )
}

# ==============================================================================
# 2. LINEAR FORWARD MODEL (misspecified EKF)
# ==============================================================================

#' State transition function for the linear perfect spring (misspecified model)
#'
#' f(x) = F_mat %*% x  where F_mat is the linearized transition matrix
#' This ignores the epsilon*q^3 nonlinearity.
#'
#' @param x Numeric vector of length 2: c(q, qdot).
#' @param k Numeric. Spring constant.
#' @param m Numeric. Mass.
#' @param c Numeric. Damping coefficient.
#' @param dt Numeric. Time step.
#' @return Numeric vector of length 2: predicted next state.
f_linear <- function(x, k, m, c, dt) {
  # Euler discretization of linear system:
  # q_new    = q + qdot * dt
  # qdot_new = qdot + (-(k/m)*q - (c/m)*qdot) * dt
  q    <- x[1]
  qdot <- x[2]
  c(q + qdot * dt,
    qdot + (-(k/m)*q - (c/m)*qdot) * dt)
}

#' Jacobian of f_linear with respect to state x
#'
#' Since f_linear is linear in x, the Jacobian is exact (not an approximation).
#' F_jac = [[1, dt], [-(k/m)*dt, 1 - (c/m)*dt]]
#'
#' @param x Numeric vector (unused since Jacobian is constant for linear model).
#' @param k,m,c,dt As above.
#' @return 2x2 Jacobian matrix.
F_jacobian_linear <- function(x, k, m, c, dt) {
  matrix(c(1,          -(k/m)*dt,
           dt,          1 - (c/m)*dt),
         nrow = 2, ncol = 2, byrow = FALSE)
  # Column-major: col1 = c(1, -(k/m)*dt), col2 = c(dt, 1-(c/m)*dt)
}

# ==============================================================================
# 3. OBSERVATION MODEL
# ==============================================================================

#' Observation function: observe position only
#' h(x) = q = x[1]
h_obs <- function(x) x[1]

#' Jacobian of observation function
#' H = [1, 0]
H_jacobian <- function(x) matrix(c(1, 0), nrow = 1)

# ==============================================================================
# 4. EXTENDED KALMAN FILTER
# ==============================================================================

#' Run the Extended Kalman Filter
#'
#' Uses the linear perfect spring as the forward model (misspecified when
#' epsilon > 0 in the true system). Returns filtered state estimates,
#' one-step-ahead predictive means and variances, and log-likelihood.
#'
#' @param obs Numeric vector of observations y_1, ..., y_T.
#' @param k,m,c Numeric. System parameters (used in forward model).
#' @param dt Numeric. Time step.
#' @param sigma_w Numeric. Process noise standard deviation.
#' @param sigma_v Numeric. Observation noise standard deviation.
#' @param x0 Numeric vector of length 2. Initial state mean.
#' @param P0 2x2 matrix. Initial state covariance.
#' @return List with:
#'   $x_filt: (T+1) x 2 matrix of filtered state means
#'   $P_filt: list of (T+1) 2x2 filtered covariance matrices
#'   $x_pred: (T+1) x 2 matrix of one-step-ahead predicted state means
#'   $P_pred: list of (T+1) 2x2 predicted covariance matrices
#'   $y_pred: numeric vector of one-step-ahead predicted observations
#'   $S_pred: numeric vector of one-step-ahead predicted observation variances
#'   $log_lik: scalar log-likelihood
ekf <- function(obs,
                k       = 2.0,
                m       = 1.0,
                c       = 0.3,
                dt      = 0.05,
                sigma_w = 0.1,
                sigma_v = 0.2,
                x0      = c(1.0, 0.0),
                P0      = diag(c(0.1, 0.1)/100)) {

  T_obs <- length(obs)

  # Process noise covariance Q
  # Noise enters through velocity (second component) only
  Q <- matrix(c(0,            0,
                0, sigma_w^2 * dt),
              nrow = 2, ncol = 2)
  # Q <- matrix(c(sigma_w^2 * dt^3 / 3,  sigma_w^2 * dt^2 / 2,
  #               sigma_w^2 * dt^2 / 2,  sigma_w^2 * dt),
  #             nrow = 2, ncol = 2)
  # Observation noise variance R (scalar, position only)
  R <- sigma_v^2

  # Storage
  x_filt <- matrix(0, nrow = T_obs, ncol = 2)
  P_filt <- vector("list", T_obs)
  x_pred <- matrix(0, nrow = T_obs, ncol = 2)
  P_pred <- vector("list", T_obs)
  y_pred <- numeric(T_obs)
  S_pred <- numeric(T_obs)
  log_lik <- 0

  # Initialize with prior at t=0 used to predict t=1
  x_curr <- x0
  P_curr <- P0

  for (t in seq_len(T_obs)) {

    # --- Prediction step ---
    # Predict state at time t from filtered state at t-1
    x_p <- f_linear(x_curr, k, m, c, dt)
    F_j <- F_jacobian_linear(x_curr, k, m, c, dt)
    P_p <- F_j %*% P_curr %*% t(F_j) + Q

    # Predicted observation
    H_j <- H_jacobian(x_p)
    y_p <- h_obs(x_p)
    S_p <- as.numeric(H_j %*% P_p %*% t(H_j) + R)  # scalar

    # Store predictions
    x_pred[t, ] <- x_p
    P_pred[[t]]  <- P_p
    y_pred[t]    <- y_p
    S_pred[t]    <- S_p

    # --- Update step ---
    innovation <- obs[t] - y_p
    K          <- P_p %*% t(H_j) / S_p  # Kalman gain (2x1 vector)
    x_curr     <- x_p + as.numeric(K) * innovation
    P_curr     <- (diag(2) - K %*% H_j) %*% P_p

    # Store filtered estimates
    x_filt[t, ] <- x_curr
    P_filt[[t]]  <- P_curr

    # Log-likelihood contribution: N(y_t; y_pred_t, S_pred_t)
    log_lik <- log_lik +
      dnorm(obs[t], mean = y_p, sd = sqrt(S_p), log = TRUE)
  }

  list(
    x_filt  = x_filt,
    P_filt  = P_filt,
    x_pred  = x_pred,
    P_pred  = P_pred,
    y_pred  = y_pred,
    S_pred  = S_pred,
    log_lik = log_lik
  )
}

# ==============================================================================
# 5. MLE FOR k
# ==============================================================================

#' Find MLE of spring constant k via numerical optimization
#'
#' Maximizes the EKF log-likelihood over k with all other parameters fixed.
#'
#' @param obs Numeric vector of observations.
#' @param k_init Numeric. Initial value for k. Default 1.0.
#' @param ... Additional arguments passed to ekf().
#' @return List with $k_mle and $log_lik.
mle_k <- function(obs, k_init = 1.0, ...) {
  neg_ll <- function(log_k) {
    k_val <- exp(log_k)
    res   <- tryCatch(
      ekf(obs, k = k_val, ...),
      error = function(e) list(log_lik = -1e15)
    )
    -res$log_lik
  }

  opt <- optim(
    par     = log(k_init),
    fn      = neg_ll,
    method  = "Brent",
    lower   = log(0.01),
    upper   = log(20)
  )

  list(k_mle = exp(opt$par), log_lik = -opt$value)
}

# ==============================================================================
# 6. COVERAGE SIMULATION
# ==============================================================================

#' Run one simulation replicate and record coverage at multiple nominal levels
#'
#' @param nominal_levels Numeric vector of nominal coverage levels.
#' @param epsilon Numeric. Nonlinearity parameter.
#' @param n_steps,dt,k,m,c,sigma_w,sigma_v,q0,qdot0 System parameters.
#' @param use_mle Logical. If TRUE, estimate k by MLE. If FALSE use true k.
#' @return Logical matrix: rows = time steps, cols = nominal levels.
#'   TRUE if true position falls within the predictive interval.
one_replicate <- function(nominal_levels = c(0.5, 0.8, 0.9, 0.95, 0.99),
                          epsilon        = 0.0,
                          n_steps        = 200,
                          dt             = 0.05,
                          k              = 2.0,
                          m              = 1.0,
                          c              = 0.3,
                          sigma_w        = 0.1,
                          sigma_v        = 0.2,
                          q0             = 1.0,
                          qdot0          = 0.0,
                          use_mle        = FALSE) {

  # Simulate true system
  sim <- simulate_spring(n_steps = n_steps, dt = dt,
                         k = k, m = m, c = c,
                         epsilon = epsilon,
                         sigma_w = sigma_w, sigma_v = sigma_v,
                         q0 = q0, qdot0 = qdot0)

  # MLE or true k for the EKF forward model
  k_ekf <- if (use_mle) {
    tryCatch(
      mle_k(sim$obs[-1], k_init = k, m = m, c = c,
            dt = dt, sigma_w = sigma_w, sigma_v = sigma_v,
            x0 = c(sim$obs[1], 0), P0 = diag(c(sigma_v^2, 0.1)))$k_mle,
      error = function(e) k
    )
  } else {
    k
  }

  # Run EKF with (possibly misspecified) linear forward model
  ekf_res <- ekf(sim$obs[-1],   # observations at t=1,...,T
                 k       = k_ekf,
                 m       = m,
                 c       = c,
                 dt      = dt,
                 sigma_w = sigma_w,
                 sigma_v = sigma_v,
                 x0      = c(sim$obs[1], 0),
                 P0      = diag(c(sigma_v^2, sigma_w^2)))

  # True positions at t=1,...,T (one-step-ahead targets)
  true_q <- sim$states[-1, 1]

  # Check coverage at each nominal level
  # Predictive distribution: N(y_pred_t, S_pred_t)
  # Note S_pred includes observation noise so we compare against
  # observed y, but for coverage of the TRUE state we use the
  # filtered state predictive: N(x_pred_t[1], P_pred_t[1,1])
  coverage <- matrix(FALSE,
                     nrow = length(true_q),
                     ncol = length(nominal_levels))
  colnames(coverage) <- as.character(nominal_levels)
  # In one_replicate, after running ekf:

  for (j in seq_along(nominal_levels)) {
    alpha    <- 1 - nominal_levels[j]
    z        <- qnorm(1 - alpha / 2)
    pred_sd  <- sqrt(sapply(ekf_res$P_pred, function(P) P[1, 1]))
    lower    <- ekf_res$x_pred[, 1] - z * pred_sd
    upper    <- ekf_res$x_pred[, 1] + z * pred_sd
    coverage[, j] <- (true_q >= lower) & (true_q <= upper)
  }

  coverage
}

#' Run full coverage simulation study
#'
#' Runs n_sim replicates for each value of epsilon and computes empirical
#' coverage at each nominal level, with Beta posterior intervals.
#'
#' @param epsilon_vals Numeric vector of nonlinearity values to test.
#' @param nominal_levels Numeric vector of nominal coverage levels.
#' @param n_sim Integer. Number of simulation replicates per epsilon.
#' @param ... Additional arguments passed to one_replicate().
#' @return Data frame with columns: epsilon, nominal, k_covered, n_obs,
#'   coverage_mean, beta_lower, beta_upper.
coverage_simulation <- function(epsilon_vals   = c(0, 0.5, 1.0, 2.0),
                                nominal_levels = c(0.5, 0.8, 0.9, 0.95, 0.99),
                                n_sim          = 200,
                                ...) {
  results <- list()

  for (eps in epsilon_vals) {
    cat(sprintf("Running epsilon = %.2f ...\n", eps))

    # Accumulate coverage counts across replicates and time steps
    # Each replicate contributes n_steps observations per nominal level
    total_covered <- numeric(length(nominal_levels))
    total_obs     <- 0

    for (s in seq_len(n_sim)) {
      cov_mat      <- one_replicate(epsilon = eps,
                                    nominal_levels = nominal_levels, ...)
      total_covered <- total_covered + colSums(cov_mat)
      total_obs     <- total_obs + nrow(cov_mat)
    }

    # Beta posterior on coverage probability at each nominal level
    # Prior: Beta(1, 1) uniform
    for (j in seq_along(nominal_levels)) {
      k_cov  <- total_covered[j]
      n_obs  <- total_obs
      # Beta(k+1, n-k+1) posterior
      results[[length(results) + 1]] <- data.frame(
        epsilon      = eps,
        nominal      = nominal_levels[j],
        k_covered    = k_cov,
        n_obs        = n_obs,
        coverage_mean = (k_cov + 1) / (n_obs + 2),  # posterior mean
        beta_lower   = qbeta(0.025, k_cov + 1, n_obs - k_cov + 1),
        beta_upper   = qbeta(0.975, k_cov + 1, n_obs - k_cov + 1)
      )
    }
  }

  do.call(rbind, results)
}

# ==============================================================================
# 7. DIRICHLET COVERAGE MODEL
# ==============================================================================

#' Compute Dirichlet posterior for coverage gaps across nominal levels
#'
#' Given observed coverage counts at m ordered nominal levels, computes
#' the Dirichlet posterior on the gap probabilities
#' delta_i = p(alpha_i) - p(alpha_{i-1}).
#'
#' @param nominal_levels Numeric vector of nominal levels, sorted ascending.
#' @param k_covered Numeric vector of coverage counts at each level.
#' @param n_obs Integer. Total number of observations (same for all levels
#'   since intervals are nested).
#' @param gamma_prior Numeric. Dirichlet prior concentration. Default 1
#'   (uniform on simplex).
#' @return List with:
#'   $alpha: Dirichlet posterior parameters (length m+1)
#'   $gap_mean: posterior mean of gap probabilities
#'   $gap_lower: 2.5th percentile of each gap (via Beta marginals)
#'   $gap_upper: 97.5th percentile of each gap
dirichlet_coverage <- function(nominal_levels,
                               k_covered,
                               n_obs,
                               gamma_prior = 1) {

  m <- length(nominal_levels)
  stopifnot(length(k_covered) == m)
  stopifnot(all(diff(nominal_levels) > 0))

  # The m+1 multinomial categories are:
  # cat 0: outside all intervals  (count = n_obs - k_covered[m])
  # cat i: inside alpha_i but not alpha_{i-1} for i=1,...,m-1
  # cat m: inside alpha_1 (innermost interval)
  #
  # Note: k_covered[i] = sum of counts in categories i, i+1, ..., m
  # So counts in each gap category:
  n_gap <- numeric(m + 1)
  n_gap[1] <- n_obs - k_covered[m]          # outside all
  for (i in seq_len(m - 1)) {
    n_gap[i + 1] <- k_covered[m - i + 1] - k_covered[m - i]
  }
  n_gap[m + 1] <- k_covered[1]              # inside innermost

  # Dirichlet posterior parameters
  alpha_post <- n_gap + gamma_prior

  # Marginal of each gap is Beta(alpha_i, sum(alpha_{-i}))
  alpha_sum <- sum(alpha_post)
  gap_mean  <- alpha_post / alpha_sum
  gap_lower <- qbeta(0.025, alpha_post, alpha_sum - alpha_post)
  gap_upper <- qbeta(0.975, alpha_post, alpha_sum - alpha_post)

  list(
    alpha      = alpha_post,
    n_gap      = n_gap,
    gap_mean   = gap_mean,
    gap_lower  = gap_lower,
    gap_upper  = gap_upper,
    nominal_levels = nominal_levels
  )
}

# ==============================================================================
# 8. PLOTTING
# ==============================================================================

#' Plot reliability curves from coverage simulation
#'
#' @param df Data frame from coverage_simulation().
#' @param epsilon_vals Numeric vector. Which epsilon values to plot.
#' @param add_diagonal Logical. Add perfect calibration line. Default TRUE.
plot_reliability <- function(df,
                             epsilon_vals  = NULL,
                             add_diagonal  = TRUE) {

  if (is.null(epsilon_vals))
    epsilon_vals <- sort(unique(df$epsilon))

  cols <- hcl.colors(length(epsilon_vals), palette = "viridis")

  plot(NULL,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nominal coverage",
       ylab = "Empirical coverage",
       main = "Reliability curves: EKF predictive intervals")

  if (add_diagonal)
    abline(0, 1, lty = 2, col = "grey60", lwd = 1.5)

  for (i in seq_along(epsilon_vals)) {
    eps <- epsilon_vals[i]
    sub <- df[df$epsilon == eps, ]
    sub <- sub[order(sub$nominal), ]

    # Shaded credible band
    polygon(c(sub$nominal, rev(sub$nominal)),
            c(sub$beta_lower, rev(sub$beta_upper)),
            col  = adjustcolor(cols[i], alpha.f = 0.2),
            border = NA)

    # Coverage mean line
    lines(sub$nominal, sub$coverage_mean,
          col = cols[i], lwd = 2)
    points(sub$nominal, sub$coverage_mean,
           col = cols[i], pch = 16, cex = 0.8)
  }

  legend("topleft",
         legend = sprintf("epsilon = %.3f", epsilon_vals),
         col    = cols, lty = 1, lwd = 2, pch = 16,
         title  = "Nonlinearity", bty = "n")
}

#' Plot Dirichlet posterior on gap probabilities
#'
#' @param dir_result Output from dirichlet_coverage().
#' @param true_gaps Optional numeric vector of true gap probabilities
#'   for a perfectly calibrated model, for reference.
plot_dirichlet_gaps <- function(dir_result, true_gaps = NULL) {
  m     <- length(dir_result$nominal_levels)
  labs  <- c(
             sprintf("(%s, %s]",
                     c("0", as.character(dir_result$nominal_levels[-m])),
                     as.character(dir_result$nominal_levels)),
             "Outside all")#,sprintf("Inside %s", dir_result$nominal_levels[1]))
  # Reorder to match gap indexing
  x     <- seq_len(m + 1)

  plot(rev(x), dir_result$gap_mean,
       ylim  = range(c(dir_result$gap_lower, dir_result$gap_upper,
                       dir_result$gap_mean)) * c(0.9, 1.1),
       pch   = 16, col = "steelblue",
       xaxt  = "n",
       xlab  = "Coverage region",
       ylab  = "Gap probability",
       main  = "Dirichlet posterior on coverage gap probabilities")

  axis(1, at = x, labels = labs, las = 2, cex.axis = 0.7)

  segments(rev(x), dir_result$gap_lower,
           rev(x), dir_result$gap_upper,
           col = "steelblue", lwd = 2)

  if (!is.null(true_gaps)) {
    points(x, true_gaps, pch = 4, col = "red", cex = 1.2)
    legend("topleft",
           legend = c("Posterior mean", "95% interval", "True (perfect cal.)"),
           col    = c("steelblue", "steelblue", "red"),
           pch    = c(16, NA, 4), lty = c(NA, 1, NA), lwd = c(NA, 2, NA),
           bty    = "n")
  }
}

# ==============================================================================
# 9. EXAMPLE USAGE
# ==============================================================================

if (TRUE) {
  sigma_v=0.1
  sigma_w = 0.05
  # --- Single simulation to check the EKF is working ---
  set.seed(3457364)
  sim <- simulate_spring(n_steps = 200, dt = 0.05,
                         k = 2.0, m = 1.0, c = 0.1,
                         epsilon = 0,
                         sigma_w = sigma_w,
                         sigma_v = sigma_v)

  ekf_res <- ekf(sim$obs[-1], k = 2.0, m = 1.0, c = 0.1,
                 dt = 0.05, sigma_w = sigma_w, sigma_v = sigma_v,
                 x0 = c(sim$states[1, 1], 0),
                 P0 = diag(c(1e-6, sigma_w^2))  # nearly known initial position
                 )

  # Plot true vs filtered position
  times <- sim$times[-1]
  plot(times, sim$states[-1, 1], type = "l", col = "black",
       xlab = "Time", ylab = "Position",
       main = "True vs EKF filtered position (epsilon=1)")
  lines(times, ekf_res$x_filt[, 1], col = "blue", lty = 2)
  lines(times, ekf_res$x_pred[, 1] +
          2 * sqrt(sapply(ekf_res$P_pred, function(P) P[1,1])),
        col = "blue", lty = 3)
  lines(times, ekf_res$x_pred[, 1] -
          2 * sqrt(sapply(ekf_res$P_pred, function(P) P[1,1])),
        col = "blue", lty = 3)
  legend("topright",
         legend = c("True", "EKF filtered", "EKF pred +/- 2sd"),
         col = c("black", "blue", "blue"), lty = c(1, 2, 3))

  # --- MLE of k ---
  mle_res <- mle_k(sim$obs[-1], k_init = 1.5, m = 1.0, c = 0.1,
                   dt = 0.05, sigma_w = sigma_w, sigma_v = sigma_v,
                   x0 = c(sim$obs[1], 0))
  cat(sprintf("True k = 2.0, MLE k = %.4f\n", mle_res$k_mle))

  # --- Coverage simulation (small run for testing) ---
  df_cov <- coverage_simulation(
    epsilon_vals   = seq(0,0.2, by=0.04),
    nominal_levels = c(seq(0.1,0.95,by=0.15),0.95),
    n_sim          = 150,
    n_steps        = 100,
    dt             = 0.05,
    k = 2.0, m = 1.0, c = 0.1,
    sigma_w = sigma_w, sigma_v = sigma_v
  )

  # Reliability curves
  plot_reliability(df_cov)

  # --- Dirichlet analysis for a single epsilon ---
  sub_eps2 <- df_cov[df_cov$epsilon == 0.04,]#max(df_cov$epsilon), ]
  sub_eps2 <- sub_eps2[order(sub_eps2$nominal), ]

  dir_res <- dirichlet_coverage(
    nominal_levels = sub_eps2$nominal,
    k_covered      = sub_eps2$k_covered,
    n_obs          = sub_eps2$n_obs[1]
  )

  # Perfect calibration gap probabilities for reference
  noms    <- c(0, sub_eps2$nominal, 1)  # include 0 and 1 boundaries
  true_g  <- diff(noms)                 # uniform if perfectly calibrated

  plot_dirichlet_gaps(dir_res, true_gaps = true_g)
}
