# ============================================================================
# engine_cavi.R
#
# PASSAGE: Layer 1 / Engine choice E3
#   Sparse-loading kernel-ordered LMC, fit by Coordinate-Ascent Variational
#   Inference (CAVI) with spike-and-slab priors on A.
#
# Model
#   Y_ij = mu_j + sum_k A_jk V_ik + epsilon_ij,
#   v_k ~ N(0, sigma_k^2 R_k(phi_k))       (Vecchia-approximated)
#   A_jk | gamma_jk = 1 ~ N(0, nu_k),
#   A_jk | gamma_jk = 0 = 0,
#   gamma_jk ~ Bern(pi_k),
#   tau_j^2 ~ InvGam(a_tau, b_tau),
#   sigma_k^2 ~ InvGam(a_sigma, b_sigma).
#
# Identifiability via kernel-scale ordering:
#   phi_1 = exp(alpha_1),   phi_k = exp(alpha_1) + sum_{l=2..k} exp(alpha_l).
#
# Variational family (mean-field):
#   q(mu, V, A, gamma, D, sigma^2)
#     = prod_j q(mu_j) prod_k q(v_k) prod_{j,k} q(A_jk, gamma_jk)
#       prod_j q(tau_j^2) prod_k q(sigma_k^2)
#
# Returns the standardized engine output (same shape as engine_pca):
#   A_hat, theta_hat, D_hat, V_scores, Qtilde_list, vecchia, ortho, ...
#
# Dependencies: Matrix, GpGp, ucminf, engine_pca.R (helpers reused)
# Version:      0.1.0 (MVP)
# ============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(GpGp)
  library(ucminf)
})


# ============================================================================
# Initialisation from PCA two-stage
# ============================================================================
.cavi_initialise <- function(Y, locs, K, m, ordering, smoothness, verbose) {
  if (verbose) cat("  CAVI initialisation: PCA two-stage warm start\n")
  init <- fit_engine_pca(
    Y = Y, locs = locs, K = K, m = m, ordering = ordering,
    smoothness = smoothness, D_orthogonalize = FALSE, verbose = FALSE)

  # The CAVI maintains q(v_k) as a Gaussian with mean m_v_k and precision matrix
  # Sigma_inv_v_k.  Initial Sigma_v_k diagonal is set from the GP fit
  # marginal variances; m_v_k initialised to the PCA factor scores in ordered
  # indexing.
  N <- init$N; G <- init$G

  ord <- init$vecchia$ord

  state <- list(
    # Variational parameters
    m_mu = init$mu_hat,                       # G-vector
    s_mu2 = rep(1e-2, G),                     # G-vector
    m_A = init$A_hat,                         # G x K (slab mean conditional on gamma=1)
    s_A2 = matrix(1e-2, G, K),                # G x K (slab var conditional on gamma=1)
    omega = matrix(0.5, G, K),                # G x K (inclusion probs)
    m_v = matrix(0, N, K),                    # N x K (in Vecchia ordering)
    diag_Sigma_v = matrix(1.0, N, K),         # N x K
    a_tau_hat = rep(NA_real_, G),
    b_tau_hat = rep(NA_real_, G),
    a_sigma_hat = rep(NA_real_, K),
    b_sigma_hat = rep(NA_real_, K),
    # Hyperparameters (EB-updated)
    nu = pmax(apply(init$A_hat, 2L, var), 1e-3),
    pi = rep(0.3, K),
    # Kernel parameters
    alpha = numeric(K),                       # cumulative-exp param
    phi = init$theta_hat[, "phi"],
    sigma2 = init$theta_hat[, "sigma2"],
    # Cached precision matrices and conditional regression info per factor
    Qtilde = init$Qtilde_list,
    vecchia_aux = vector("list", K),          # built from build_vecchia_precision()
    # Constants
    sigma_mu2 = 1e6,
    a_tau = 1e-3, b_tau = 1e-3,
    a_sigma = 1e-3, b_sigma = 1e-3
  )

  # Reorder PCA scores to Vecchia ordering
  state$m_v <- init$V_scores[ord, , drop = FALSE]

  # Initialise omega from standardised |A| (genes with strong loadings get
  # high inclusion prob)
  tau_init <- pmax(init$D_hat, 1e-6)
  std_load <- abs(init$A_hat) / sqrt(tau_init)
  state$omega <- plogis(std_load - 1.5)
  state$omega <- pmin(pmax(state$omega, 1e-4), 1 - 1e-4)

  # Initial tau2: keep PCA estimates
  state$a_tau_hat <- rep(state$a_tau + nrow(Y) / 2, G)
  state$b_tau_hat <- state$b_tau + 0.5 * init$D_hat *
                     (rep(state$a_tau_hat[1], G) - state$a_tau) * 2 / nrow(Y) * nrow(Y)
  # Simplified: just set tau_inv = 1/D_hat directly
  state$b_tau_hat <- init$D_hat * state$a_tau_hat

  # Initial sigma2 similarly
  state$a_sigma_hat <- rep(state$a_sigma + nrow(Y) / 2, K)
  state$b_sigma_hat <- state$sigma2 * state$a_sigma_hat

  # Enforce kernel ordering by sorting phi (and the rest accordingly)
  phi_order <- order(state$phi)
  state$phi      <- state$phi[phi_order]
  state$sigma2   <- state$sigma2[phi_order]
  state$nu       <- state$nu[phi_order]
  state$pi       <- state$pi[phi_order]
  state$m_A      <- state$m_A[, phi_order, drop = FALSE]
  state$s_A2     <- state$s_A2[, phi_order, drop = FALSE]
  state$omega    <- state$omega[, phi_order, drop = FALSE]
  state$m_v      <- state$m_v[, phi_order, drop = FALSE]
  state$Qtilde   <- state$Qtilde[phi_order]
  state$a_sigma_hat <- state$a_sigma_hat[phi_order]
  state$b_sigma_hat <- state$b_sigma_hat[phi_order]
  # Cumulative-exp parameterisation
  state$alpha[1] <- log(max(state$phi[1], 1e-6))
  if (K >= 2L) {
    for (k in 2:K) {
      state$alpha[k] <- log(max(state$phi[k] - state$phi[k - 1], 1e-6))
    }
  }

  list(state = state, engine_pca_init = init)
}


# ============================================================================
# Helper: rebuild Qtilde_k from current phi_k
# ============================================================================
.cavi_rebuild_Qtilde <- function(state, locs_ord, NNarray, smoothness) {
  K <- length(state$phi)
  for (k in seq_len(K)) {
    qk <- .build_vecchia_precision(
      locs_ord = locs_ord, NNarray = NNarray,
      phi = state$phi[k], smoothness = smoothness)
    state$Qtilde[[k]] <- qk$Q_tilde
    state$vecchia_aux[[k]] <- qk
  }
  state
}


# ============================================================================
# Per-cycle CAVI updates
# ============================================================================

# ---- Update mu_j -------------------------------------------------------------
.update_mu <- function(state, Y_ord, tau_inv) {
  N <- nrow(Y_ord)
  EA <- state$omega * state$m_A
  fitted_spatial <- state$m_v %*% t(EA)   # N x G
  resid <- Y_ord - fitted_spatial          # N x G
  prec_mu <- N * tau_inv + 1 / state$sigma_mu2
  state$s_mu2 <- 1 / prec_mu
  state$m_mu  <- state$s_mu2 * (tau_inv * colSums(resid))
  state
}

# ---- Update v_k -------------------------------------------------------------
.update_v <- function(state, Y_ord, tau_inv) {
  N <- nrow(Y_ord); K <- ncol(state$m_v)
  EA  <- state$omega * state$m_A                              # G x K
  EAA <- state$omega * (state$m_A^2 + state$s_A2)             # G x K
  mu_row <- matrix(state$m_mu, nrow = N, ncol = ncol(Y_ord), byrow = TRUE)

  for (k in seq_len(K)) {
    # Partial residual: Y - mu - sum_{l != k} A_l V_l
    EA_minus_k <- EA;  EA_minus_k[, k] <- 0
    Y_partial <- Y_ord - mu_row - state$m_v %*% t(EA_minus_k)

    # Linear coefficient c_ik = sum_j tau_inv_j * EA_jk * Y_partial_ij
    c_k <- as.numeric(Y_partial %*% (EA[, k] * tau_inv))      # length N

    # Quadratic coefficient d_k (scalar)
    d_k <- sum(EAA[, k] * tau_inv)

    # Variational precision
    sigma_inv_k <- state$a_sigma_hat[k] / state$b_sigma_hat[k]
    Sigma_inv_v_k <- sigma_inv_k * state$Qtilde[[k]] + Diagonal(N, max(d_k, 1e-10))

    # Sparse Cholesky and solve
    chol_factor <- tryCatch(
      Cholesky(Sigma_inv_v_k, perm = FALSE, LDL = FALSE, super = FALSE),
      error = function(e) {
        warning("Cholesky failed for v_k update; using diagonal ridge")
        Cholesky(Sigma_inv_v_k + Diagonal(N, 1e-6), perm = FALSE)
      })

    state$m_v[, k] <- as.numeric(solve(chol_factor, c_k))

    # Diagonal of Sigma_v_k via dense inverse (v1 fallback)
    # TODO v2: replace with Takahashi selected inversion via sparseinv pkg
    # For moderate N this is acceptable; for N > ~10000 it's prohibitive.
    Sigma_v_k_diag <- tryCatch({
      inv_dense <- solve(as.matrix(Sigma_inv_v_k))
      diag(inv_dense)
    }, error = function(e) {
      # Fallback: per-row solve (slower but more robust)
      diag(N) |> apply(2L, function(e_i) {
        as.numeric(solve(chol_factor, e_i))[which(e_i == 1)]
      })
    })
    state$diag_Sigma_v[, k] <- Sigma_v_k_diag
  }
  state
}

# ---- Update (A_jk, gamma_jk) -------------------------------------------------
.update_A <- function(state, Y_ord, tau_inv) {
  N <- nrow(Y_ord); G <- ncol(Y_ord); K <- ncol(state$m_A)
  EA  <- state$omega * state$m_A
  mu_row <- matrix(state$m_mu, nrow = N, ncol = G, byrow = TRUE)

  for (k in seq_len(K)) {
    EA_minus_k <- EA;  EA_minus_k[, k] <- 0
    # Y_partial: N x G
    Y_partial <- Y_ord - mu_row - state$m_v %*% t(EA_minus_k)

    EV_k  <- state$m_v[, k]                                   # N-vector
    EVV_k <- state$diag_Sigma_v[, k] + EV_k^2                 # N-vector

    # f_jk = tau_inv_j * sum_i V_ik * Y_partial_ij
    # In matrix form: tau_inv * t(Y_partial) %*% EV_k
    f_k <- as.numeric(crossprod(Y_partial, EV_k)) * tau_inv   # length G

    # g_jk = tau_inv_j * sum_i E[V_ik^2]
    g_k <- sum(EVV_k) * tau_inv                                # length G

    # Slab posterior
    s_A2_k <- 1 / (g_k + 1 / state$nu[k])
    m_A_k  <- s_A2_k * f_k

    # Log Bayes factor (numerically safe)
    log_BF <- 0.5 * (log(pmax(s_A2_k, 1e-300)) - log(state$nu[k])) +
              0.5 * m_A_k^2 / pmax(s_A2_k, 1e-300)
    log_BF <- pmin(pmax(log_BF, -50), 50)
    logit_pi <- log(state$pi[k] / max(1 - state$pi[k], 1e-10))
    omega_k <- plogis(logit_pi + log_BF)

    state$m_A[, k]   <- m_A_k
    state$s_A2[, k]  <- s_A2_k
    state$omega[, k] <- pmin(pmax(omega_k, 1e-4), 1 - 1e-4)
  }
  state
}

# ---- Update tau_j^2 ---------------------------------------------------------
.update_tau2 <- function(state, Y_ord) {
  N <- nrow(Y_ord); G <- ncol(Y_ord); K <- ncol(state$m_A)
  EA  <- state$omega * state$m_A                              # G x K
  EAA <- state$omega * (state$m_A^2 + state$s_A2)             # G x K
  mu_row <- matrix(state$m_mu, nrow = N, ncol = G, byrow = TRUE)
  fitted_full <- mu_row + state$m_v %*% t(EA)                 # N x G

  resid_sq <- colSums((Y_ord - fitted_full)^2)                # length G

  # Variance correction term:
  # sum_i sum_k [E[A_jk^2]E[V_ik^2] - EA_jk^2 EV_ik^2]
  EVV <- state$diag_Sigma_v + state$m_v^2                     # N x K
  EAA_minus_EA2 <- EAA - EA^2                                 # G x K (>= 0)
  EVV_minus_EV2 <- state$diag_Sigma_v                         # N x K (>= 0)
  # Term 1: sum_i sum_k E[A_jk^2] E[V_ik^2] - sum_i sum_k EA_jk^2 EV_ik^2
  # = sum_k [(sum_i E[V_ik^2]) * E[A_jk^2] - (sum_i EV_ik^2) * EA_jk^2]
  sum_EVV <- colSums(EVV)                                     # K
  sum_EV2 <- colSums(state$m_v^2)                             # K
  var_corr_AV <- as.numeric(EAA %*% sum_EVV - (EA^2) %*% sum_EV2)  # length G

  var_corr_mu <- N * state$s_mu2                              # length G

  S_j <- resid_sq + var_corr_mu + var_corr_AV
  S_j <- pmax(S_j, 1e-8)

  state$a_tau_hat <- rep(state$a_tau + N / 2, G)
  state$b_tau_hat <- state$b_tau + 0.5 * S_j
  state
}

# ---- Update sigma_k^2 -------------------------------------------------------
.update_sigma2 <- function(state) {
  N <- nrow(state$m_v); K <- ncol(state$m_v)
  for (k in seq_len(K)) {
    Qk <- state$Qtilde[[k]]
    mv <- state$m_v[, k]
    # T_k = m_v^T Q m_v + tr(Q Sigma_v)
    quad <- as.numeric(crossprod(mv, Qk %*% mv))
    # tr(Q Sigma_v) approximated using diagonal Sigma_v (v1 simplification)
    # Exact would use sparse selected inversion at Q's pattern
    tr_term <- sum(diag(Qk) * state$diag_Sigma_v[, k])

    T_k <- quad + tr_term
    state$a_sigma_hat[k] <- state$a_sigma + N / 2
    state$b_sigma_hat[k] <- state$b_sigma + 0.5 * max(T_k, 1e-8)
    state$sigma2[k]      <- state$b_sigma_hat[k] / max(state$a_sigma_hat[k] - 1, 1)
  }
  state
}

# ---- EB updates for pi_k, nu_k ----------------------------------------------
.update_pi_nu <- function(state) {
  K <- ncol(state$omega)
  for (k in seq_len(K)) {
    state$pi[k] <- max(mean(state$omega[, k]), 1e-3)
    state$pi[k] <- min(state$pi[k], 1 - 1e-3)
    denom <- sum(state$omega[, k])
    if (denom < 1e-6) {
      state$nu[k] <- max(state$nu[k], 1e-3)
    } else {
      numer <- sum(state$omega[, k] * (state$m_A[, k]^2 + state$s_A2[, k]))
      state$nu[k] <- max(numer / denom, 1e-6)
    }
  }
  state
}

# ---- L-BFGS on alpha to update phi ------------------------------------------
.update_alpha <- function(state, locs_ord, NNarray, smoothness, max_steps = 10L) {
  K <- length(state$alpha)
  current_alpha <- state$alpha

  obj <- function(alpha) {
    # phi from cumulative-exp
    phi_k <- cumsum(exp(alpha))
    # Sum across factors of L_k(phi_k)
    val <- 0
    for (k in seq_len(K)) {
      qk <- .build_vecchia_precision(locs_ord, NNarray,
                                     phi = phi_k[k], smoothness = smoothness)
      mv <- state$m_v[, k]
      quad <- as.numeric(crossprod(mv, qk$Q_tilde %*% mv))
      tr_term <- sum(diag(qk$Q_tilde) * state$diag_Sigma_v[, k])
      T_k <- quad + tr_term
      sigma_inv_k <- state$a_sigma_hat[k] / state$b_sigma_hat[k]
      # ELBO contribution from factor k (up to constants):
      #   (1/2) log|Q_tilde| - (1/2) sigma_inv_k * T_k
      val <- val + 0.5 * qk$log_det_Q - 0.5 * sigma_inv_k * T_k
    }
    -val   # ucminf minimises
  }

  fit <- tryCatch(
    ucminf::ucminf(par = current_alpha, fn = obj,
                   control = list(maxeval = max_steps, trace = 0)),
    error = function(e) list(par = current_alpha, convergence = -99))

  state$alpha <- fit$par
  state$phi <- cumsum(exp(state$alpha))
  state
}

# ---- Sign-fix and diagnostics -----------------------------------------------
.sign_fix <- function(state) {
  K <- ncol(state$m_A)
  for (k in seq_len(K)) {
    # Use the gene with largest omega as anchor
    anchor <- which.max(state$omega[, k])
    if (length(anchor) && state$m_A[anchor, k] < 0) {
      state$m_A[, k]  <- -state$m_A[, k]
      state$m_v[, k]  <- -state$m_v[, k]
      # squared quantities (s_A2, EVV) unchanged
    }
  }
  state
}


# ============================================================================
# ELBO computation
# ============================================================================
.compute_elbo <- function(state, Y_ord, tau_inv) {
  N <- nrow(Y_ord); G <- ncol(Y_ord); K <- ncol(state$m_A)
  EA  <- state$omega * state$m_A
  EAA <- state$omega * (state$m_A^2 + state$s_A2)
  EVV <- state$diag_Sigma_v + state$m_v^2

  mu_row <- matrix(state$m_mu, nrow = N, ncol = G, byrow = TRUE)
  fitted_full <- mu_row + state$m_v %*% t(EA)
  resid_sq <- colSums((Y_ord - fitted_full)^2)

  # Variance corrections (same expressions as in tau2 update)
  sum_EVV <- colSums(EVV)
  sum_EV2 <- colSums(state$m_v^2)
  var_corr_AV <- as.numeric(EAA %*% sum_EVV - (EA^2) %*% sum_EV2)
  var_corr_mu <- N * state$s_mu2

  S_j <- resid_sq + var_corr_mu + var_corr_AV

  # E[log p(Y|.)]: -1/2 sum_{i,j} [ tau_inv_j * E[(Y-..)^2] + log tau^2 ]
  log_tau_terms <- log(state$b_tau_hat) - digamma(state$a_tau_hat)
  ll_data <- -0.5 * sum(tau_inv * S_j) - 0.5 * N * sum(log_tau_terms)

  # E[log p(v_k)]: -0.5 sigma_inv_k T_k - 0.5 sum_i log d_ki - N/2 E[log sigma_k^2]
  ll_v <- 0
  for (k in seq_len(K)) {
    sigma_inv_k <- state$a_sigma_hat[k] / state$b_sigma_hat[k]
    qk <- state$Qtilde[[k]]
    quad <- as.numeric(crossprod(state$m_v[, k], qk %*% state$m_v[, k]))
    tr_term <- sum(diag(qk) * state$diag_Sigma_v[, k])
    T_k <- quad + tr_term
    # log|Q| from Vecchia
    log_det_Q <- if (!is.null(state$vecchia_aux[[k]])) {
      state$vecchia_aux[[k]]$log_det_Q
    } else {
      determinant(as.matrix(qk), logarithm = TRUE)$modulus
    }
    log_sigma_term <- log(state$b_sigma_hat[k]) - digamma(state$a_sigma_hat[k])
    ll_v <- ll_v + 0.5 * log_det_Q - 0.5 * sigma_inv_k * T_k -
            0.5 * N * log_sigma_term
  }

  # E[log p(A, gamma)]: spike-and-slab
  ll_A <- 0
  for (k in seq_len(K)) {
    ll_A <- ll_A +
      sum(state$omega[, k] * (-0.5 * log(2 * pi * state$nu[k]) -
                              (state$m_A[, k]^2 + state$s_A2[, k]) /
                              (2 * state$nu[k]))) +
      sum(state$omega[, k] * log(state$pi[k] / max(1 - state$pi[k], 1e-12)) +
          log(max(1 - state$pi[k], 1e-12)))
  }

  # Entropies (negative because ELBO = E[log p] - E[log q])
  # H(q(mu)): 0.5 sum log(2 pi e s_mu^2)
  H_mu <- 0.5 * sum(log(2 * pi * exp(1) * state$s_mu2))
  # H(q(A, gamma)): Bernoulli x Gaussian
  H_A <- 0
  om <- state$omega
  H_A <- H_A - sum(om * log(om) + (1 - om) * log(1 - om))
  H_A <- H_A + sum(om * 0.5 * log(2 * pi * exp(1) * state$s_A2))
  # H(q(v_k)): -log|Q_v| where Q_v = Sigma_inv_v_k.  Approximate via
  # log|Sigma_v_k| ~ sum(log(diag_Sigma_v_k)) (diagonal approximation).
  H_v <- 0.5 * sum(log(2 * pi * exp(1) * pmax(state$diag_Sigma_v, 1e-12)))
  # H(q(tau^2)): inverse gamma entropy
  H_tau <- sum(state$a_tau_hat + log(state$b_tau_hat) +
               lgamma(state$a_tau_hat) -
               (1 + state$a_tau_hat) * digamma(state$a_tau_hat))
  H_sigma <- sum(state$a_sigma_hat + log(state$b_sigma_hat) +
                 lgamma(state$a_sigma_hat) -
                 (1 + state$a_sigma_hat) * digamma(state$a_sigma_hat))

  elbo <- ll_data + ll_v + ll_A + H_mu + H_A + H_v + H_tau + H_sigma
  elbo
}


# ============================================================================
# Main CAVI engine fit
# ============================================================================
fit_engine_cavi <- function(Y, locs, K,
                            m = 20L,
                            ordering = "maxmin",
                            smoothness = 1.5,
                            max_iter = 100L,
                            tol = 1e-5,
                            alpha_update_every = 5L,
                            D_orthogonalize = TRUE,
                            verbose = TRUE) {
  Y <- as.matrix(Y); locs <- as.matrix(locs)
  stopifnot(nrow(Y) == nrow(locs))
  N <- nrow(Y); G <- ncol(Y)

  if (verbose) {
    cat(sprintf("fit_engine_cavi: N=%d, G=%d, K=%d, m=%d\n", N, G, K, m))
  }

  # Initialise
  init <- .cavi_initialise(Y, locs, K, m, ordering, smoothness, verbose)
  state <- init$state

  # Reorder Y to Vecchia ordering once
  ord <- init$engine_pca_init$vecchia$ord
  Y_ord <- Y[ord, , drop = FALSE]
  locs_ord <- init$engine_pca_init$vecchia$locs_ord
  NNarray <- init$engine_pca_init$vecchia$NNarray

  # CAVI loop
  elbo_history <- numeric(max_iter)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    tau_inv <- state$a_tau_hat / state$b_tau_hat
    state <- .update_mu(state, Y_ord, tau_inv)
    state <- .update_v(state, Y_ord, tau_inv)
    state <- .update_A(state, Y_ord, tau_inv)
    state <- .update_tau2(state, Y_ord)
    state <- .update_sigma2(state)
    state <- .update_pi_nu(state)

    if (it %% alpha_update_every == 0L) {
      state <- .update_alpha(state, locs_ord, NNarray, smoothness)
      state <- .cavi_rebuild_Qtilde(state, locs_ord, NNarray, smoothness)
    }

    state <- .sign_fix(state)

    elbo_history[it] <- .compute_elbo(state, Y_ord, tau_inv)
    if (verbose && it %% 5 == 0L) {
      cat(sprintf("  iter %3d  ELBO = %.4f\n", it, elbo_history[it]))
    }

    if (it > 1L) {
      rel_change <- abs(elbo_history[it] - elbo_history[it - 1L]) /
                    max(abs(elbo_history[it - 1L]), 1e-8)
      if (rel_change < tol) {
        converged <- TRUE
        if (verbose) cat(sprintf("  converged at iter %d (rel change %.2e)\n",
                                 it, rel_change))
        elbo_history <- elbo_history[seq_len(it)]
        break
      }
    }
  }

  # Final outputs: package as standard engine triple
  A_hat <- state$omega * state$m_A
  tau2_hat <- state$b_tau_hat / pmax(state$a_tau_hat - 1, 1)
  theta_hat <- cbind(
    sigma2 = state$sigma2,
    phi    = state$phi,
    eta2   = rep(0.01, K)   # CAVI absorbs nugget into D, no separate eta
  )

  # Optional D-orthogonalisation post-processing
  ortho <- NULL
  if (D_orthogonalize) {
    od <- .d_orthogonalize(A_hat, tau2_hat)
    A_hat <- od$A_new
    state$m_v <- state$m_v %*% solve(t(od$rotation))
    ortho <- list(rotation = od$rotation, lambda_A = od$lambda)
  }

  # Inverse-reorder V_scores back to original indexing
  V_scores_orig <- matrix(0, N, K)
  V_scores_orig[ord, ] <- state$m_v

  out <- list(
    A_hat = A_hat,
    theta_hat = theta_hat,
    D_hat = tau2_hat,
    V_scores = V_scores_orig,
    Qtilde_list = state$Qtilde,
    vecchia = init$engine_pca_init$vecchia,
    X = NULL,
    mu_hat = state$m_mu,
    B_hat = NULL,
    ortho = ortho,
    fit_diagnostics = list(
      converged = converged, n_iter = length(elbo_history),
      elbo_history = elbo_history,
      smoothness = smoothness,
      omega = state$omega,
      pi = state$pi, nu = state$nu
    ),
    engine_type = "cavi_sparse_kernel_ordered",
    N = N, G = G, K = K
  )
  class(out) <- c("passage_engine", "list")
  out
}
