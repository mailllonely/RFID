
log_likelihood_Y <- function(Y, PhiTheta, alpha, p, lambda) {
    .Call('_RFID_log_likelihood_Y', PACKAGE = 'RFID', Y, PhiTheta, alpha, p, lambda)
}

log_likelihood_Y_purity <- function(Y, PhiTheta, alpha, p, lambda) {
    .Call('_RFID_log_likelihood_Y_purity', PACKAGE = 'RFID', Y, PhiTheta, alpha, p, lambda)
}

likelihood_mat <- function(Y, PhiTheta, alpha, p, lambda) {
    .Call('_RFID_likelihood_mat', PACKAGE = 'RFID', Y, PhiTheta, alpha, p, lambda)
}

log_likelihood_alpha <- function(alpha, g0, h0) {
    .Call('_RFID_log_likelihood_alpha', PACKAGE = 'RFID', alpha, g0, h0)
}

log_likelihood_Phi <- function(Phi, eta) {
    .Call('_RFID_log_likelihood_Phi', PACKAGE = 'RFID', Phi, eta)
}

log_likelihood_Theta <- function(Theta, r, c) {
    .Call('_RFID_log_likelihood_Theta', PACKAGE = 'RFID', Theta, r, c)
}

rmultinomial <- function(N, probs) {
    .Call('_RFID_rmultinomial', PACKAGE = 'RFID', N, probs)
}

update_ell <- function(Phi, Theta, ELL) {
    .Call('_RFID_update_ell', PACKAGE = 'RFID', Phi, Theta, ELL)
}

rand_normal <- function(mean, stddev) {
    .Call('_RFID_rand_normal', PACKAGE = 'RFID', mean, stddev)
}

CrtRng <- function(n, r) {
    .Call('_RFID_CrtRng', PACKAGE = 'RFID', n, r)
}

update_ELL <- function(PhiTheta, x) {
    .Call('_RFID_update_ELL', PACKAGE = 'RFID', PhiTheta, x)
}

update_z <- function(alpha, PhiTheta, Y) {
    .Call('_RFID_update_z', PACKAGE = 'RFID', alpha, PhiTheta, Y)
}

update_m <- function(alpha, z, lambda) {
    .Call('_RFID_update_m', PACKAGE = 'RFID', alpha, z, lambda)
}

update_alpha <- function(g0, h0, m, p) {
    .Call('_RFID_update_alpha', PACKAGE = 'RFID', g0, h0, m, p)
}

update_Phi <- function(eta, ell_ik, Theta, p) {
    .Call('_RFID_update_Phi', PACKAGE = 'RFID', eta, ell_ik, Theta, p)
}

update_Theta <- function(r, c, p, ell_jk, Phi) {
    .Call('_RFID_update_Theta', PACKAGE = 'RFID', r, c, p, ell_jk, Phi)
}

update_p <- function(a0, b0, n, Theta, alpha) {
    .Call('_RFID_update_p', PACKAGE = 'RFID', a0, b0, n, Theta, alpha)
}

update_c <- function(e0, f0, r, Theta) {
    .Call('_RFID_update_c', PACKAGE = 'RFID', e0, f0, r, Theta)
}

update_gamma0 <- function(a0, b0, p_tilde, gamma0, ell_k) {
    .Call('_RFID_update_gamma0', PACKAGE = 'RFID', a0, b0, p_tilde, gamma0, ell_k)
}

update_c0 <- function(e0, f0, gamma0, r) {
    .Call('_RFID_update_c0', PACKAGE = 'RFID', e0, f0, gamma0, r)
}

update_ell_tilde <- function(ell_jk, r) {
    .Call('_RFID_update_ell_tilde', PACKAGE = 'RFID', ell_jk, r)
}

update_r <- function(c0, gamma0, p_tilde_j, ell_tilde_k) {
    .Call('_RFID_update_r', PACKAGE = 'RFID', c0, gamma0, p_tilde_j, ell_tilde_k)
}

