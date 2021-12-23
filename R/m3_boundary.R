#' Get One-Sided OB boundary
#'
#'
#' @export
m3_ob_bound <- function(tp = 0.65, alpha = 0.0125, type = "asOF") {
    ob_bound <- getDesignGroupSequential(sided            = 1,
                                         alpha            = alpha,
                                         informationRates = c(tp, 1),
                                         typeOfDesign     = type)

    list(ob         = ob_bound$criticalValues,
         alpha      = ob_bound$stageLevels,
         alpha_cumu = ob_bound$alphaSpent)
}


#' Boundary 2 for all grids
#'
#'
#'
#' @export
#'
m3_b2_grid_single <- function(samples, tp = 0.5, steps = 50, ..., n_cores = 5) {
    ob  <- m3_ob_bound(tp = tp, ...)$ob
    rst <- parallel::mclapply(seq(ob[1], ob[2], length.out = steps),
                              function(k) {
                                  rst <- m3_rej_interim(samples, ob, k)
                                  c(k, rst)
                              }, mc.cores = n_cores)

    rst <- simplify2array(rst)
    cbind(tp, t(rst))
}

#' Estimate b2 based on asymptotic distribution
#'
#'
#'
#' @export
#'
m3_b2_grid <- function(tps     = seq(0.5, 0.8, by = 0.05),
                       n_large = 1000000,
                       ...) {
    rst_all <- NULL
    for (tp in tps) {
        print(tp)
        cov_sig <- m3_asym_cov(tp)
        smps    <- mvrnorm(n_large, mu = c(0, 0, 0), Sigma = cov_sig)
        rst_b2  <- m3_b2_grid_single(smps, tp = tp, ...)
        rst_all <- rbind(rst_all, rst_b2)
    }

    bs_grid <- m3_add_testnames(rst_all, c("InfoFrac", "Bound")) %>%
        mutate(InfoFrac = factor(InfoFrac))

    b2_est  <- m3_est_b2_from_grid(bs_grid)

    ## return
    list(b2_grid = bs_grid,
         b2_est  = b2_est)
}



#' Find Boundary 2
#'
#'
#' @export
#'
m3_est_b2_from_grid <- function(rst_grid, alpha = 0.0125, ...) {
    rst_grid %>%
        filter(Type == "R11" &
               Prob <= alpha) %>%
        group_by(InfoFrac) %>%
        arrange(Bound) %>%
        slice(n = 1)

}


#' Get boundary of b2
#'
#'
#' @export
#'
m3_calculate_b2 <- function(tp, est_b2 = rbind(c(0.65, 2.58),
                                               c(0.70, 2.43),
                                               c(0.75, 2.37)),
                            type = "asP",
                            ...) {

    stopifnot(tp > 0 & tp < 1)

    if (is.null(est_b2)) {
        rst <- m3_ob_bound(tp, type = type, ...)
        return(rst$ob[1])
    }

    inx <- which(tp >= est_b2[, 1])
    if (0 == length(inx))
        return(Inf)

    inx <- max(inx)
    if (inx == nrow(est_b2))
        return(est_b2[nrow(est_b2), 2])

    x1 <- est_b2[inx, 1]
    y1 <- est_b2[inx, 2]
    x2 <- est_b2[inx + 1, 1]
    y2 <- est_b2[inx + 1, 2]

    rst <- y1 + (y2 - y1) * (tp - x1) / (x2 - x1)
    return(rst)

}
