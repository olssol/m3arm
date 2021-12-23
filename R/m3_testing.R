#' Test result
#'
#' @export
#'
m3_rej_interim <- function(smps_zs, ob, b2) {
    f_r <- function(zs) {
        r    <- rep(NA, 11)
        r[1] <- zs[1] > ob[1]     ## inter rej at ori bound arm 2
        r[2] <- zs[1] > b2        ## inter rej at alt bound arm 2
        r[3] <- zs[3] > ob[2]     ## fin   rej at org bound arm 2
        r[4] <- zs[2] > ob[1]     ## inter rej at ori bound arm 1

        r[5] <- (!r[1]) & r[2]    ## inter arm 2 btn org B and alt B
        r[6] <- r[4]    & r[5]    ## arm 1 rej; arm 2 fail ori B, rej alt B
        r[7] <- (!r[1]) & r[3]    ## arm 2 fail int, rej final
        r[8] <- (!r[4]) & r[7]    ## arm 1/2 fail inter; arm 2 rej at final

        r[9]  <- r[1] | r[3]      ## gcd arm 2 rej
        r[10] <- r[1] | r[8]      ## gcd+termination arm 2 rej
        r[11] <- r[6] | r[10]     ## gcd+termination+alt B arm 2 rej

        r
    }

    b2  <- min(b2, ob[1])
    rst <- apply(rbind(smps_zs), 1, f_r)
    apply(rst, 1, mean)
}


#' Test for an individual trial
#'
#' @export
#'
m3_rej_trial_single <- function(zs, tp,
                                ob_type = "asOF",
                                b2_type = c("asP", "optimal"),
                                alpha = 0.0125, ...) {

    b2_type <- match.arg(b2_type)

    ob <- m3_ob_bound(tp = tp, alpha = alpha, type = ob_type)$ob
    b2 <- switch(b2_type,
                 asP     = m3_calculate_b2(tp, est_b2 = NULL,
                                           type = "asP", alpha = 0.125, ...),
                 optimal = m3_calculate_b2(tp, alpha = 0.125, ...))

    rej_12 <- m3_rej_interim(zs[1:3],        ob, b2)
    rej_13 <- m3_rej_interim(zs[c(2, 1, 4)], ob, b2)

    rst <- rbind(c(2, tp, rej_12),
                 c(3, tp, rej_13))
    rst
}

#' Test for all trials
#'
#' @export
#'
m3_rej_trial <- function(rst_trial, ...) {
    rst <- apply(rst_trial, 1, function(x) {
        zs  <- x[c("z21", "z31", "z22", "z32")]
        tp  <- x["tp"]
        rst <- m3_rej_trial_single(zs, tp, ...)
        data.frame(rst)
    })
    rst <- rbindlist(rst)
    rst <- m3_add_testnames(rst,
                            c("Arm", "InfoFrac"))

    rst_summary <- rst %>%
        group_by(Arm, Type) %>%
        summarize(InfoFrac = mean(InfoFrac),
                  RejProb  = mean(Prob))

    list(all     = rst,
         summary = rst_summary)
}
