#' Log-rank test
#'
#'
#'
#' @export
#'
m3_surv_logrank <- function(dta_surv, test_arms = rbind(c(1, 2), c(1, 3))) {
    f_test <- function(dta) {
        rst     <- survdiff(Surv(obs_time, obs_event) ~ arm, data = dta)
        z_score <- sqrt(rst$chisq)

        if (rst$obs[2] > rst$exp[2])
            z_score <- -z_score

        z_score
    }

    rst <- NULL
    for (i in seq_len(nrow(test_arms))) {
        cur_dta <- dta_surv %>%
            filter(arm %in% test_arms[i, ])

        cur_rst <- f_test(cur_dta)
        rst     <- c(rst, cur_rst)
    }

    rst
}


#' Asymptotic co-variance matrix
#'
#' Get covariance of [z_2(t) z_2(1) z_1(t)]
#'
#' @export
#'
m3_asym_cov <- function(ts = 0.5) {
    rst       <- diag(3)
    rst[1, 2] <- rst[2, 1] <- 1 / 2
    rst[1, 3] <- rst[3, 1] <- sqrt(ts)
    rst[2, 3] <- rst[3, 2] <- sqrt(ts) / 2
    rst
}
