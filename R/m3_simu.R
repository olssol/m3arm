#' Simulate patients
#'
#'
#' @export
#'
m3_simu_pts <- function(median_mth_pla   = 5,
                        arms_hr          = c(1, 1),
                        arms_randr       = c(1, 1),
                        drop_annual      = 0.05,
                        enroll_duration  = 342 / 21,
                        ntot             = 342,
                        ntot_event       = 157) {

    hazard_placebo  <- - log(0.5)  / median_mth_pla
    hazard_drop     <- - log(1 - drop_annual) / 12
    rand_ratio      <- c(1, arms_randr)
    rand_ratio      <- rand_ratio / sum(rand_ratio)
    rand_arm        <- apply(rmultinom(ntot, 1, rand_ratio),
                             2,
                             function(x) which(1 == x))

    hz_arm          <- hazard_placebo * c(1, arms_hr)
    hz_all          <- hz_arm[rand_arm]
    rand_enroll     <- runif(ntot, 0, enroll_duration)
    rand_censor     <- rexp(ntot, hazard_drop)
    rand_event      <- rexp(ntot, hz_all)

    rst_dta <- data.frame(arm      = rand_arm,
                          t_enroll = rand_enroll,
                          t_censor = rand_censor,
                          t_event  = rand_event) %>%
        rowwise() %>%
        mutate(obs_time    = min(t_censor, t_event),
               obs_event   = t_censor > t_event,
               obs_caltime = obs_time + t_enroll) %>%
        data.frame()

    ## stop when total number of events occurred
    cut_caltime <- m3_get_cut_caltime(rst_dta,
                                      n_event = ntot_event,
                                      arms    = 1:3)

    ## censor patients after number of events occurred
    rst <- m3_censor_cut(rst_dta, cut_caltime)

    ##
    rst
}

#' Find calendar time for given number of events
#'
#'
#' @export
#'
m3_get_cut_caltime <- function(dta_surv, n_event = Inf, arms = 1:3) {
    rst <- dta_surv %>%
        filter(arm %in% arms &
               1 == obs_event) %>%
        arrange(obs_caltime) %>%
        slice(min(n(), n_event))

    rst$obs_caltime
}

#' Censoring patients
#'
#'
#'
#' @export
#'
m3_censor_cut <- function(dta_surv, cut_caltime) {
    rst <- apply(dta_surv, 1, function(x) {
        obs_caltime <- x["obs_caltime"]
        if (cut_caltime < obs_caltime) {
            obs_time    <- cut_caltime - x["t_enroll"]
            obs_event   <- 0
            obs_caltime <- cut_caltime
        } else {
            obs_time    <- x["obs_time"]
            obs_event   <- x["obs_event"]
        }

        c(obs_time, obs_event, obs_caltime)
    })

    dta_surv[, c("obs_time", "obs_event", "obs_caltime")] <- t(rst)
    dta_surv %>%
        filter(obs_time > 0) %>%
        data.frame()
}

#' Simulate a trial and get z-sccores
#'
#' Get z-scores for one interim analysis
#'
#' @export
#'
m3_simu_trial_single <- function(..., ntot_event   = 157,
                                 int_nevent_arm_pt = 111,
                                 int_nevent_arm_p  = 70,
                                 tp = NULL) {

    f_e <- function(dta) {
        rst <- m3_count_events(dta)
        rst <- as.matrix(rst[ , 2])
    }

    simu_dta <- m3_simu_pts(..., ntot_event = ntot_event)

    ## interim analysis
    ## find interim time point
    if (!is.null(tp)) {
        cut_caltime <- m3_get_cut_caltime(simu_dta,
                                          ceiling(j * ntot_event),
                                          arms = 1:3)
    } else {
        cut_caltime_12 <- m3_get_cut_caltime(simu_dta,
                                             int_nevent_arm_pt,
                                             arms = c(1, 2))

        cut_caltime_13 <- m3_get_cut_caltime(simu_dta,
                                             int_nevent_arm_pt,
                                             arms = c(1, 3))

        cut_caltime_11 <- m3_get_cut_caltime(simu_dta,
                                             int_nevent_arm_p,
                                             arms = c(1))

        cut_caltime <- min(cut_caltime_11,
                           max(cut_caltime_12, cut_caltime_13))
    }

    interim_dta <- m3_censor_cut(simu_dta, cut_caltime)
    tp          <- sum(interim_dta$obs_event) / ntot_event


    cur_rst  <- NULL
    for (j in list(interim_dta, simu_dta)) {
        cur_zscore  <- m3_surv_logrank(j,
                                       test_arms = rbind(c(1, 2),
                                                         c(1, 3)))
        cur_rst     <- c(cur_rst, cur_zscore)
    }

    c(cur_rst, tp, f_e(interim_dta), f_e(simu_dta))
}

#' Simulate a trial and get z-sccores
#'
#' Get z-scores for one interim analysis
#'
#' @export
#'
m3_simu_trials <- function(nreps = 100000, ..., n_cores = 5) {
    rst <- parallel::mclapply(1:nreps,
                              function(k) {
                                  if (0 == k %% 500)
                                      print(k)
                                  m3_simu_trial_single(...)
                              }, mc.cores = n_cores)

    rst <- simplify2array(rst)
    rst <- t(rst)
    colnames(rst) <- c("z21", "z31", "z22", "z32",
                       "tp",
                       "d11", "d21", "d31",
                       "d12", "d22", "d32")
    rst
}
