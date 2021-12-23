#' count events
#'
#'
m3_count_events <- function(dta) {
    dta %>%
        group_by(arm) %>%
        summarize(n_event = sum(obs_event))
}

#' add test names
#'
#'
m3_add_testnames <- function(dta, vec_name) {
    colnames(dta) <- c(vec_name,
                       paste("R", 1:11, sep = ""))
    rst <- data.frame(dta) %>%
        gather(Type, Prob,
               R1, R2, R3, R4, R5,
               R6, R7, R8, R9, R10,
               R11)
    rst
}
