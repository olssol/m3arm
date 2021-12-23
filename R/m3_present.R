#' Present grid search results
#'
#'
#'
#' @export
#'
m3_plot_b2_grid <- function(rst_b2_grid, types = c("R11"), alpha = 0.0125) {
    ggplot(data = rst_b2_grid %>% filter(Type %in% types),
           aes(x = Bound, y = Prob)) +
        geom_line(aes(group = InfoFrac, col = InfoFrac)) +
        geom_hline(yintercept = alpha) +
        theme_bw() +
        facet_wrap(~ Type)
}
