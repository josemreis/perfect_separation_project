require(tidyverse)

make_separation_plot <- function(n = 100, beta = 0.5, seed = 1234) {
  
  set.seed(seed)
  alpha <- 0 # centered data
  df <- tibble(x = rnorm(n, 0, 1)) %>%
    mutate(pi = plogis((alpha + beta * x)), ## P(Y = 1| X) = logit^-1 \pi(x) = \frac{exp(\alpha + x * \beta)}{1 + exp(\alpha + x * \beta)}
           y = rbinom(n, 1, prob = pi))
  
  p <- df %>%
    ggplot(aes(x = x, y = y)) + 
    geom_point() +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = .5, fill = "palegreen", alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = .5, fill = "#F8766D", alpha = 0.5)  +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = .5, fill = "#F8766D", alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = Inf, ymax = .5, fill = "palegreen", alpha = 0.5) +
    stat_smooth(method = "glm", 
                method.args = list(family = "binomial"),
                se = FALSE) +
    labs(subtitle = bquote(beta == .(beta)))
  
  return(p) 
}

params <- c(2, 4, 12, 22)
plots <- map(params, ~make_separation_plot(n = 50, beta = .x, seed = 1234))
cowplot::plot_grid(plotlist = plots)

