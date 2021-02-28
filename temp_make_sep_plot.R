### packages
packs <- c("tidyverse", "furrr")
for (pack in packs) {
  
  if (pack %in% installed.packages()){
    
    require(pack, character.only = TRUE)
    
  } else {
    
    install.packages(pack)
    require(pack, character.only = TRUE)
    
  }
}
## turn off scientific notation
options(scipen=999)

### function for generating logistic regression data
logistic_sim <- function(n = 100, beta1 = 0.5, seed = 1234){
  
  set.seed(seed)
  alpha <- 0 # centered data
  df <- tibble(x = rnorm(n, 0, 1)) %>%
    mutate(pi = plogis((alpha + beta * x)), ## P(Y = 1| X) = logit^-1 \pi(x) = \frac{exp(\alpha + x * \beta)}{1 + exp(\alpha + x * \beta)}
           y = rbinom(n, 1, prob = pi))
  
  return(df)
}
### function for making the "separation plots"
make_separation_plot <- function(n = 100, beta = 0.5, seed = 1234) {
  
  df <- logistic_sim(n = n, beta = beta, seed = seed)
  
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
    labs(subtitle = bquote(beta == .(beta) ~ "; " ~ n == .(n)),
         y = bquote(pi[(x)] == beta~X))
  
  return(p) 
}

# small sample and separation
params <- c(0.5, 1, 2, 4)
plots <- map(params, ~make_separation_plot(n = 25, beta = .x, seed = 123))
p1 <- cowplot::plot_grid(plotlist = plots)
p1
