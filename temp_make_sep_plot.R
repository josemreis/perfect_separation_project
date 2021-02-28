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
# large sample and separation
params <- c(0.5, 1, 2, 4, 6)
plots <- map(params, ~make_separation_plot(n = 500, beta = .x, seed = 123))
cowplot::plot_grid(plotlist = plots)

# separation as a consequence of the data
params <- c(2, 6, 14, 22)
plots <- map(params, ~make_separation_plot(n = 100, beta = .x, seed = 12))
cowplot::plot_grid(plotlist = plots)



### simulate a highly correlated population
POP_N <- 100000
SAMPLE_N <- 50
REPS <- POP_N/SAMPLE_N
BETA <- 5
population <- logistic_sim(n = POP_N, beta = BETA, seed = 123)
## plot
population %>%
  ggplot(aes(x = x, y = y)) + 
  geom_point() +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = .5, fill = "palegreen", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = .5, fill = "#F8766D", alpha = 0.5)  +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = .5, fill = "#F8766D", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = Inf, ymax = .5, fill = "palegreen", alpha = 0.5) +
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"),
              se = FALSE) +
  labs(subtitle = bquote(beta == .(BETA) ~ "; " ~ n == .(POP_N)),
       y = bquote(pi[(x)] == beta~X))
## count the number of samples with complete separation out of 100 samples of size 25
sep_fun <- function(df, quasi_tresh = 0.90, kosmidis = FALSE) {
  
  ## check if there is complete separation: Rule based
  out <- df %>%
    mutate(overlap_complete = if_else(x < 0 & y == 0 | x > 0 & y == 1,1,0),
           overlap_quasi = if_else(x <= 0 & y == 0 | x >= 0 & y == 1,1,0)) %>%
    mutate(n_complete = sum(overlap_complete),
              n_quasi = sum(overlap_quasi)) %>%
    mutate(complete_separation = n_complete == nrow(.),
           quasi_complete_separation = n_quasi > nrow(.) * quasi_tresh)

  ## check if there is separation using Kosmidis and Konis algorithm
  if (kosmidis & out$quasi_complete_separation[1] == TRUE) {
    
    fit <- glm(y ~ x, family = binomial(link = "logit"), method = "detect_separation", data = out)
    out$kosmidis_sep <- coef(fit)["x"] %in% c(Inf, -Inf)
    
  } else {
    
    out$kosmidis_sep <- NA
    
  }
  
  return(out)  
}

set.seed(1234)
plan("multicore")
### draw several samples and assess whether there is identification: (i) rules based; (ii) if (i) assess using kosmidis et al algorithm
samples <- replicate(REPS, sample_n(population, size = SAMPLE_N), simplify = FALSE) %>%
  furrr::future_map2_dfr(., 1:REPS, ~ sep_fun(df = .x, kosmidis = TRUE) %>%
           mutate(sample_id = .y)) 

### identify the samples containing separation
sep_samples<- samples %>%
  group_by(sample_id) %>%
  summarise(rule_complete = mean(complete_separation),
            rule_quasi = mean(quasi_complete_separation),
            kosmidis_sep = mean(kosmidis_sep)) %>%
  filter((rule_complete == 1 | rule_quasi == 1) & kosmidis_sep == 1)
### double check with kosmidis

