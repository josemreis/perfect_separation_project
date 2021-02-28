### packages
packs <- c("tidyverse", "furrr", "detectseparation")
for (pack in packs) {
  
  if (pack %in% installed.packages()){
    
    require(pack, character.only = TRUE)
    
  } else {
    
    install.packages(pack)
    require(pack, character.only = TRUE)
    
  }
}

### function for generating logistic regression data
logistic_sim <- function(n = 100, beta1 = 0.5, seed = 1234){
  
  set.seed(seed)
  alpha <- 0 # centered data
  df <- tibble(x = rnorm(n, 0, 1)) %>%
    mutate(pi = plogis((alpha + beta * x)), ## P(Y = 1| X) = logit^-1 \pi(x) = \frac{exp(\alpha + x * \beta)}{1 + exp(\alpha + x * \beta)}
           y = rbinom(n, 1, prob = pi))
  
  return(df)
}
## count the number of samples with complete separation out of 100 samples of size 25
sep_fun <- function(df, quasi_tresh = 0.92, kosmidis = FALSE) {
  
  ## check if there is complete separation: Rule based
  out <- df %>%
    mutate(overlap_complete = if_else(x < 0 & y == 0 | x > 0 & y == 1,1,0),
           overlap_quasi = if_else(x <= 0 & y == 0 | x >= 0 & y == 1,1,0)) %>%
    mutate(n_complete = sum(overlap_complete),
           n_quasi = sum(overlap_quasi)) %>%
    mutate(complete_separation = n_complete == nrow(.),
           tresh_quasi_separation = n_quasi > nrow(.) * quasi_tresh)
  
  ## check if there is separation using Kosmidis and Konis algorithm
  if (kosmidis & out$tresh_quasi_separation[1] == TRUE) {
    
    fit <- glm(y ~ x, family = binomial(link = "logit"), method = "detect_separation", data = out)
    out$kosmidis_sep <- coef(fit)["x"] %in% c(Inf, -Inf)
    
  } else {
    
    out$kosmidis_sep <- FALSE
    
  }
  
  return(out)  
}

generate_separation_data <- function(population_size = 10000, sample_size = 50, beta = 4, seed = NULL){
  
  ### simulate population data
  if (length(seed) == 0) {
    
    random_seed <- sample(1:99999999, 1) 
    
  } else {
    
    random_seed <- seed
    
  }
  ### simulate population level data with a high regression coefficient
  population <- logistic_sim(n = population_size, beta = beta, seed = random_seed)
  ### draw several samples and assess whether there is identification: (i) rules based; (ii) if (i) assess using kosmidis et al algorithm
  set.seed(random_seed)
  reps <- population_size/sample_size ## repeated sampling: n repetitions
  ### paralellization
  no_cores <- availableCores() - 1
  plan(multicore, workers = no_cores)
  ### prepare the samples data, check if there is separation, and clean up
  samples <- replicate(reps, sample_n(population, size = sample_size), simplify = FALSE) %>% ## sampling
    furrr::future_map2_dfr(., 1:reps, ~ sep_fun(df = .x, kosmidis = TRUE) %>% ## assess separation
                             mutate(sample_id = .y)) %>%
    mutate(true_beta = beta, ## add the value of the true beta parameter, sample size, as well as the sample seed used for obtaining the population and sample data 
           population_size = population_size,
           sample_size = sample_size,
           random_seed = random_seed,
           population_id = paste(random_seed, beta, sep = "_"), ## add beta to population id case there were overlapping seeds
           sample_id = paste(sample_id, population_id, sep = "_"),
           separation = if_else(complete_separation == TRUE, "1", "0"),
           quasi_separation = if_else(complete_separation == FALSE & kosmidis_sep == TRUE, "1", "0")) %>%
    select(population_id, sample_id, population_size, sample_size, true_beta, y, x, pi, separation, quasi_separation, random_seed)
  return(samples)
}

### generate the data
## prepare the input grid: varying population parameters and sample sizes
set.seed(1234)
betas <- sample(seq(from = 0, to = 6, by = 0.01), 20)
sample_size <- c(25, 50, 75, 100)
simulation_grid <- expand_grid(betas, sample_size)
### generate the populations and draw pop/sample_size samples
sample_list <- map2(simulation_grid$betas, simulation_grid$sample_size, ~generate_separation_data(beta = .x, sample_size  = .y))
## turn to one df
all_samples <-  bind_rows(sample_list)
## separate all data from samples with separation
sep_samples <- all_samples %>%
  filter(separation == "1" | quasi_separation == "1")
### export
if (!dir.exists("data")){
  dir.create("data")
}
## all samples
write_csv(all_samples, gzfile("data/all_samples.csv.gz"))
## samples with separation issues
write_csv(sep_samples, gzfile("data/separation_samples.csv.gz"))

