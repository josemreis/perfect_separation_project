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
logistic_sim <- function(n = 100, beta = 0.5, seed = 1234){
  
  ####################################################################################################################
  ### simulate population data for a logistic regression y|x ~ logit^{-1} \pi(x)
  # model details:
      # y|x ~ Binomial(1, \pi) \\
      # \pi = logit^{-1} \pi = \frac{exp(\alpha + \beta X)}{1 + exp(\alpha + \beta X)} \\
      # X ~ N(0,1) \\
      # \alpha = 0
  # n: int; size of the population data to simulate. 
  # beta: numerical; true beta parameter value used to generate the population data, i.e. \pi for y ~ Binomial(1, \pi)
  # seed: int; random seed user defined
  #####################################################################################################################
  
  set.seed(seed)
  alpha <- 0 # centered data
  df <- tibble(x = rnorm(n, 0, 1)) %>%
    mutate(pi = plogis((alpha + beta * x)), 
           y = rbinom(n, 1, prob = pi))
  
  return(df)
}

### function for detecting separation
sep_fun <- function(df, quasi_tresh = 0.92, kosmidis = FALSE) {
  
  ##################################################################################################################################################################
  # assess whether the logistic regression data contains separation issues
  # df: data frame object
  # quasi_tresh: numerical; treshhold to consider quasi-complete separation; proportion of overlap. This will be double checked using Kosmidis and Konis algorithm
  # kosmidis: logical. Use Kosmidis and Konis algorithm?
  ##################################################################################################################################################################
  
  ## check if there is complete separation: Rule based
  out <- df %>%
    mutate(overlap_complete = if_else(x < 0 & y == 0 | x > 0 & y == 1,1,0),
           overlap_quasi = if_else(x <= 0 & y == 0,1,0)) %>%
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

generate_separation_data <- function(population_size = 10000, sample_size = 50, beta = 4, seed = NULL) {
  
  ####################################################################################################################
  ### simulate population data and take #pop/sample_n samples of data for y|x ~ logit^{-1} \pi(x)
  # df: int; size of the population data to simulate. Uses logistic_sic, see above.
  # sample_size: int; size of the samples to draw from the population
  # beta: numerical; true beta parameter value used to generate the population data, i.e. \pi for y ~ Binomial(1, \pi)
  # seed: int; random seed can be user defined or, if NULL, generated randomly within the function
  #####################################################################################################################
  
  ### simulate population data
  if (length(seed) == 0) {
    
    random_seed <- sample(1:99999999, 1) 
    
  } else {
    
    random_seed <- seed
    
  }
  ### simulate population level data 
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

