### set up and helpers ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
### packages
packs <- c("devtools", "tidyverse", "furrr", "logistf", "rstanarm", "tidybayes", "detectseparation")
for (pack in packs) {
  
  if (pack %in% installed.packages()){
    
    require(pack, character.only = TRUE)
    
  } else {
    
    t <- try(install.packages(pack))
    
    if (class(t) == "try-error") {
      
      # get from github "cran" mirror
      library(devtools)
      install_github(sprintf("cran/%s", pack))
      
    }
    require(pack, character.only = TRUE)
    
  }  
  
}

### load simulation data
load_data_list <- function() {
  ## define the data types
  col_types <- cols(
    population_id = col_character(),
    sample_id = col_character(),
    population_size = col_double(),
    sample_size = col_double(),
    true_beta = col_double(),
    y = col_double(),
    x = col_double(),
    pi = col_double(),
    separation = col_double(),
    quasi_separation = col_double(),
    random_seed = col_double(),
    prop_sep_sample = col_double()
  )
  ## read all datasets
  files <- list.files("data", full.names = TRUE)
  data <- map(files, ~try(read_csv(.x, col_types = col_types))) %>%
    set_names(files)
  return(data)
  
}


## filter out dfs below the relevant treshold
filter_tresh <- function(df, sep_treshold =  1000) {
  if (nrow(df) > 0) {
    samples <- df %>% distinct(sample_id)
  } else {
    samples <- data.frame()
  }
  msg <- "Sample treshold is set to %s, %s has %s. %s."
  if (nrow(samples) > sep_treshold) {
    sprintf(msg, sep_treshold, unique(df$population_id), nrow(samples), " Fitting the models.")
    out <- TRUE
  } else {
    sprintf(msg, sep_treshold, unique(df$population_id), nrow(samples), " Skipping it")
    out <- FALSE
  }
  return(out)
}

## make_dir
make_dir <- function(dirs) {
  for (dir_name in dirs) {
    if (!dir.exists(dir_name)){
      dir.create(dir_name)
    }
  }
}

## get samples
get_samples <- function(df, n = 1000, seed = 1234, all = FALSE) {
  if (all == FALSE) {
    ## get a random sample of n sample ids
    set.seed(seed)
    sampled_ids <- sample(unique(df$sample_id), n)
    ## fetch
    pulled <- df %>%
      filter(sample_id %in% sampled_ids)
  } else {
    pulled <- df
  }
  ## split the datasets and turn to list of dfs
  sample_list <- pulled %>%
    group_split(sample_id)
  ## return
  return(sample_list)
}

#### data generation -------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
sep_fun <- function(df, quasi_tresh = 0.90, kosmidis = TRUE) {
  
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
### simulate the population and sample
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
    
    random_seed <- sample(.Random.seed, 1) 
    
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
           population_id = paste(abs(random_seed), beta, sep = "_"), ## add beta to population id case there were overlapping seeds. Will take the absolute value of the seed for id purposes (original seed stored)
           sample_id = paste(sample_id, population_id, sep = "_"),
           separation = if_else(complete_separation == TRUE, "1", "0"),
           quasi_separation = if_else(complete_separation == FALSE & kosmidis_sep == TRUE, "1", "0")) %>%
    select(population_id, sample_id, population_size, sample_size, true_beta, y, x, pi, separation, quasi_separation, random_seed)
  return(samples)
}

### generate the data
make_simulation_datasets <- function(beta = beta, sample_size  = 25, population_size = 1000000, seed = NULL, new_dataset = FALSE) {
  
  ####################################################################################################################
  ### Generate the simulation data
  # sample_size: int; size of the samples to draw from the population
  # beta: numerical; true beta parameter value used to generate the population data, i.e. \pi for y ~ Binomial(1, \pi)
  # seed: int; random seed can be user defined or, if NULL, generated randomly within the function
  # population size: size of the population used to extract the sample_n/pop_n samples
  # new_dataset: delete the dataset, existing, and make a new one
  #####################################################################################################################
  
  ## filename
  filename <- sprintf("data/separation_samples_%s_%s.csv.gz", beta, sample_size)
  ## make a new dataset?
  if (new_dataset) {
    try(file.remove(filename))
  } 
  
  if (!file.exists(filename)) {
    
    ### generate the populations and draw pop/sample_size samples
    sample_list <- generate_separation_data(beta = beta, sample_size  = sample_size, population_size = 1000000)
    ## turn to one df
    all_samples <-  bind_rows(sample_list)
    ## separate all data from samples with separation
    sep_samples <- all_samples %>%
      filter(separation == "1" | quasi_separation == "1")
    ## check how manylength(unique(sep_samples$sample_id))/length(unique(all_samples$sample_id))
    n_sep <- length(unique(sep_samples$sample_id))
    n <- length(unique(all_samples$sample_id))
    print(n_sep)
    print(n_sep/n)
    ## add as columns
    sep_samples$prop_sep_sample <- n_sep/n
    # clear up some space
    rm(all_samples)
    ### export
    if (!dir.exists("data")){
      dir.create("data")
    }
    
    ## samples with separation issues
    write_csv(sep_samples, gzfile(filename))
    ## return the dataset
    out <- sep_samples 
    
  } else {
    
    out <- read_csv(filename)
    
  }
}
#### models ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
### calculate the coverage probability
coverage_probability <- function(ci_upper, ci_lower, true_beta){
  ## check if the true beta is contained in the relevant confidence interval
  coverage <- map2_dbl(ci_lower, ci_upper, 
                       ~if_else(.x <= true_beta & .y >= true_beta, 1, 0))
  ## calculate the proportion of intervals containint the true beta
  cov_prob <- mean(coverage)
  return(cov_prob)
}

### fit glm
fit_glm <- function(samples, link_function = c("logit", "probit", "identity")) {
  ################################################################################
  # fit simulated samples with separation issues and evaluate performance using MLE (logit or probit link functions)
  # confint method: "wald"
  # NA since they do not exist
  # inference performance metrics
  # * bias: E[\^{\beta}] - \beta
  # * bias %: 100 \times (\frac{\^{\beta}}{\beta} - 1)
  # * variance of estimates: var(\^{\beta})
  # * Mean-squared error: E[bias(\^{\beta})^2] var(\^{\beta}) + [Bias(\^{\beta})]^2
  # * Monte Carlo Standard Error Error: \sqrt{var(\^{\beta})}; see: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjA6-D81pbvAhVO4YUKHXJ_DJ8QFjAAegQIBxAD&url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fpmc%2Farticles%2FPMC3337209%2F&usg=AOvVaw1MmJaBCHVx0Ig2oiFGlCwg
  # * coverage_probability: \frac{1}{n_sim}\sum_{i=1}^nI(\^{\beta}_{i,\lower} \leq \beta \geq \beta}_{i,\upper})
  # Args:
  # 1) list containing data frames
  # 2) link: chr, link function
  # output:
  # tibble
  ################################################################################
  ## trials
  n <- length(samples)
  ## create place holders
  est_beta <- est_alpha <- est_se <- ci_upper <- ci_lower <- numeric(length = n)
  for (i in seq_along(samples)) {
    df <- samples[[i]]
    if (link_function != "identity") {
      fit <- glm(formula = y ~ x, family = binomial(link = link_function), data = df)
    } else {
      ## glm with identity link kept throwing non convergence related errors, run the LPM using a normal ols
      fit <- lm(formula = y ~ x, data = df)
    }
    est_beta[i] <- coef(fit)[2] # fetch the estimated beta
    est_alpha[i] <-coef(fit)[1] # fetch the estimated intercept
    est_se[i] <- sqrt(diag(vcov(fit))) # fetch the standard-errors (sqrt of the diagonal of the var-cov matrix)
    ci_upper[i] <- NA
    ci_lower[i] <- NA
  }
  ## assign
  # static values
  true_beta <- unique(df$true_beta)
  ev_beta <- mean(est_beta)
  median_beta <- median(est_beta)
  # compute the inference performance stats and assign to df
  out <- tibble(
    sample_size = unique(df$sample_size),
    n_trials = n,
    pop_id = unique(df$population_id),
    estimation_method = case_when(
      link_function == "identity" ~ "MLE lpm",
      link_function == "logit" ~ "MLE logit",
      link_function == "probit" ~ "MLE probit"
    ),
    point_estimate = "MLE",
    ci_method = "wald",
    ci_alpha = 0.95,
    coverage_probability = NA_integer_, # coverage_probability(ci_lower = ci_lower, ci_upper = ci_upper, true_beta = true_beta)
    true_beta = true_beta,
    ev_beta = ev_beta,
    median_beta = median_beta,
    bias = ev_beta - true_beta,
    percent_bias = 100*(ev_beta/true_beta - 1),
    empirical_variance_beta = var(est_beta),
    empirical_sd_beta = sd(est_beta),
    mse = var(est_beta) + (ev_beta - true_beta)^2,
    mcse = sd(est_beta)/sqrt(n)
  )
  
  return(out)
  
}
### fit pml
fit_f <- function(samples) {
  ################################################################################
  # fit simulated samples with separation issues and evaluate performance using penalized maximum likelihood
  # confint method: "profile"
  # inference performance metrics
  # * bias: E[\^{\beta}] - \beta
  # * bias %: 100 \times (\frac{\^{\beta}}{\beta} - 1)
  # * variance of estimates: var(\^{\beta})
  # * Mean-squared error: E[bias(\^{\beta})^2] var(\^{\beta}) + [Bias(\^{\beta})]^2
  # * Monte Carlo Standard Error Error: \sqrt{var(\^{\beta})}; see: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjA6-D81pbvAhVO4YUKHXJ_DJ8QFjAAegQIBxAD&url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fpmc%2Farticles%2FPMC3337209%2F&usg=AOvVaw1MmJaBCHVx0Ig2oiFGlCwg
  # * coverage_probability: \frac{1}{n_sim}\sum_{i=1}^nI(\^{\beta}_{i,\lower} \leq \beta \geq \beta}_{i,\upper})
  # Args:
  # 1) list containing data frames
  # output:
  # tibble
  ################################################################################
  ## trials
  n <- length(samples)
  ## create place holders
  est_beta <- est_alpha <- est_se <- ci_upper <- ci_lower <- numeric(length = n)
  for (i in seq_along(samples)) {
    df <- samples[[i]]
    fit <- logistf(y ~ x, pl = TRUE, firth = TRUE, data=df)
    est_beta[i] <- coef(fit)[2] # fetch the estimated beta
    est_alpha[i] <-coef(fit)[1] # fetch the estimated intercept
    est_se[i] <- sqrt(diag(vcov(fit))) # fetch the standard-errors (sqrt of the diagonal of the var-cov matrix)
    confints <- confint(fit, level = 0.95)
    ci_upper[i] <- confints[2,2]
    ci_lower[i] <- confints[2,1]
  }
  ## assign
  # static values
  true_beta <- unique(df$true_beta)
  ev_beta <- mean(est_beta)
  median_beta <- median(est_beta)
  # compute the inference performance stats and assign to df
  out <- tibble(
    sample_size = unique(df$sample_size),
    n_trials = n,
    pop_id = unique(df$population_id),
    estimation_method = "Penalized Maximum Likelihood Estimation",
    point_estimate = "PMLE",
    ci_method = "profile",
    ci_alpha = 0.95,
    coverage_probability = coverage_probability(ci_lower = ci_lower, ci_upper = ci_upper, true_beta = true_beta),
    true_beta = true_beta,
    ev_beta = ev_beta,
    median_beta = median_beta,
    bias = ev_beta - true_beta,
    percent_bias = 100*(ev_beta/true_beta - 1),
    empirical_variance_beta = var(est_beta),
    empirical_sd_beta = sd(est_beta),
    mse = var(est_beta) + (ev_beta - true_beta)^2,
    mcse = sd(est_beta)/sqrt(n)
  )
  
  return(out)
  
}

### fit bayes
fit_b <- function(samples, beta_prior = student_t(df = 1, location = 0, scale = 2.5), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = TRUE) {
  ################################################################################
  # fit simulated samples with separation issues and evaluate performance using bayesian approaches
  # since we will obtain a posterior distribution of coefficients, we will compute the performance statistic using the posterior mode (MAP)
  # confint method: "PHighest density interval" 
  # inference performance metrics
  # * bias: E[\^{\beta}] - \beta
  # * bias %: 100 \times (\frac{\^{\beta}}{\beta} - 1)
  # * variance of estimates: var(\^{\beta})
  # * Mean-squared error: var(\^{\beta}) + [Bias(\^{\beta})]^2
  # * Monte Carlo Standard Error Error: \sqrt{var(\^{\beta})}; see: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjA6-D81pbvAhVO4YUKHXJ_DJ8QFjAAegQIBxAD&url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fpmc%2Farticles%2FPMC3337209%2F&usg=AOvVaw1MmJaBCHVx0Ig2oiFGlCwg
  # * coverage_probability: \frac{1}{n_sim}\sum_{i=1}^nI(\^{\beta}_{i,\lower} \leq \beta \geq \beta}_{i,\upper})
  # Args:
  # 1) list containing data frames
  # 2) prior distribution for the alpha and beta coefficient
  # output:
  # tibble
  ################################################################################
  
  ## trials
  n <- length(samples)
  ## create place holders
  est_beta <- est_alpha <- est_se <- ci_upper <- ci_lower <- numeric(length = n)
  for (i in seq_along(samples)) {
    df <- samples[[i]]
    ## if applicable, transform x so it is centered and scaled to have a mean o 0 and sd of 0.5
    if (gelman_transform) {
      df <- df %>%
        mutate(x = (x - mean(x)) * (0.5/sd(x)))
    }
    ## fit
    cores <- future::availableCores() - 1
    fit <- stan_glm(y ~ x, data = df, 
                    family = binomial(link = "logit"), 
                    prior = beta_prior, prior_intercept = alpha_prior,  
                    cores = cores, seed = 12345, algorithm="sampling")
    
    ## get point estimates and confidence intervals 
    beta_df <- fit %>% spread_draws(x) %>% mode_hdi(.width = .95)
    ## assign
    est_beta[i] <- beta_df$x[1]
    est_se[i] <- se(fit)
    ci_upper[i] <- beta_df$.upper[1]
    ci_lower[i] <- beta_df$.lower[1]
  }
  
  ## assign
  # static values
  true_beta <- unique(df$true_beta)
  ev_beta <- mean(est_beta)
  median_beta <- median(est_beta)
  prior_info <- head(unlist(beta_prior), -1)
  # compute the inference performance stats and assign to df
  out <- tibble(
    sample_size = unique(df$sample_size),
    n_trials = n,
    pop_id = unique(df$population_id),
    estimation_method = sprintf("Bayesian logit model, beta ~ %s(%s)", prior_info, paste(prior_info[-1], collapse = ", "))[1],
    point_estimate = "MAP",
    ci_method = "highest density posterior interval",
    ci_alpha = 0.95,
    coverage_probability = coverage_probability(ci_lower = ci_lower, ci_upper = ci_upper, true_beta = true_beta),
    true_beta = true_beta,
    ev_beta = ev_beta,
    median_beta = median_beta,
    bias = ev_beta - true_beta,
    percent_bias = 100*((ev_beta/true_beta) - 1),
    empirical_variance_beta = var(est_beta),
    empirical_sd_beta = sd(est_beta),
    mse = var(est_beta) + (ev_beta - true_beta)^2,
    mcse = sd(est_beta)/sqrt(n)
  )
  
  return(out)
  
}

### run the simulation
run_simulation <- function(bottom_up = FALSE, sep_treshold =  1000, random_samples = FALSE, n_samples = NULL, seed = 1234) {
  ########################################################################################################################
  # Run the simulation experiment
  # It stores the results in "results/simulation/output" in the project directory using the population id as base filename
  # args:
  # * bottom_up; (logical) given a list of samples should we start from bottom up? Usefull for running two machines parallel
  # * sep_treshhold; (int) minimum number of samples with perfect separation for a population to be considered for the simulation
  # * random_samples: (logical) use all separated samples or a random sample?
  # * n_samples: (int) if yes to the previous, how many?
  # * seed: (int) random seed for sampling the separation samples
  ###############################################################################################################################
  
  ## load the data
  df_list <- load_data_list()
  ## progress
  pb = txtProgressBar(min = 0, max = length(df_list), initial = 0)
  ## define the order of the simulation experiments
  if (bottom_up) {
    indicators <- rev(seq_along(df_list))
  } else {
    indicators <- seq_along(df_list)
  }
  ## population level loop
  for (pop_i in indicators) {
    stepi <- pop_i
    setTxtProgressBar(pb,stepi)
    ## extract the df with the samples associated with a population simulation
    pop <- df_list[[pop_i]]
    ## if missing create the output repos
    make_dir(c("results", "results/simulation", "results/simulation/output"))
    ## if empty or below treshold, remove file and skip
    filter_test <- filter_tresh(pop, sep_treshold = sep_treshold)
    if (filter_test == FALSE) {
      next
    }
    ## get samples if applicable else just turn the samples df into list with samples as items
    complete <- if_else(random_samples == FALSE, TRUE, FALSE)
    samples <- get_samples(df = pop, n = n_samples, seed = seed, all = complete)
    ## run the simulations
    # make the filenames
    model_names <- c("glm_logit", "glm_identity", "glm_probit", "pml", "bayes_cauchy", "bayes_student_t", "bayes_normal_weakly_informative")
    filenames <- map_chr(model_names, ~paste0(sprintf("results/simulation/output/%s_results_",.x), str_remove(str_split(names(df_list[pop_i]), "data/")[[1]][2], ".gz"))) %>%
      set_names(model_names)
    # run
    for (model in model_names) {
      file <- filenames[model]
      if (!file.exists(file)) {
        cat(sprintf("\n> Fitting %s model on the the samples with separation from population %s\n\n", model, samples[[1]]$population_id[1]))
        if (model ==  "glm_logit") {
          out <- try(fit_glm(samples = samples, link_function = "logit"))
        } else if (model == "glm_identity") {
          out <- try(fit_glm(samples = samples, link_function = "identity"))
        } else if (model == "glm_probit") {
          out <- try(fit_glm(samples = samples, link_function = "probit"))
        } else if (model == "pml") {
          out <- try(fit_f(samples = samples))
        } else if (model == "bayes_cauchy") {
          out <- try(fit_b(samples = samples, beta_prior = cauchy(location = 0, scale = 2.5), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = TRUE))
        } else if (model == "bayes_student_t") {
          out <- try(fit_b(samples = samples, beta_prior = student_t(df = 7), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = FALSE))
        } else {
          out <- try(fit_b(samples = samples, beta_prior = normal(location = 0, 2.5), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = FALSE))
        }
        ## export
        if (class(out) != "try-error") {
          write_csv(out, file)
        }
        
      }
      
    }
    
  }
}

### make dataset
make_results_dataset <- function(dataset_path = "results/simulation/simulation_results.csv") {
  
  if (!file.exists(dataset_path)) {
   
    ## list the output files
    files <- list.files("results/simulation/output", full.names = TRUE)
    ## read, combine and clean up
    out <- map_df(files, read_csv) %>%
      mutate(estimation_method = str_replace(estimation_method, "NA, ", ""))
    ## export
    write_csv(out, dataset_path)
    
  } else {
    
    out <- read_csv(dataset_path)
    
  }
  
  return(out)
  
}
#### run -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### generate the simulation data
## prepare the input grid: varying population parameters and sample sizes
if (!dir.exists("data")) {
  set.seed(1234)
  betas <- sample(seq(from = 2, to = 6, by = 0.001), 100)
  sample_size <- c(25, 50, 75)
  simulation_grid <- expand_grid(betas, sample_size)
  ### generate the populations N = 1.000.000 and draw pop/sample_size samples given the population parameter. Store the ones with perfect separation
  sample_list <- map2(simulation_grid$betas, simulation_grid$sample_size, ~make_simulation_datasets(beta = .x, sample_size  = .y, population_size = 10000000, seed = NULL, new_dataset = FALSE))
}
### fit the models
## select the order of model fitting conditional on the machine
order <- if_else(path.expand('~') == "/home/jmr", TRUE, FALSE)
## run
run_simulation(bottom_up = order, sep_treshold =  1000, random_samples = TRUE, n_samples = 1000, seed = 1234)
data <- make_results_dataset()