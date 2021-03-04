#### set up and helpers
### ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
### packages
packs <- c("devtools", "tidyverse", "furrr", "logistf", "rstanarm", "tidybayes")
for (pack in packs) {
  
  if (pack %in% installed.packages()){
    
    require(pack, character.only = TRUE)
    
  } else {
    
    t <- try(install.packages(pack))
    
    if (class(t) == "try-error") {
      
      # get from a github mirror
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

#### models
### ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
fit_glm <- function(samples) {
  ################################################################################
  # fit simulated samples with separation issues and evaluate performance using MLE
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
  # output:
  # tibble
  ################################################################################
  ## trials
  n <- length(samples)
  ## create place holders
  est_beta <- est_alpha <- est_se <- ci_upper <- ci_lower <- numeric(length = n)
  for (i in seq_along(samples)) {
    df <- samples[[i]]
    fit <- glm(formula = y ~ x, family = binomial("logit"), data = df)
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
    estimation_method = "Maximum Likelihood Estimation",
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
    mse = mean((est_beta - true_beta)^2),
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
    mse = mean((est_beta - true_beta)^2),
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
  # * Mean-squared error: E[bias(\^{\beta})^2] var(\^{\beta}) + [Bias(\^{\beta})]^2
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
    mse = mean((est_beta - true_beta)^2),
    mcse = sd(est_beta)/sqrt(n)
  )
  
  return(out)
  
}

### run the simulation
run_simulation <- function(sep_treshold =  1000, random_samples = FALSE, n_samples = NULL, seed = 1234) {
  ## load the data
  df_list <- load_data_list()
  ## progress
  pb = txtProgressBar(min = 0, max = length(df_list), initial = 0) 
  ## population level loop
  for (pop_i in seq_along(df_list)) {
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
    model_names <- c("glm", "pml", "bayes_cauchy", "bayes_student_t", "bayes_normal_informative")
    filenames <- map_chr(model_names, ~paste0(sprintf("results/simulation/output/%s_results_",.x), str_remove(str_split(names(df_list[pop_i]), "data/")[[1]][2], ".gz"))) %>%
      set_names(model_names)
    # run
    for (model in model_names) {
      file <- filenames[model]
      if (!file.exists(file)) {
        cat(sprintf("\n> Fitting %s model on the the samples with separation from population %s\n\n", model, samples[[1]]$population_id[1]))
        if (model ==  "glm") {
          out <- fit_glm(samples = samples)
        } else if (model == "pml") {
          out <- fit_f(samples = samples)
        } else if (model == "bayes_cauchy") {
          out <- fit_b(samples = samples, beta_prior = cauchy(location = 0, scale = 2.5), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = TRUE)
        } else if (model == "bayes_student_t") {
          out <- fit_b(samples = samples, beta_prior = student_t(df = 7), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = FALSE)
        } else {
          out <- fit_b(samples = samples, beta_prior = normal(location = 1, 2), alpha_prior = student_t(df = 1, location = 0, scale = 10), gelman_transform = FALSE) 
        }
        ## export
        write_csv(out, file)
      }
      
    }

  }
}

#### run
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
run_simulation(sep_treshold =  1000, random_samples = TRUE, n_samples = 1000, seed = 1234)