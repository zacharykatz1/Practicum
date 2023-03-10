# To specify current index for dataset choice
args <- commandArgs(TRUE)
currind <- as.integer(args[1])
print(currind)

# May not be necessary
# getwd()
# .libPaths( c( "ifs/home/msph/biostat/zak2132/R/x86_64-pc-linux-gnu-library/3.5", 
# .libPaths()))
# .libPaths()

# Load relevant libraries
library(tidyverse)
library(rsample)
library(glmnet)
library(caret)
library(earth)
library(BART)
library(ranger)
library(mgcv)
library(SuperLearner)
library(qgcomp)
library(grf)
library(bkmr)

# Load data
simulated_data <- readRDS('/ifs/scratch/msph/biostat/zak2132/data/simulated_data.RDS')
quantiles <- readRDS('/ifs/scratch/msph/biostat/zak2132/data/quantiles.RDS')

# Parameters
cores = 400
M = 400 # Number of datasets per scenario
length = M / cores # Number of datasets per core
N = 1000 # Number of observations per dataset

# Relevant functions
# Functions for quantiles to use for estimation
metal_IQR_table = function(scenario, metal){
  
  # Average across datasets of each quantile for given metal in given scenario
  m.25 = quantiles %>% 
    filter(quantile == 0.25 & scenario == {{scenario}}) %>% 
    dplyr::select(c(metal + 2)) %>% 
    t() %>% 
    mean()
  
  m.75 = quantiles %>% 
    filter(quantile == 0.75 & scenario == {{scenario}}) %>% 
    dplyr::select(c(metal + 2)) %>% 
    t() %>% 
    mean()
  
  # For interquartile range change
  # Create tables with all vars at median except for given exposure to obtain potential outcomes (counterfactuals)
  q1.m = apply(filter(quantiles, quantile == 0.5 & scenario == {{scenario}}), 2, mean)
  q1.m[[c(metal + 2)]] = m.25
  q2.m = q1.m
  q2.m[[c(metal + 2)]] = m.75
  
  table = rbind(q1.m, q2.m) %>% as.data.frame() %>% dplyr::select(-scenario, -quantile)
  
  return(table)
  
}

metal_IQR_mix_table = function(scenario){
  
  # Average across datasets of each quantile for given metal in given scenario
  m.25 = quantiles %>% 
    filter(quantile == 0.25 & scenario == {{scenario}}) %>% 
    dplyr::select(M1:M10)
  
  m.75 = quantiles %>% 
    filter(quantile == 0.75 & scenario == {{scenario}}) %>% 
    dplyr::select(M1:M10)
  
  c.50 = quantiles %>% 
    filter(quantile == 0.50 & scenario == {{scenario}}) %>% 
    dplyr::select(C1:C5)
  
  q.1.means = apply(m.25, 2, mean)
  q.2.means = apply(m.75, 2, mean)
  c.means = apply(c.50, 2, mean)
  
  table = rbind(q.1.means, q.2.means) %>% as.data.frame()
  confound = rbind(c.means, c.means)
  table = cbind(table, confound)
  
  return(table)
  
}

# Function to bootstrap N times on a given model
boot_fnct = function(model_type, bootstrap_data, boot_count){
  
  set.seed(2132)
  
  boot = bootstrap_data %>% 
    dplyr::select(M1:Y) %>%
    bootstraps(times = boot_count) %>% 
    rowwise() %>% 
    mutate(data_sample = (list(analysis(splits)))) %>% 
    dplyr::select(id, data_sample) %>% 
    ungroup() %>%
    mutate(models = map(data_sample, model_type = model_type, model_fit_fnct)) %>% 
    dplyr::select(-data_sample)
  
  return(boot)
  
}

# Function to estimate causal effect using fitted model
predict_diff_fnct = function(model_type, fitted_model, q1, q2){
  
  if (model_type %in% c("linear","elastic net", "MARS", "GAM", "Random Forest", "GAM2", "Causal Forest")){
    
    predicted = predict(fitted_model, q2) - predict(fitted_model, q1)
    
  }
  
  else if (model_type %in% c("BART")){
    
    # For BART, SE medians from MC draws; is there a better way?
    predicted = quantile(predict(fitted_model, q2), 0.5) - quantile(predict(fitted_model, q1), 0.5)
    
  }
  
  else if (model_type %in% c("SuperLearner")){
    
    predicted = predict(fitted_model, q2, type = "response")$pred - predict(fitted_model, q1, type = "response")$pred
    
  }
  
  else if (model_type %in% c("Quantile G-Computation", "AMC G-Computation")){
    
    exposures = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
    
    predicted = predict(object = fitted_model, expnms = exposures, newdata = q2) - predict(object = fitted_model, expnms = exposures, newdata = q1)
    
  }
  
  else if (model_type %in% c("BKMR1")){
    
    predicted = 
      mean(
        SamplePred(fitted_model,
                   Znew = q2[,variable.names(fitted_model$Z)],
                   Xnew = as.matrix(q2[,variable.names(fitted_model$X)]))
        - 
          SamplePred(fitted_model,
                     Znew = q1[,variable.names(fitted_model$Z)],
                     Xnew = as.matrix(q1[,variable.names(fitted_model$X)]))
      )
    
  }
  
  else if (model_type %in% c("BKMR2")){
    
    predicted = 
      mean(
        SamplePred(fitted_model,
                   Znew = q2[,variable.names(fitted_model$Z)])
        - 
          SamplePred(fitted_model,
                     Znew = q1[,variable.names(fitted_model$Z)])
      )
    
  }
  
  return(predicted)
  
}

# Function to estimate causal effects using fitted model on bootstraps
predict_diff_boots = function(model_type, boots, q1, q2){
  
  bstar = NULL 
  n = dim(boots)[1];
  
  for (mod in 1:n) {
    
    if (model_type == "elastic net"){
      
      quant1 = predict(boots[[2]][[mod]], q1)
      quant1 = quant1[length(quant1)] # Pull for smallest lambda
      quant2 = predict(boots[[2]][[mod]], q2)
      quant2 = quant2[length(quant2)] # Pull for smallest lambda
      diff = tibble(diff = quant2 - quant1)[[1]]
      
    }
    
    else if (model_type == "BART"){
      
      # Use medians; is there a better method?
      quant1 = quantile(predict(boots[[2]][[mod]], q1), 0.5)
      quant2 = quantile(predict(boots[[2]][[mod]], q2), 0.5)
      diff = tibble(diff = quant2 - quant1)[[1]]
      
    }
    
    else {
      
      diff = tibble(diff = predict_diff_fnct(
        model_type,
        fitted_model = boots[[2]][[mod]],
        q1,
        q2)[[1]]
      )
      
    }
    
    bstar = rbind(bstar, diff)
  } # Next draw
  
  return(bstar)
  
}

# Function to filter to correct dataset
dataset_filter = function(scenario){
  
  # Filter for scenario and then to each individual dataset for looping
  first_index = (scenario * M) - (M - 1) 
  last_index = first_index + (M - 1)
  indices = seq(first_index, last_index, by = 1)
  filtered_datasets = simulated_data[indices]
  
  return(filtered_datasets)
  
}

# Modeling functions
# Function to fit models
model_fit_fnct = function(model_type, data){
  
  if (model_type %in% c("linear")){
    
    fit = lm(Y ~ ., 
             data = data)
    
  }
  
  else if (model_type %in% c("elastic net")){
    
    x_data = data %>% 
      dplyr::select(-Y) %>% 
      as.matrix()
    
    y_data = data %>% 
      dplyr::select(Y) %>% 
      as.matrix()
    
    # Force covariates into model (specify covariates not to drop; 0 for vars to keep)
    p_fac = rep(1, dim(x_data)[2])
    p_fac[11:15] = 0 # Keep confounders in model
    
    # Vector of values identifying what fold each observation is in for 10-fold cross-validation
    fold_id = sample(x = 1:10, size = length(y_data), replace = TRUE)
    
    # Fit model to training data with cross-validation
    # Alpha = 1 by default
    cv_fit = cv.glmnet(
      x = x_data,
      y = y_data,
      family = "gaussian",
      maxit = 5000,
      penalty.factor = p_fac,
      foldid = fold_id 
    )
    
    # Need to hard code cross-validation over range of alpha values
    for (i in seq(0, 1, by = 0.1)){
      cv_temp = cv.glmnet(
        x = x_data,
        y = y_data,
        family = "gaussian",
        maxit = 5000,
        penalty.factor = p_fac,
        foldid = fold_id,
        alpha = i
      )
      
      if(min(cv_temp$cvm) < min(cv_fit$cvm)) 
      {cv_fit = cv_temp}
      
    }
    
    # Optimal beta coefficients
    selected_beta_coefs = array(t(as.matrix(coef(cv_fit,
                                                 s = cv_fit$lambda.min))))
    
    # Final trained model for prediction
    # Returns elastic net model with coefficients across grid of values for regularization parameter lambda
    fit = glmnet(
      x = x_data,
      y = y_data,
      family = "gaussian",
      init = selected_beta_coefs,
      iter = 0
    )
    
  }
  
  else if (model_type %in% c("MARS")){
    
    train_data = data %>% 
      as.matrix()
    
    # Specify confounders to force into the model
    covars = c("C1", "C2", "C3", "C4", "C5") 
    
    # Implement 5-fold CV using caret package
    myControl = trainControl(
      method = "cv",
      number = 5,
      summaryFunction = defaultSummary
    )
    
    # Set hyperparameter tuning grid, may want to consider different values
    hyper_grid = expand.grid(
      degree = c(1,2), # Limit to up to second-order polynomials
      nprune = seq(5, 50, length.out = 10) %>% floor()
    )
    
    # Fit model
    fit = train(
      Y ~ .,
      data = train_data,
      method = "earth",
      linpreds = covars, # Enter covariates linearly
      thresh = 0, # Include a predictor even if it has very little predictive power
      penalty = -1, # Do not discard any terms in backwards pass
      trControl = myControl,
      tuneGrid = hyper_grid
    )
    
  }
  
  else if (model_type %in% c("GAM")){
    
    x_train = data %>% 
      dplyr::select(M1:C5) %>% 
      as.matrix()
    
    y_train = data %>% 
      dplyr::select(Y) %>% 
      as.matrix() %>% 
      as.numeric()
    
    # Cross-validation for smoothing functions
    ctrl = trainControl(method = "cv", number = 10)
    
    # Grid for hyperparameters
    # select = TRUE/FALSE for variable selection using GCV
    grid = data.frame(method = "GCV.Cp", select = c(TRUE,FALSE))
    
    # Fit GAM
    fit = train(
      x = x_train, # By default, GAM adds s() term for predictors having > 10 values, linear for those that do not
      y = y_train,
      method = "gam", # MGCV implementation of GAM
      tuneGrid = grid,
      trControl = ctrl
    )
    
  }
  
  else if (model_type %in% c("BART")){
    
    x_train = data %>%
      dplyr::select(M1:C5) %>%
      as.matrix()
    
    y_train = data %>%
      dplyr::select(Y) %>%
      as.matrix()
    
    fit = mc.wbart( # Use mc.wbart for parallel computation (just regular wbart without)
      x.train = x_train,
      y.train = y_train,
      x.test = x_train,
      ntree = 50, # Number of trees in the sum
      nskip = 250, # Number of MCMC iterations to be treated as burn in
      ndpost = 1000, # Number of posterior draws returned
      keepevery = 250,
      seed = 2132,
      mc.cores = 4 # Parallel computing
    )
    
  }
  
  else if (model_type %in% c("Random Forest")){
    
    train_data = data %>% 
      as.matrix()
    
    ctrl = trainControl(method = "cv") 
    
    rf_grid = expand.grid(
      mtry = 1:10, # of predictors at each split
      splitrule = "variance", # RSS
      min.node.size = 1:6 # Controls size of tree
    )
    
    fit = train(
      Y ~ ., # Assumes no interaction terms?
      data = train_data,
      method = "ranger",
      tuneGrid = rf_grid,
      trControl = ctrl
    )
    
  }
  
  else if (model_type %in% c("GAM2")){
    
    # Fit GAM using MGCV
    fit = mgcv::gam(
      # k = 4 is max polynomial for smoothing basis
      Y ~ ti(M1, k = 4) + ti(M2, k = 4) + ti(M3, k = 4) + ti(M4, k = 4) + ti(M5, k = 4) + ti(M6, k = 4) + ti(M7, k = 4) + ti(M8, k = 4) + ti(M9, k = 4) + ti(M10, k = 4) + 
        # Add two-way interactions for tensor smoothing
        ti(M1, M2, k = 4) + ti(M1, M3, k = 4) + ti(M1, M4, k = 4) + ti(M1, M5, k = 4) + ti(M1, M6, k = 4) + ti(M1, M7, k = 4) + ti(M1, M8, k = 4) + ti(M1, M9, k = 4) + ti(M1, M10, k = 4) + ti(M2, M3, k = 4) + ti(M2, M4, k = 4) + ti(M2, M5, k = 4) + ti(M2, M6, k = 4) + ti(M2, M7, k = 4) + ti(M2, M8, k = 4) + ti(M2, M9, k = 4) + ti(M2, M10, k = 4) + ti(M3, M4, k = 4) + ti(M3, M5, k = 4) + ti(M3, M6, k = 4) + ti(M3, M7, k = 4) + ti(M3, M8, k = 4) + ti(M3, M9, k = 4) + ti(M3, M10, k = 4) + ti(M4, M5, k = 4) + ti(M4, M6, k = 4) + ti(M4, M7, k = 4) + ti(M4, M8, k = 4) + ti(M4, M9, k = 4) + ti(M4, M10, k = 4) + ti(M5, M6, k = 4) + ti(M5, M7, k = 4) + ti(M5, M8, k = 4) + ti(M5, M9, k = 4) + ti(M5, M10, k = 4) + ti(M6, M7, k = 4) + ti(M6, M8, k = 4) + ti(M6, M9, k = 4) + ti(M6, M10, k = 4) + ti(M7, M8, k = 4) + ti(M7, M9, k = 4) + ti(M7, M10, k = 4) + ti(M8, M9, k = 4) + ti(M8, M10, k = 4) + ti(M9, M10, k = 4) + 
        # Add covariates
        C1 + C2 + C3 + C4 + C5,
      data = data,
      method = "REML"
    )
    
  }
  
  else if (model_type %in% c("SuperLearner")){
    
    x_train = data %>% 
      dplyr::select(M1:C5) %>% 
      as.data.frame()
    
    y_train = data %>% 
      dplyr::select(Y) %>% 
      as.matrix() %>% 
      as.numeric()
    
    # Select variety of candidate learners, some parametric, some non-parametric
    SL_library_chosen = c("SL.glmnet", "SL.randomForest", "SL.glm", "SL.ranger", "SL.xgboost", "SL.svm", "SL.earth", "SL.gbm", "SL.ranger")
    
    # Select loss function for meta learner and estimate risk
    # Here, use non-negative least squares loss
    loss_chosen = "method.NNLS"
    
    # Fit WITHOUT cross-validation (try with CV later to improve stability / robustness?)
    fit = SuperLearner(
      Y = y_train,
      X = x_train,
      SL.library = SL_library_chosen,
      family = "gaussian",
      method = loss_chosen,
      verbose = FALSE
    )
    
  }
  
  else if (model_type %in% c("Quantile G-Computation")){
    
    data = data %>% 
      dplyr::select(M1:Y) %>% 
      as.data.frame()
    
    exposures = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
    
    # Fit g-computation model with all two-way interactions (including self-interactions, i.e. quadratic terms)
    fit = qgcomp.boot(
      Y ~ . + .^2,
      data = data,
      family = gaussian(),
      expnms = exposures,
      q = 4, # Number of quantiles to create indicators representing exposure vars (can increase if we want)
      B = 200, # Number of bootstrap iterations
      degree = 2, # Polynomial bases for marginal model
      seed = 2132
    )
    
    # https://github.com/alexpkeil1/qgcomp/blob/main/vignettes/qgcomp-vignette.Rmd 
    # Note: Why don't I get weights from the `boot` functions?
    # Users often use the `qgcomp.*.boot` functions because the want to marginalize over confounders or fit a non-linear joint exposure function. In both cases, the overall exposure response will no longer correspond to a simple weighted average of model coefficients, so none of the `qgcomp.*.boot` functions will calculate weights. In most use cases, the weights would vary according to which level of joint exposure you're at, so it is not a straightforward proposition to calculate them (and you may not wish to report 4 sets of weights if you use the default `q=4`). If you fit an otherwise linear model, you can get weights from a `qgcomp.*.noboot` which will be very close to the weights you might get from a linear model fit via `qgcomp.*.boot` functions, but be explicit that the weights come from a different model than the inference about joint exposure effects.
    
  }
  
  else if (model_type %in% c("Causal Forest")){
    
    # Causal random forest (`grf`) using regression_forest, which trains a regression forest that can be used to estimate the conditional mean function mu(x) = E[Y | X = x]
    # https://grf-labs.github.io/grf/REFERENCE.html#:~:text=For%20regression%20forests%2C%20the%20prediction,status%20of%20the%20neighbor%20examples.
    # Note: should I use quantile_forest instead?
    
    x_train = data %>% 
      dplyr::select(M1:C5) %>% 
      as.matrix()
    
    y_train = data %>% 
      dplyr::select(Y) %>% 
      as.matrix()
    
    # Fit regression forest model using mostly defaults
    fit = regression_forest(
      X = x_train,
      Y = y_train,
      seed = 2132
      
    )
    
  }
  
  else if (model_type %in% c("AMC G-Computation")){
    
    data = data %>% 
      distinct() # To make it work with bootstrapping (other contrasts issue re: number of levels)
    
    x_train = data %>% 
      dplyr::select(M1:C5)
    
    amc_result = amc(
      data = x_train
    )
    
    output = amc_result$output
    
    oldseed = .Random.seed
    
    output = cbind(data %>% dplyr::select(Y), output)
    
    exposures = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
    
    # Obtain joint mixture effect using same method as Quantile-based G-Computation
    fit = qgcomp.boot(
      Y ~ . + .^2,
      data = output,
      family = gaussian(),
      expnms = exposures,
      q = 4, # Number of quantiles to create indicators representing exposure vars (can increase if we want)
      B = 200, # Number of bootstrap iterations
      degree = 2, # Polynomial bases for marginal model
      seed = 2132
    )
    
  }
  
  else if (model_type %in% c("BKMR1")){
    
    data = data %>% 
      distinct()
    
    covars = data %>% 
      dplyr::select(C1:C5) %>% 
      as.matrix()
    
    exposures = data %>% 
      dplyr::select(M1:M10) %>% 
      as.matrix()
    
    y_train = data %>% 
      dplyr::select(Y) %>% 
      as.matrix()
    
    # Use knots to improve speed
    knots100 = fields::cover.design(
      exposures,
      nd = 100
    )$design
    
    # Fit BKMR using MCMC
    fit = kmbayes(
      y = y_train,
      Z = exposures,
      X = covars,
      iter = 1000, # Change to 10000; run time was too long
      varsel = TRUE, # Variable selection,
      est.h = TRUE
    )
    
  }
  
  else if (model_type %in% c("BKMR2")){
    
    data = data %>% 
      distinct()
    
    exposures_and_covars = data %>% 
      dplyr::select(M1:C5) %>% 
      as.matrix()
    
    y_train = data %>% 
      dplyr::select(Y) %>% 
      as.matrix()
    
    # Use knots to improve speed
    knots100 = fields::cover.design(
      exposures_and_covars,
      nd = 100
    )$design
    
    # Fit BKMR using MCMC
    fit = kmbayes(
      y = y_train,
      Z = exposures_and_covars,
      iter = 1000, # Change to 10000; run time was too long
      varsel = TRUE, # Variable selection,
      est.h = TRUE
    )
    
  }
  
  return(fit)
  
}

model_fnct = function(scenarios_start, scenarios_end, metals_count, boot_count, model_type){
  
  estimates = data.frame()
  
  for (i in scenarios_start:scenarios_end){
    
    datasets = data.frame()
    
    filtered_datasets = dataset_filter(i)
    
    for (j in c(1:length) + ((currind - 1) * length)){

      t1_indiv = Sys.time()
      
      data = filtered_datasets[[j]] %>% as.data.frame() %>% dplyr::select(M1:Y)
      
      model_fits = model_fit_fnct(model_type, data)
      
      boot_fits = boot_fnct(model_type, data, boot_count)
      
      metals = data.frame()
      
      for (k in 1:metals_count){
        
        if (model_type %in% c("linear", "GAM", "GAM2", "SuperLearner", "Quantile G-Computation", "AMC G-Computation", "BKMR1", "BKMR2")){
          
          # Set Q1 and Q2 for estimand
          q1 = metal_IQR_table(i, k)[1, ] %>% dplyr::select(M1:C5) # Consider reverting back to without curly brackets?
          q2 = metal_IQR_table(i, k)[2, ] %>% dplyr::select(M1:C5) # Consider reverting back to without curly brackets?
          
          # Set Q1 and Q1 for aggregate mixture profile
          q1_mix = metal_IQR_mix_table(i)[1, ] %>% dplyr::select(M1:C5)
          q2_mix = metal_IQR_mix_table(i)[2, ] %>% dplyr::select(M1:C5)
          
          # Apply fitted model to original data to estimate causal effect for individual metals and for mixture profile
          estimated_diff = predict_diff_fnct(model_type, model_fits, q1, q2)[[1]]
          estimated_diff_mix = predict_diff_fnct(model_type, model_fits, q1_mix, q2_mix)[[1]]
          
          # For each bootstrap, predict effect estimate for individual metals and for mixture profile
          boot_estimated_diffs = predict_diff_boots(model_type, boot_fits, q1, q2)
          boot_estimated_diffs_mix = predict_diff_boots(model_type, boot_fits, q1_mix, q2_mix)
          
        }
        
        else if (model_type == "elastic net"){
          
          # Set Q1 and Q2 for estimand
          q1 = metal_IQR_table(i, k)[1, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          q2 = metal_IQR_table(i, k)[2, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          
          # Set Q1 and Q1 for aggregate mixture profile
          q1_mix = metal_IQR_mix_table(i)[1, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          q2_mix = metal_IQR_mix_table(i)[2, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          
          # Apply fitted model to original data to estimate causal effect
          estimated_diff = predict_diff_fnct(model_type, model_fits, q1, q2)
          estimated_diff = estimated_diff[length(estimated_diff)] # Take only last prediction at point of convergence
          estimated_diff_mix = predict_diff_fnct(model_type, model_fits, q1_mix, q2_mix)
          estimated_diff_mix = estimated_diff_mix[length(estimated_diff_mix)]
          
          # For each bootstrap, predict effect estimate
          boot_estimated_diffs = predict_diff_boots(model_type, boot_fits, q1, q2)
          boot_estimated_diffs_mix = predict_diff_boots(model_type, boot_fits, q1_mix, q2_mix)
          
        }
        
        else if (model_type %in% c("MARS", "Random Forest", "Causal Forest")){
          
          # Set Q1 and Q2 for estimand
          q1 = metal_IQR_table(i, k)[1, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          q2 = metal_IQR_table(i, k)[2, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          
          # Set Q1 and Q1 for aggregate mixture profile
          q1_mix = metal_IQR_mix_table(i)[1, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          q2_mix = metal_IQR_mix_table(i)[2, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          
          # Apply fitted model to original data to estimate causal effect
          estimated_diff = predict_diff_fnct(model_type, model_fits, q1, q2)[[1]]
          estimated_diff_mix = predict_diff_fnct(model_type, model_fits, q1_mix, q2_mix)[[1]]
          
          # For each bootstrap, predict effect estimate
          boot_estimated_diffs = predict_diff_boots(model_type, boot_fits, q1, q2)
          boot_estimated_diffs_mix = predict_diff_boots(model_type, boot_fits, q1_mix, q2_mix)
          
        }
        
        else if (model_type == "BART"){
          
          # Set Q1 and Q2 for estimand
          q1 = metal_IQR_table(i, k)[1, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          q2 = metal_IQR_table(i, k)[2, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          
          # Set Q1 and Q1 for aggregate mixture profile
          q1_mix = metal_IQR_mix_table(i)[1, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          q2_mix = metal_IQR_mix_table(i)[2, ] %>% dplyr::select(M1:C5) %>% as.matrix()
          
          # Apply fitted model to original data to estimate causal effect
          # Using median difference for BART; is there a preferable method, e.g. yhat.train.mean or yhat.test.mean?
          estimated_diff = predict_diff_fnct(model_type, model_fits, q1, q2)[[1]]
          estimated_diff_mix = predict_diff_fnct(model_type, model_fits, q1_mix, q2_mix)[[1]]
          
          # For each bootstrap, predict effect estimate
          boot_estimated_diffs = predict_diff_boots(model_type, boot_fits, q1, q2)
          boot_estimated_diffs_mix = predict_diff_boots(model_type, boot_fits, q1_mix, q2_mix)
          
        }
        
        # Find quantiles for confidence intervals
        ci_ll = quantile(boot_estimated_diffs %>% as.matrix(), 0.025)
        ci_ul = quantile(boot_estimated_diffs %>% as.matrix(), 0.975)
        ci_ll_mix = quantile(boot_estimated_diffs_mix %>% as.matrix(), 0.025)
        ci_ul_mix = quantile(boot_estimated_diffs_mix %>% as.matrix(), 0.975)
        
        # Find variance
        variance_indiv = var(boot_estimated_diffs)
        variance_mix = var(boot_estimated_diffs_mix)

      	# Calculate time difference
      	t2_indiv = Sys.time()
        
        row = bind_cols(
          scenario = i,
          dataset = j,
          metal = k,
          point_estimate = estimated_diff,
          variance = variance_indiv,
          CI_lower = ci_ll,
          CI_upper = ci_ul
        )
        
        metals = bind_rows(metals, row)
        
      }
      
      row_mix = bind_cols(
        scenario = i,
        dataset = j,
        metal = 11, # Must be numeric; consider metal 11 as aggregate mixture
        point_estimate = estimated_diff_mix,
        variance = variance_mix,
        CI_lower = ci_ll_mix,
        CI_upper = ci_ul_mix
      )
      
      metals = bind_rows(metals, row_mix)
  
      datasets = bind_rows(datasets, metals)
      
      datasets$time = (t2_indiv - t1_indiv)[1]
      
    }

    estimates = bind_rows(estimates, datasets)
    
  }
  
  output = cbind(estimates, model_type)
  return(output)
  
}

# Run model on a few scenarios with 100 bootstraps
t1 = Sys.time()
gam_df = model_fnct(scenarios_start = 8, scenarios_end = 8, metals_count = 10, boot_count = 100, model_type = "GAM")
t2 = Sys.time()

# Output file
write.csv(gam_df, file = paste0("gam1_estimates_scen8_", currind, ".csv"))

# Check time to run script
time_mins = t2 - t1
time_mins
