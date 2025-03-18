
#' The two predictive knockoff filters we presented in Section 2.3
#'
#'
#' @param X data.frame (or tibble) the input X, nrow(X) are the number of examples and ncol(X) the number of variables (features)
#' @param X_tilde data.frame (or tibble) the knockoffs of X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial").
#' @param t treatment allocation binary vector with \code{length(t) = nrow(X)}. 
#' @param fdr_nominal target false discovery rate. Can be a vector of multiple thresholds.
#' @param family should be "gaussian" if y is numeric, but "binomial" if y is a binary factor variable.
#'
#' @return a vector of selected indices

# linear interaction KO filter (presented in Section 2.3.1)
linear_model_predictive_filter = function(X, X_tilde, y, t, fdr_nominal, family){
  
  p = dim(X)[2]
  
  # Combine the original Xs and their interactions T:X
  X_original_interactions = matrix(0, dim(X)[1], dim(X)[2])
  for (col in 1:p){
    X_original_interactions[,col]=  X[,col]*t
  }
  
  # Generate interactions with knockoffs the knockoffs Xs and their interactions T:X
  X_knockoffs_interactions = matrix(0, dim(X)[1], dim(X)[2])
  for (col in 1:p){
    X_knockoffs_interactions[,col] =  X_tilde[,col]*t
  }
  
  # Combine data
  X_combined = cbind(X, X_tilde, X_original_interactions, X_knockoffs_interactions, t)
  
  
  # Derive the variable importances
  X_combined = as.matrix(X_combined)
  y = as.matrix(y)
  models_glmnet_cv = glmnet::cv.glmnet(X_combined, y, family = family, nfolds=10, alpha = 1.0,  standardize=TRUE)
  
  importance_scores = coef(models_glmnet_cv, s = "lambda.1se")[2:(dim(X_combined)[2]+1)] # Ignore intercept
  
  # Derive W statistic only for the interaction terms
  W = abs(importance_scores[(2*p+1):(3*p)]) - abs(importance_scores[(3*p+1):(4*p)])
  
  # Calculate the threshold
  threshold = knockoff::knockoff.threshold(W, fdr=fdr_nominal)
  # Find the selected features
  selected_features = which(W >= threshold)
  
  return(selected_features)
  
}

causal_forest_predictive_filter = function(X, X_tilde, y, t, fdr_nominal, family = NULL){
  print('check1')
  p = dim(X)[2]
  print('check2')
  
  # Build causal forests using initial features and knockoffs
  c.forest = grf::causal_forest(cbind(X, X_tilde), y, t, mtry = 2*p)
  print('check3')
  # Derive the variable importance scores for initial featurs and knockoffs
  # We are using the default parameters for decay.exponent and max.depth
  importance_scores = grf::variable_importance(c.forest, decay.exponent = 2, max.depth = 5)
  print('check4')
  # Calculate the W statistic as the difference
  W = abs(importance_scores[1:p]) - abs(importance_scores[(p+1):(2*p)])
  print('check5')
  # Calculate the threshold
  threshold = knockoff::knockoff.threshold(W, fdr=fdr_nominal)
  print('check5')
  # Find the selected features
  selected_features = which(W >= threshold)
  print('check6')
  
  # return(selected_features)
  return(selected_features)
}









#' Simulates the scenarios of Section 3
#'
#'
#' @param scenario an integer that represents the scenario we explore in the paper, eg. scenario=1, simulates the scenario of Section 3.4.1, scenario=2, simulates the scenario of Section 3.4.2 etc
#' @param model an integer that represents the model from which we simulate the data, eg model=1, model=2 etc 
#' @param predictive_amplitude a number that captures the strength of the predictive signal, e.g. the parameter theta_pred of the paper
#' @param sample_size the number of exmamples (e.g. patients)
#' @param num_features the number of variables (features, covarietes)
#' 
#' @return the simulated data 
generate_scenarios_predictive = function(scenario, model, sample_size, num_features, predictive_amplitude, predictive_amplitude_high=-1, predictive_amplitude_mod=-1, predictive_amplitude_low=-1, first_split=-1, second_split=-1) {
  
  
  if(scenario == 'S1')
  {
    # Sample T:
    prob_t = 0.50 
    t = rbinom(sample_size,1,prob_t)
    # Sample Xs
    rho = 0.50
    sigma_z =  toeplitz(rho^(0:(num_features - 1)))
    X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
    # Generate y
    if(model=='M1')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(11,20,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M2')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(1,10,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M3')
    {
      nonzero_pred = seq(1,5,1)
      nonzero_prog = seq(1,10,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M4')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(1,5,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M5')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = 0
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M6') # we are considering high-dimensional scenarios. in the first case, S_pred = S_prog
    {
      nonzero_pred = seq(1,100,1)
      nonzero_prog = seq(1,100,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  
  
  if(scenario == 'S2')
  {
    # Sample t
    prob_t = 0.50 
    t = rbinom(sample_size,1,prob_t)
    # Sample Xs
    rho = 0.50
    sigma_z =  toeplitz(rho^(0:(num_features - 1)))
    X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
    # Generate y 
    if(model=='M1')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = (apply(exp(matrix(runif(sample_size*length(nonzero_prog),0.50,1.50), sample_size, length(nonzero_prog))*X[,nonzero_prog]),1,sum))
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude !=0){
        pred_part = predictive_amplitude*t*(apply(exp(kronecker(t(runif(length(nonzero_pred),0.50,1.50)), matrix(1, sample_size, 1))*X[,nonzero_pred]),1,sum))
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if(model=='M2')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = X[,nonzero_prog[1]] * X[,nonzero_prog[2]] + X[,nonzero_prog[3]]*X[,nonzero_prog[4]] + X[,nonzero_prog[5]]*X[,nonzero_prog[6]]+ X[,nonzero_prog[7]]*X[,nonzero_prog[8]] + X[,nonzero_prog[9]]*X[,nonzero_prog[10]] 
      if(predictive_amplitude == 0){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude != 0 ){
        pred_part = X[,nonzero_pred[1]] * X[,nonzero_pred[2]] + X[,nonzero_pred[3]]*X[,nonzero_pred[4]] + X[,nonzero_pred[5]]*X[,nonzero_pred[6]]+ X[,nonzero_pred[7]]*X[,nonzero_pred[8]] + X[,nonzero_pred[9]]*X[,nonzero_pred[10]]
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if(model=='M3')
    {
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = X[,nonzero_prog[1]] + X[,nonzero_prog[2]]^2 + X[,nonzero_prog[3]]^3 + X[,nonzero_prog[4]]*X[,nonzero_prog[5]] + X[,nonzero_prog[6]]*X[,nonzero_prog[7]] + X[,nonzero_prog[8]]*X[,nonzero_prog[9]]*X[,nonzero_prog[10]]
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude != 0){
        pred_part = X[,nonzero_pred[1]] + X[,nonzero_pred[2]]^2 + X[,nonzero_pred[3]]^3 + X[,nonzero_pred[4]]*X[,nonzero_pred[5]] + X[,nonzero_pred[6]]*X[,nonzero_pred[7]] + X[,nonzero_pred[8]]*X[,nonzero_pred[9]]*X[,nonzero_pred[10]]
        pred_part = predictive_amplitude*t*pred_part
        mu =  0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  
  
  if(scenario == 'S3')
  {
    # Sample t
    prob_t = 0.50 
    t = rbinom(sample_size,1,prob_t)
    # Sample Xs
    rho = 0.50
    sigma_z =  toeplitz(rho^(0:(num_features - 1)))
    X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
    # Generate y
    if(model=='M1')
    {
      nonzero_pred = seq(1,5,1)
      nonzero_prog = seq(1,5,1)
      prog_part =   X[,nonzero_prog[1]] + X[,nonzero_prog[2]] + X[,nonzero_prog[3]] + X[,nonzero_prog[4]] + X[,nonzero_prog[5]]
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      
      if(predictive_amplitude != 0){
        pred_part = (X[,nonzero_pred[1]]>0) + (X[,nonzero_pred[2]]<0) + (X[,nonzero_pred[3]]<0.675) + (X[,nonzero_pred[4]]>0.675) + (X[,nonzero_pred[5]]<0.675)
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    
    if(model=='M2')
    {
      nonzero_pred = seq(1,5,1)
      nonzero_prog = seq(1,5,1)
      prog_part =   X[,nonzero_prog[1]] + X[,nonzero_prog[2]] + X[,nonzero_prog[3]] + X[,nonzero_prog[4]] + X[,nonzero_prog[5]]
      if(predictive_amplitude == 0){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      if(predictive_amplitude != 0){
        pred_part = abs(X[,nonzero_pred[1]])* (X[,nonzero_pred[2]]>0.675) + abs(X[,nonzero_pred[3]])*(X[,nonzero_pred[4]]<0)*(X[,nonzero_pred[5]]>0)
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    if(model=='M3')
    {
      nonzero_pred_high = seq(76,100,1)
      nonzero_pred_mod = seq(26,75,1)
      nonzero_pred_low = seq(1,25,1)
      nonzero_prog = seq(1,100,1)
      nonzero_pred = seq(1,100,1)
      prog_part =   X[,nonzero_prog[1]] + X[,nonzero_prog[2]] + X[,nonzero_prog[3]] + X[,nonzero_prog[4]] + X[,nonzero_prog[5]]
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      # this line doesnt actually do anything here
      if(predictive_amplitude != 0){
        # pred_part = (X[,nonzero_pred[1]]>0) + (X[,nonzero_pred[2]]<0) + (X[,nonzero_pred[3]]<0.675) + (X[,nonzero_pred[4]]>0.675) + (X[,nonzero_pred[5]]<0.675)
        pred_part_high <- predictive_amplitude_high * t * rowSums(X[, nonzero_pred_high])
        pred_part_mod  <- predictive_amplitude_mod  * t * rowSums(X[, nonzero_pred_mod])
        pred_part_low  <- predictive_amplitude_low  * t * rowSums(X[, nonzero_pred_low])
        pred_part = pred_part_high + pred_part_mod + pred_part_low
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
     if(model=='M3_dash')
    {
      nonzero_pred_high = seq(second_split+1,100,1)
      nonzero_pred_mod = seq(first_split+1,second_split,1)
      nonzero_pred_low = seq(1,first_split,1)
      nonzero_prog = seq(1,100,1)
      nonzero_pred = seq(1,100,1)
      prog_part =   X[,nonzero_prog[1]] + X[,nonzero_prog[2]] + X[,nonzero_prog[3]] + X[,nonzero_prog[4]] + X[,nonzero_prog[5]]
      if(predictive_amplitude == 0 ){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      # this line doesnt actually do anything here
      if(predictive_amplitude != 0){
        # pred_part = (X[,nonzero_pred[1]]>0) + (X[,nonzero_pred[2]]<0) + (X[,nonzero_pred[3]]<0.675) + (X[,nonzero_pred[4]]>0.675) + (X[,nonzero_pred[5]]<0.675)
        pred_part_high <- predictive_amplitude_high * t * rowSums(X[, nonzero_pred_high])
        pred_part_mod  <- predictive_amplitude_mod  * t * rowSums(X[, nonzero_pred_mod])
        pred_part_low  <- predictive_amplitude_low  * t * rowSums(X[, nonzero_pred_low])
        pred_part = pred_part_high + pred_part_mod + pred_part_low
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    if (predictive_amplitude == 0) 
      {nonzero_pred = 0}
  }
  
  
  if(scenario == 'S4')
  {
    
    if(model=='M1')
    { 
      
      # Sample t
      prob_t = 0.50 
      t = rbinom(sample_size,1,prob_t)
      # Sample Xs
      rho = 0.50
      sigma_z =  toeplitz(rho^(0:(num_features - 1)))
      X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
      # Generate y
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(11,20,1)
      prog_part = apply(X[,nonzero_prog],1,sum)
      pred_part = predictive_amplitude*t*(apply(X[,nonzero_pred],1,sum))
      mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
      y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
    }
    
    if(model=='M2')
    {     
      
      # Sample t
      prob_t = 0.50 
      t = rbinom(sample_size,1,prob_t)
      # Sample Xs
      rho = 0.50
      sigma_z =  toeplitz(rho^(0:(num_features - 1)))
      X = data.frame(mvtnorm::rmvnorm(sample_size, rep(0, num_features), sigma_z, "chol"))
      # Generate y
      nonzero_pred = seq(1,10,1)
      nonzero_prog = seq(6,15,1)
      prog_part = X[,nonzero_prog[1]] * X[,nonzero_prog[2]] + X[,nonzero_prog[3]]*X[,nonzero_prog[4]] + X[,nonzero_prog[5]]*X[,nonzero_prog[6]]+ X[,nonzero_prog[7]]*X[,nonzero_prog[8]] + X[,nonzero_prog[9]]*X[,nonzero_prog[10]] 
      if(predictive_amplitude == 0){
        mu = 0.5*t + as.numeric(as.matrix((prog_part)))
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
      if(predictive_amplitude != 0 ){
        pred_part = X[,nonzero_pred[1]] * X[,nonzero_pred[2]] + X[,nonzero_pred[3]]*X[,nonzero_pred[4]] + X[,nonzero_pred[5]]*X[,nonzero_pred[6]]+ X[,nonzero_pred[7]]*X[,nonzero_pred[8]] + X[,nonzero_pred[9]]*X[,nonzero_pred[10]] 
        pred_part = predictive_amplitude*t*pred_part
        mu = 0.5*t + as.numeric(as.matrix((prog_part)) + as.matrix((pred_part)) )
        y = mu + matrix(rnorm(sample_size, mean = 0, sd = 1), sample_size, 1)
      }
    }
    if (predictive_amplitude == 0) {nonzero_pred = 0}
  }
  
  output = c()
  output$t = t
  output$X = X
  output$y = y[,1]
  
  output$predictive = nonzero_pred
  
  
  return(output)
}

# Function to compute False Discovery Proportion (FDP)
fdp <- function(selected_features, true_predictive_features) {
  if (length(selected_features) == 0) {
    return(0)  # If no features are selected, FDP is 0 by convention
  }
  
  false_positives = setdiff(selected_features, true_predictive_features)  # Features that were wrongly selected
  fdp = length(false_positives) / length(selected_features)  # FDP formula
  
  return(fdp)
}

# Function to compute Power
power <- function(selected_features, true_predictive_features) {
  if (length(true_predictive_features) == 0) {
    return(NA)  # Avoid division by zero if no predictive features exist
  }
  
  true_positives = intersect(selected_features, true_predictive_features)  # Correctly selected predictive features
  power = length(true_positives) / length(true_predictive_features)  # Power formula
  
  return(power)
}
