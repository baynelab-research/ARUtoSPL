library(dplyr)
library(tidyr)
library(stringr)
library(nimble)
library(coda)
# Removed parallel library as it is no longer required

setwd("~/Desktop/cali/csvs/")
load("wontbreak.Rdata")
set.seed(123)

#############################################################################
## 7) MCMC SETTINGS
#############################################################################
nChains <- 2
nIter   <- 25000
burnin  <- 5000
thin    <- 5
K <- 4

##############
mvCode <- nimbleCode({

  alpha ~ dnorm(0, 0.001)

  beta ~ dnorm(0, 0.001)

  beta_temp    ~ dnorm(0, 0.001)
  beta_wind_sp ~ dnorm(0, 0.001)
  beta_rh      ~ dnorm(0, 0.001)

  for(mm in 1:max_model) {
    Omega[mm, 1:32, 1:32] ~ dwish(R[1:32, 1:32], df)
  }

  for(i in 1:n) {
    for(j in 1:32) {
      mu[i,j] <- alpha +
                 beta * aruvalues[i,j] +
                 beta_temp * temp[i] +
                 beta_wind_sp * wind_sp[i] +
                 beta_rh * rh[i]
    }

    y[i,1:32] ~ dmnorm(mu[i,1:32], Omega[model[i], 1:32, 1:32])
  }

})

#############################################################################
## 6) CROSS-VALIDATION SPLIT
#############################################################################

nAll <- ml$n

fold_ids <- sample(rep(1:K, length.out = nAll))

cv_folds <- vector("list", K)
for(k in seq_len(K)) {
  test_idx  <- which(fold_ids == k)
  train_idx <- setdiff(seq_len(nAll), test_idx)
  cv_folds[[k]] <- list(train_index = train_idx,
                        test_index  = test_idx)
}

#############################################################################
## 8) HELPER FUNCTION: BUILD & RUN NIMBLE MODEL (MULTIPLE CHAINS, SEQUENTIALLY)
#############################################################################
run_nimble_model <- function(code, constants, data,
                             nChains, nIter, burnin, thin,
                             baseSeed = 555) {
  initsFun <- function(chain) {
    inits <- list(
      alpha = 0,
      alpha_model = rep(0, constants$max_model),
      beta_model = rep(0, constants$max_model),
      beta_temp = 0,
      beta_wind_sp = 0,
      beta_rh = 0,
      Omega = array(rep(diag(1, 32), constants$max_model),
                    dim = c(constants$max_model, 32, 32))
    )
    for(mm in 1:constants$max_model) {
      inits$Omega[mm, , ] <- diag(1, 32)
    }
    return(inits)
  }
  
  # 1) Set up the model and compile using chain 1 as a placeholder
  nm <- nimbleModel(code, data = data, constants = constants,
                    inits = initsFun(1))
  Cm <- compileNimble(nm)
  
  # 2) Configure and compile the MCMC
  spec <- buildMCMC(nm)
  Cspec <- compileNimble(spec, project = nm)
  
  # 3) Run all chains sequentially
  chains <- list()
  for(chain in 1:nChains) {
    cat("Running MCMC for chain", chain, "...\n")
    chain_result <- runMCMC(Cspec,
                            inits   = initsFun(chain),
                            niter   = nIter,
                            nburnin = burnin,
                            thin    = thin,
                            setSeed = baseSeed + chain - 1)
    chains[[chain]] <- chain_result
  }
  
  # 4) Convert the results to an mcmc.list for further analysis
  mcmcList <- mcmc.list(lapply(chains, as.mcmc))
  return(mcmcList)
}

#############################################################################
## 9) FINAL MODEL ON ALL DATA
#############################################################################
cat("======== Fitting final model on full data (no missing) ========\n")

constants_full <- list(
  n         = ml$n,
  max_model = max(ml$model),
  R         = diag(1, 32),
  df        = 33,
  model     = ml$model,
  temp      = ml$temp,
  wind_sp   = ml$wind_sp,
  rh        = ml$rh
)

data_full <- list(
  y = ml$y,
  aruvalues = matrix(as.numeric(ml$aruvalues), 
                     nrow = nrow(ml$aruvalues), 
                     ncol = ncol(ml$aruvalues))
)

final_mcmcList <- run_nimble_model(
  code      = mvCode,
  constants = constants_full,
  data      = data_full,
  nChains   = nChains,
  nIter     = nIter,
  burnin    = burnin,
  thin      = thin
)

cat("Final model done.\n")
save.image("full_simp_done.Rdata")

#############################################################################
## Gelman diag (final model)
#############################################################################
chainMatList_final <- lapply(final_mcmcList, as.matrix)
paramNames_final   <- colnames(chainMatList_final[[1]])
gdMat_final <- matrix(NA, nrow = length(paramNames_final), ncol = 2,
                      dimnames = list(paramNames_final, c("PointEst","UpperCI")))
for(p in seq_along(paramNames_final)) {
  pName <- paramNames_final[p]
  chainParamList <- lapply(chainMatList_final, function(x) x[,pName,drop=FALSE])
  paramMcmcList  <- mcmc.list(lapply(chainParamList, as.mcmc))
  gd <- gelman.diag(paramMcmcList, autoburnin=FALSE)
  gdMat_final[p,1] <- gd$psrf[1]
  gdMat_final[p,2] <- gd$psrf[2]
}
cat("\nGelman diag for final model stored in 'gdMat_final'.\n")

#############################################################################
## 10) 5-FOLD CROSS-VALIDATION (updated with weather covariates)
#############################################################################
foldRMSE <- numeric(K)
gelmanResults <- vector("list", K)

for(k in seq_len(K)) {
  cat("\n======== Starting fold", k, "========\n")
  
  train_idx <- cv_folds[[k]]$train_index
  test_idx  <- cv_folds[[k]]$test_index
  
  # Mark test as NA
  y_train <- ml$y
  y_train[test_idx, ] <- NA
  
  constants_cv <- list(
    n         = ml$n,
    max_model = max(ml$model),
    R         = diag(1, 32),
    df        = 33,
    model     = ml$model,
    temp      = ml$temp,
    wind_sp   = ml$wind_sp,
    rh        = ml$rh
  )
  data_cv <- list(
    y         = y_train,
    aruvalues = matrix(as.numeric(ml$aruvalues), 
                       nrow = nrow(ml$aruvalues), 
                       ncol = ncol(ml$aruvalues))
  )
  
  # Fit model
  mcmcList_cv <- run_nimble_model(
    code      = mvCode,
    constants = constants_cv,
    data      = data_cv,
    nChains   = nChains,
    nIter     = nIter,
    burnin    = burnin,
    thin      = thin
  )
  
  # Gelman diag per parameter
  chainMatList_cv <- lapply(mcmcList_cv, as.matrix)
  paramNames_cv <- colnames(chainMatList_cv[[1]])
  
  gdMat_cv <- matrix(NA, nrow = length(paramNames_cv), ncol = 2,
                     dimnames = list(paramNames_cv, c("PointEst","UpperCI")))
  
  for(p in seq_along(paramNames_cv)) {
    pName <- paramNames_cv[p]
    chainParamList <- lapply(chainMatList_cv, function(x) x[,pName,drop=FALSE])
    paramMcmcList  <- mcmc.list(lapply(chainParamList, as.mcmc))
    gd <- gelman.diag(paramMcmcList, autoburnin=FALSE)
    gdMat_cv[p,1] <- gd$psrf[1]
    gdMat_cv[p,2] <- gd$psrf[2]
  }
  gelmanResults[[k]] <- gdMat_cv
  
  # Posterior mean predictions
  postMeans <- colMeans(as.matrix(mcmcList_cv))
  
  alpha_est <- postMeans["alpha"]
  beta_temp_est    <- postMeans["beta_temp"]
  beta_wind_sp_est <- postMeans["beta_wind_sp"]
  beta_rh_est      <- postMeans["beta_rh"]
  
  max_m <- constants_cv$max_model
  alpha_model_est <- numeric(max_m)
  beta_model_est  <- numeric(max_m)
  
  for(mm in seq_len(max_m)) {
    alpha_model_est[mm] <- postMeans[paste0("alpha_model[", mm, "]")]
    beta_model_est[mm]  <- postMeans[paste0("beta_model[", mm, "]")]
  }
  
  # Predict for test
  testY <- ml$y[test_idx, , drop=FALSE]
  nTest <- length(test_idx)
  
  predY <- matrix(NA, nrow=nTest, ncol=32)
  
  for(tt in seq_len(nTest)) {
    iRow <- test_idx[tt]
    mID  <- ml$model[iRow]
    for(j in 1:32) {
      predY[tt,j] <- alpha_est +
                     alpha_model_est[mID] +
                     beta_model_est[mID] * ml$aruvalues[iRow, j] +
                     beta_temp_est * ml$temp[iRow] +
                     beta_wind_sp_est * ml$wind_sp[iRow] +
                     beta_rh_est * ml$rh[iRow]
    }
  }
  
  # RMSE
  diffsq <- (testY - predY)^2
  foldRMSE[k] <- sqrt(mean(diffsq))
  
  cat("Fold", k, "RMSE =", foldRMSE[k], "\n")
}

cat("\nCross-validation RMSE by fold:\n")
print(foldRMSE)
cat("Mean CV RMSE =", mean(foldRMSE), "\n")
save.image("x_sim_done.Rdata")
cat("\nAll done.\n")
