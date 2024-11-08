#!/usr/bin/R
## Senan Hogan-Hennessy, 04 October 2023
## Define the CCI estimator, in a function.

################################################################################
## Define the CCI estimation weights.

cci_weights.est <- function(Y, D, Z, X, count.trees = 2000){
    # Requires generalised Random Forests, https://grf-labs.github.io/grf/
    library(grf)
    # @param Y is an outcome, type vector, continuous
    # @param D is the treatment, type vector, binary
    # @param Z is the instrument, type vector, binary
    # @param X is the characteritics, type matrix, continuous.
    
    # Estimate the total conditional effect
    totaleffect.forest <- causal_forest(X, Y, Z, num.trees = count.trees)
    
    ## Non-parametrically estimate inst propensity + compliance.
    # Estimate Pr(Z = 1 | X) instrument propensity
    instpropensity.forest <- probability_forest(X, as.factor(Z),
        num.trees = count.trees)
    # Estimate Pr(complier | X) = E[D | X, Z = 1] - E[D | X, Z = 0] by honest GRF
    complier.forest <- causal_forest(X, D, Z, num.trees = count.trees)
    nevertaker.forest <- probability_forest(X[Z == 1,],
        as.factor(1 - D[Z == 1]), num.trees = count.trees)
    alwaystaker.forest <- probability_forest(X[Z == 0,],
        as.factor(D[Z == 0]), num.trees = count.trees)
    # Estimate E[Y | X, Z, D] by honest GRF
    Y.forest <- boosted_regression_forest(
        cbind(Z, D, X), Y, num.trees = count.trees)
        # Predict across X
    totaleffect.prediction <- predict(totaleffect.forest, X)$predictions
    inst.propensity <- predict(instpropensity.forest, X)$predictions[, 2]
    complier.score <- predict(complier.forest, X)$predictions
    nevertaker.score <- predict(nevertaker.forest, X)$predictions[, 2]
    alwaystaker.score <- predict(alwaystaker.forest, X)$predictions[, 2]
    hat_Y_1_1 <- predict(Y.forest, newdata = cbind(1, 1, X))$predictions
    hat_Y_1_0 <- predict(Y.forest, newdata = cbind(1, 0, X))$predictions
    hat_Y_0_1 <- predict(Y.forest, newdata = cbind(0, 1, X))$predictions
    hat_Y_0_0 <- predict(Y.forest, newdata = cbind(0, 0, X))$predictions
    hat_Y_1 <- predict(Y.forest, newdata = cbind(1, D, X))$predictions
    hat_Y_0 <- predict(Y.forest, newdata = cbind(0, D, X))$predictions
    # Define the object to return.
    cci.weights <- data.frame(
        totaleffect = totaleffect.prediction,
        inst.propensity = inst.propensity,
        complier.score = complier.score,
        nevertaker.score = nevertaker.score,
        alwaystaker.score = alwaystaker.score,
        hat_Y_1_1 = hat_Y_1_1,
        hat_Y_1_0 = hat_Y_1_0,
        hat_Y_0_1 = hat_Y_0_1,
        hat_Y_0_0 = hat_Y_0_0,
        hat_Y_1 = hat_Y_1,
        hat_Y_0 = hat_Y_0)
    # Return the named object, dataframe of weights.
    return(cci.weights)
}


################################################################################
## Define the CCI estimator.

cci_point.est <- function(Y, D, Z, X, cci.weights, subset = NULL){
    # If specified, subset data.  CCI weights pre-calculated on entire data.
    if (!is.null(subset)){
        cci.weights <- cci.weights[subset, ]
        Y <- Y[subset]
        D <- D[subset]
        Z <- Z[subset]
        X <- X[subset, ]
    }
    # Separately get the scores + pre-computed k-weights
    totaleffect.est <- cci.weights$totaleffect
    inst.propensity <- cci.weights$inst.propensity
    complier.score <- cci.weights$complier.score
    nevertaker.score <- cci.weights$nevertaker.score
    alwaystaker.score <- cci.weights$alwaystaker.score
    hat_Y_1_1 <- cci.weights$hat_Y_1_1
    hat_Y_1_0 <- cci.weights$hat_Y_1_0
    hat_Y_0_1 <- cci.weights$hat_Y_0_1
    hat_Y_0_0 <- cci.weights$hat_Y_0_0
    hat_Y_1 <- cci.weights$hat_Y_1
    hat_Y_0 <- cci.weights$hat_Y_0
    # First estimate the reduced form total effect, E[Y | X, Z = 0] - E[Y | X, Z = 0]
    total_effect.est <- totaleffect.est
    # Resulting average mechanism effect estimator
    mechanism_effect.est <- (inst.propensity * (hat_Y_1_1 - hat_Y_1_0) +
        (1 - inst.propensity) * (hat_Y_0_1 - hat_Y_0_0))
    ame_lower.est <- (
        inst.propensity * (1 - nevertaker.score) * (hat_Y_1_1 - hat_Y_1) +
        (1 - inst.propensity) * (1 - alwaystaker.score) * (hat_Y_0 - hat_Y_0_0))
    # Compare the effect estimates
    # print(c("Average mechanism effect D -> Y:", mean(mechanism_effect.est)))
    # print(c("LAME, D -> Y:", mean(lame.est)))
    # Resulting direct effect estimator.
    direct_effect.est <- total_effect.est - (complier.score * mechanism_effect.est)
    # print(c("Average direct effect Z -> Y:", mean(direct_effect.est)))
    # print(c("Average indirect effect D(Z) -> Y:",
    #     mean(complier.score * lame.est)))
    # Compute the bias in IV estimation, Wald est minus local mechanism effect
    iv.est <- totaleffect.est / complier.score
    # Naive OLS est.
    ols.est <- coef(lm(Y ~ 1 + D + X))["D"]
    # print(c("naive IV est of D -> Y:", mean(iv.est)))
    # Define the object to return.
    return.list <- c(
        mean(total_effect.est),
        mean(mechanism_effect.est),
        weighted.mean(mechanism_effect.est, complier.score),
        mean(ame_lower.est),
        mean(direct_effect.est),
        mean(complier.score * mechanism_effect.est),
        mean(iv.est),
        ols.est)
    return.names <- c(
        "total_effect",
        "mechanism_effect",
        "LAME",
        "mechanism_lower",
        "direct_effect",
        "indirect_effect",
        "iv_est",
        "ols_est")
    return.list <- setNames(return.list, return.names)
    # Return the named object, point estimates.
    return(return.list)
}

#TODO here:
#TODO 1. Bootstrap the above
#TODO 2. Code a figure for est: OLS, naive IV, CCI.
#TODO 3. calculate the conditional bounds.


################################################################################
## Bootstrap the compliance weighting CCI estimator.

cci_bootstrap.est <- function(Y, D, Z, X, count.trees = 2000, boot.count = 1000){
    # @param Y is an outcome, type vector, continuous
    # @param D is the treatment, type vector, binary
    # @param Z is the instrument, type vector, binary
    # @param X is the characteritics, type matrix, continuous.
    # Define lists to return.
    total_effect.list <- c()
    LAME.list <- c()
    mechanism_effect.list <- c()
    mechanism_lower.list <- c()
    direct_effect.list <- c()
    indirect_effect.list <- c()
    iv_est.list <- c()
    ols_est.list <- c()
    # Perform each bootstrap calculation:
    for (i in seq(1, boot.count)){
        print(paste0("Calculating boot estimate ",
            i - 1, " out of ", boot.count, ", progress: ", (i - 1) / boot.count))
        # Bootstrap the data.
        boot.indices <- sample(seq(1, length(Y)), length(Y), replace = TRUE)
        Y.boot <- Y[boot.indices]
        D.boot <- D[boot.indices]
        Z.boot <- Z[boot.indices]
        X.boot <- X[boot.indices, ]
        # Calculte the boot estimates.
        cci_boot.weights <- cci_weights.est(
            Y.boot, D.boot, Z.boot, X.boot, count.trees = count.trees)
        cci_boot.est <- cci_point.est(
            Y.boot, D.boot, Z.boot, X.boot, cci_boot.weights)
        # Store the boot estimates.
        total_effect.list <- c(total_effect.list,
            cci_boot.est["total_effect"])
        LAME.list <- c(LAME.list,
            cci_boot.est["LAME"])
        mechanism_effect.list <- c(mechanism_effect.list,
            cci_boot.est["mechanism_effect"])
        mechanism_lower.list <- c(mechanism_lower.list,
            cci_boot.est["mechanism_lower"])
        direct_effect.list <- c(direct_effect.list,
            cci_boot.est["direct_effect"])
        indirect_effect.list <- c(indirect_effect.list,
            cci_boot.est["indirect_effect"])
        iv_est.list <- c(iv_est.list,
            cci_boot.est["iv_est"])
        ols_est.list <- c(ols_est.list,
            cci_boot.est["ols_est"])
        # Save the memory, thanks to R no auto-GC
        rm(cci_boot.weights, cci_boot.est)
        gc()
    }
    # Define the data to return
    return.data <- data.frame(
        total_effect.list,
        LAME.list,
        mechanism_effect.list,
        mechanism_lower.list,
        direct_effect.list,
        indirect_effect.list,
        iv_est.list,
        ols_est.list)
    names(return.data) <- c(
        "total_effect",
        "LAME",
        "mechanism_effect",
        "mechanism_lower",
        "direct_effect",
        "indirect_effect",
        "iv_est",
        "ols_est")
    # Return a dataframe of the bootstrap estimates.
    return(return.data)
}
