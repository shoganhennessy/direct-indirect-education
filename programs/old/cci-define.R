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

    ## Non-parametrically estimate inst propensity + compliance.
    # Estimate Pr(Z = 1 | X) instrument propensity
    instpropensity.forest <- probability_forest(X, as.factor(Z),
        num.trees = count.trees)
    # Estimate Pr(AT | X) = E[D | X, Z = 0] by honest GRF
    alwaystaker.forest <- probability_forest(X[Z == 0, ], as.factor(D[Z == 0]),
        num.trees = count.trees)
    # Estimate Pr(NT | X) = 1 - E[D | X, Z = 1] by honest GRF
    nevertaker.forest <- probability_forest(X[Z == 1, ], as.factor(1 - D[Z == 1]),
        num.trees = count.trees)
    # Predict across X
    inst.propensity <- predict(instpropensity.forest, X)$predictions[, 2]
    alwaystaker.score <- predict(alwaystaker.forest, X)$predictions[, 2]
    nevertaker.score <- predict(nevertaker.forest, X)$predictions[, 2]
    complier.score <- (1 - nevertaker.score - alwaystaker.score)
    # Define the k_{1,1} weight for Y(1,1), i.e. Y_i | Z_i = D_i = 1
    k_1_1 <- (Z * D) / (complier.score + alwaystaker.score)
    # Define the k_{1,0} weight for Y(1,0), i.e. Y_i | Z_i = 1, D_i = 0
    k_1_0 <- (Z * (1 - D)) / (nevertaker.score * inst.propensity)
    # Define the k_{0,1} weight for Y(0,1), i.e. Y_i | Z_i = 0, = D_i = 1
    k_0_1 <- ((1 - Z) * D) / (alwaystaker.score * (1 - inst.propensity))
    # Define the k_{0,0} weight for Y(0,0), i.e. Y_i | Z_i = 0, D_i = 0
    k_0_0 <- ((1 - Z) * (1 - D)) / (complier.score + nevertaker.score)
    # Define the object to return.
    cci.weights <- data.frame(
        inst.propensity = inst.propensity,
        alwaystaker.score = alwaystaker.score,
        nevertaker.score = nevertaker.score,
        complier.score = complier.score,
        k_1_1 = k_1_1,
        k_1_0 = k_1_0,
        k_0_1 = k_0_1,
        k_0_0 = k_0_0)
    # Return the named object, dataframe of weights.
    return(cci.weights)
}


################################################################################
## Define the direct compliance weighting CCI estimator.

cci_direct.est <- function(Y, D, Z, X, cci.weights, subset = NULL){
    # If specified, subset data.  CCI weights pre-calculated on entire data.
    if (!is.null(subset)){
        cci.weights <- cci.weights[subset, ]
        Y <- Y[subset]
        D <- D[subset]
        Z <- Z[subset]
        X <- X[subset, ]
    }
    # Separately get the scores + pre-computed k-weights
    inst.propensity <- cci.weights$inst.propensity
    alwaystaker.score <- cci.weights$alwaystaker.score
    nevertaker.score <- cci.weights$nevertaker.score
    complier.score <- cci.weights$complier.score
    k_1_1 <- cci.weights$k_1_1
    k_1_0 <- cci.weights$k_1_0
    k_0_1 <- cci.weights$k_0_1
    k_0_0 <- cci.weights$k_0_0
    # First estimate the reduced form total effect, E[Y | Z = 0] - E[Y | Z = 0]
    total_effect.est <- coef(lm(Y ~ 1 + Z + X))["Z"]
    # Resulting average mechanism effect estimator, from weight estimator.
    mechanism_effect.est <- mean(inst.propensity) * (
            weighted.mean(Y, k_1_1) - weighted.mean(Y, k_1_0)) +
        mean(1 - inst.propensity) * (
            weighted.mean(Y, k_0_1) - weighted.mean(Y, k_0_0))
    # Show the effect estimates
    print(c("Average mechanism effect D -> Y:", mean(mechanism_effect.est)))
    # Resulting direct effect estimator, averaged across compliance groups.
    direct_effect.est <- mean(complier.score) * (mean(inst.propensity) * (
            weighted.mean(Y, k_1_1) - weighted.mean(Y, k_0_1)) +
        mean(1 - inst.propensity) * (
            weighted.mean(Y, k_1_0) - weighted.mean(Y, k_0_0))) +
        mean(alwaystaker.score) * (
            weighted.mean(Y, k_1_1) - weighted.mean(Y, k_0_1)) +
        mean(nevertaker.score) * (
            weighted.mean(Y, k_1_0) - weighted.mean(Y, k_0_0))
    print(c("Average direct effect Z -> Y:", direct_effect.est))
    # Compute the bias in IV estimation, Wald est minus local mechanism effect
    wald.est <- (mean(Y[Z == 1]) - mean(Y[Z == 0])) / (
        mean(D[Z == 1]) - mean(D[Z == 0]))
    iv_bias.est <- wald.est - mechanism_effect.est
    print(c("Bias in IV est of LAME D -> Y:", iv_bias.est))

    # Define the object to return.
    return.list <- c(mechanism_effect.est, direct_effect.est,
        total_effect.est, iv_bias.est)
    return.names <- c("mechanism_effect", "direct_effect",
        "total_effect", "iv_bias")
    return.list <- setNames(return.list, return.names)
    # Return the named object, point estimates.
    return(return.list)
}

#TODO, write a function here for bootstrapping the above, for SEs.
