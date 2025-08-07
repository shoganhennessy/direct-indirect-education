
# Senan Hogan-Hennessy, 6 June 2025.
# forked from original code by Hyeokmoon Kweon, h.kweon@vu.nl

######################################################################################
# function GORIV: Runs ORIV 

# Note: the function depends on "fixest", "data.table" packages

# Inputs:
# fm: formula object to specify a model. Do not include a PGI here, but only covariates. 
# PGI1: character object of first PGI name
# PGI2: character object of second PGI name
# data: data.frame or data.table 
# IID: character object of indivdiual ID column name
# FID: character object of family ID column name, Default to NULL. If provided, SE will be clustered at family, too. Has to be supplied for within-family estimation.  
# resid: logical, indicating whether the target outcome should be residualised and standardised first. Default to TRUE
# within: logical, indicating within-family estimation. Default to FALSE

# Load required packages.
library(data.table)
library(fixest)

# Output:
# fixest object from "fixest" package
# see https://www.rdocumentation.org/packages/fixest/versions/0.8.4/topics/feols

GORIV <- function(fm, PGI1, PGI2, IID, data, FID = NULL,
        IV.list = NULL,
        control.list = NULL){
    # Put data into data.table format.
    data <- data.table(data)
    # scale PGI (only do if not using random component IVs).
    # compute scaling factor
    R <- sqrt(data[, cor(get(PGI1), get(PGI2))])
    # scale PGIs.
    data[, (PGI1) := get(PGI1) / R]
    data[, (PGI2) := get(PGI2) / R]
    cat("Scaling factor =", R, "\n")
    # stack the data
    N = nrow(data)
    data = rbind(data, data)
    data$rep = c(rep(0,N), rep(1,N))
    data[, PGI_MAIN := ifelse(rep==0, get(PGI1), get(PGI2))]
    data[, PGI_IV := ifelse(rep==0, get(PGI2), get(PGI1))]
    # Set up the stacked PGI control variables (parental values).
    if ((!is.null(control.list))) {
        control1 <- control.list[1]
        control2 <- control.list[2]
        data[, PGIcontrol1 := ifelse(rep==0, get(control1), get(control2))]
    }
    # Set up the stacked IVs (random PGI components).
    if ((!is.null(IV.list))) {
        IV1 <- IV.list[1]
        IV2 <- IV.list[2]
        data[, IV1 := ifelse(rep==0, get(IV2), get(IV1))]
    }
    # Define the formulae to estimate.
    formula.long <- paste0(fm[2], "~", fm[3])
    # Add stacked controls, if specified.
    if ((!is.null(control.list))) {
        formula.long <- paste0(formula.long, " + PGIcontrol1 ")
    }
    # Add fam fixed effects, if specified
    if (!is.null(FID)){
        formula.long <- paste0(formula.long, " | ", FID, "^")
    }
    else {
        formula.long <- paste0(formula.long, " | ")
    }
    # Stack the random components as instruments (if specified)
    if (!is.null(IV.list)){
        formula.long <- paste0(formula.long, "rep | PGI_MAIN ~ IV1")
    }
    else {
        formula.long <- paste0(formula.long, "rep | PGI_MAIN ~ PGI_IV")
    }
    cat("Actual number of observations =", N, "\n")
    # Run estimation
    if (!is.null(FID)){
        return(feols(as.formula(formula.long), data=data, se="twoway", cluster=c(FID, IID)))
    } else {
        return(feols(as.formula(formula.long), data=data, se="cluster", cluster=IID))
    }
}


### Example ###

# UKB <- fread("../TEMP/UKB_temp.csv")

# pc <- paste0("pc", 1:20)
# pc <- paste0(pc, collapse="+")
# fm = as.formula(paste0("edu ~ male*factor(yob) + geno + ", pc))

# GORIV(fm, "EA_UKB_PGI", "EA_23_PGI", data=UKB, IID="n_eid", FID="familyID")
# GORIV(fm, "EA_UKB_PGI", "EA_23_PGI", data=UKB, IID="n_eid", FID="familyID", resid=FALSE)

# GORIV(fm, "EA_UKB_PGI", "EA_23_PGI", data=UKB, IID="n_eid", FID="familyID", within=TRUE)
# GORIV(fm, "EA_UKB_PGI", "EA_23_PGI", data=UKB, IID="n_eid", FID="familyID", within=TRUE, resid=FALSE)
GORIV_interaction <- function(fm, interact.var, PGI1, PGI2, IID, data,
        FID = NULL,
        IV.list = NULL,
        control.list = NULL){
    # Put data into data.table format.
    data <- data.table(data)
    # scale PGI (only do if not using random component IVs).
    # compute scaling factor
    R <- sqrt(data[, cor(get(PGI1), get(PGI2))])
    # scale PGIs.
    data[, (PGI1) := get(PGI1) / R]
    data[, (PGI2) := get(PGI2) / R]
    cat("Scaling factor =", R, "\n")
    # stack the data
    N = nrow(data)
    data = rbind(data, data)
    data$rep = c(rep(0,N), rep(1,N))
    data[, PGI_MAIN := ifelse(rep==0, get(PGI1), get(PGI2))]
    data[, PGI_IV := ifelse(rep==0, get(PGI2), get(PGI1))]
    # Set up the stacked PGI control variables (parental values).
    if ((!is.null(control.list))) {
        control1 <- control.list[1]
        control2 <- control.list[2]
        data[, PGIcontrol1 := ifelse(rep==0, get(control1), get(control2))]
    }
    # Set up the stacked IVs (random PGI components).
    if ((!is.null(IV.list))) {
        IV1 <- IV.list[1]
        IV2 <- IV.list[2]
        data[, IV1 := ifelse(rep==0, get(IV2), get(IV1))]
    }
    # Define the formulae to estimate.
    formula.long <- paste0(fm[2], "~", fm[3])
    # Add stacked controls, if specified.
    if ((!is.null(control.list))) {
        formula.long <- paste0(formula.long, " + PGIcontrol1 ")
    }
    # Add fam fixed effects, if specified
    if (!is.null(FID)){
        formula.long <- paste0(formula.long, " | ", FID, "^")
    }
    else {
        formula.long <- paste0(formula.long, " | ")
    }
    # Stack the random components as instruments (if specified)
    if (!is.null(IV.list)){
        formula.long <- paste0(formula.long, "rep | PGI_MAIN *", interact.var,
            " ~ IV1 * ", interact.var)
    }
    else {
        formula.long <- paste0(formula.long, "rep | PGI_MAIN *", interact.var,
            " ~ PGI_IV* ", interact.var)
    }
    cat("Actual number of observations =", N, "\n")
    # Run estimation
    if (!is.null(FID)){
        return(feols(as.formula(formula.long), data=data, se="twoway", cluster=c(FID, IID)))
    } else {
        return(feols(as.formula(formula.long), data=data, se="cluster", cluster=IID))
    }
}




GORIV_mediate <- function(outcome, mediator, fm,
        PGI1, PGI2, IID, data, FID = NULL,
        IV.list = NULL,
        control.list = NULL){
    firststage.reg <- GORIV(formula(paste0(mediator, "~", fm)),
        PGI1, PGI2, IID, data,
        FID, IV.list, control.list)
    print(summary(firststage.reg))
    total.reg <- GORIV(formula(paste0(outcome, "~", fm)),
        PGI1, PGI2, IID, data,
        FID, IV.list, control.list)
    print(summary(total.reg))
    # Put data into data.table format.
    data <- data.table(data)
    # scale PGI (only do if not using random component IVs).
    # compute scaling factor
    R <- sqrt(data[, cor(get(PGI1), get(PGI2))])
    # scale PGIs.
    data[, (PGI1) := get(PGI1) / R]
    data[, (PGI2) := get(PGI2) / R]
    cat("Scaling factor =", R, "\n")
    # stack the data
    N = nrow(data)
    data = rbind(data, data)
    data$rep = c(rep(0,N), rep(1,N))
    data[, PGI_MAIN := ifelse(rep==0, get(PGI1), get(PGI2))]
    data[, PGI_IV := ifelse(rep==0, get(PGI2), get(PGI1))]
    # Set up the stacked PGI control variables (parental values).
    if ((!is.null(control.list))) {
        control1 <- control.list[1]
        control2 <- control.list[2]
        data[, PGIcontrol1 := ifelse(rep==0, get(control1), get(control2))]
    }
    # Set up the stacked IVs (random PGI components).
    if ((!is.null(IV.list))) {
        IV1 <- IV.list[1]
        IV2 <- IV.list[2]
        data[, IV1 := ifelse(rep==0, get(IV2), get(IV1))]
    }
    # Define the formulae to estimate.
    formula.long <- paste0(outcome, "~", fm)
    # Add stacked controls, if specified.
    if ((!is.null(control.list))) {
        formula.long <- paste0(formula.long, " + PGIcontrol1 ")
    }
    # Add fam fixed effects, if specified
    if (!is.null(FID)){
        formula.long <- paste0(formula.long, " | ", FID, "^")
    }
    else {
        formula.long <- paste0(formula.long, " | ")
    }
    # Stack the random components as instruments (if specified)
    if (!is.null(IV.list)){
        formula.long <- paste0(formula.long,
            "rep | PGI_MAIN + ", mediator, " + PGI_MAIN:", mediator,
                " ~ IV1 +", mediator, " + IV1:", mediator)
    }
    else {
        formula.long <- paste0(formula.long, "rep | PGI_MAIN ~ PGI_IV")
    }
    cat("Actual number of observations =", N, "\n")
    # Run estimation
    if (!is.null(FID)){
        secondstage.reg <- feols(as.formula(formula.long), data=data, se="twoway", cluster=c(FID, IID))
    } else {
        secondstage.reg <- feols(as.formula(formula.long), data=data, se="cluster", cluster=IID)
    }
    print(summary(secondstage.reg))
    # Collate the mediation effects.
    firststage.est <- coeftable(firststage.reg)["fit_PGI_MAIN", "Estimate"]
    total.est <- coeftable(total.reg)["fit_PGI_MAIN", "Estimate"]
    direct.est <- coeftable(secondstage.reg)["fit_PGI_MAIN", "Estimate"]
    indirect.est <- coeftable(secondstage.reg)[paste0("fit_", mediator), "Estimate"]
    interaction.est <- coeftable(secondstage.reg)[paste0("fit_PGI_MAIN:", mediator), "Estimate"]
    ade.est <- direct.est + interaction.est * mean(data[[mediator]])
    aie.est <- firststage.est * (
        indirect.est + interaction.est * mean(data[["PGI_MAIN"]]))
    print(c(total.est, firststage.est, direct.est, indirect.est, interaction.est,
        ade.est, aie.est, ade.est, aie.est / total.est))
}
