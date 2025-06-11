
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

GORIV <- function(fm, PGI1, PGI2, data, IID,
        FID = NULL, resid = FALSE, within = FALSE){
    # Put data into data.table format.
    data <- data.table(data)
    data <- na.omit(data[, .SD, .SDcols=unique(c(all.vars(fm), IID, FID, PGI1, PGI2))])
    if (within) {
        data <- data[get(FID) %in% data[, .N, by=get(FID)][N>=2][[1]] ]
    }

    # Residualise / Set up model
    if (resid==TRUE){
        data[, Y := lm(fm, data)$residuals]
        data[, Y := scale(Y)]
        if (within==TRUE){
            FM = paste0("Y ~ 1 | ", FID, "^rep | PGI_MAIN ~ PGI_IV")
        } else {
            FM = "Y ~ 1 | rep | PGI_MAIN ~ PGI_IV"
        }

    }
    else if (resid==FALSE){
        fm = as.character(fm)
        if (within==TRUE){
            FM = paste0(fm[2], "~", fm[3], "| ", FID, "^rep | PGI_MAIN ~ PGI_IV")
        } else {
            FM = paste0(fm[2], "~", fm[3], "| rep | PGI_MAIN ~ PGI_IV")
        }
    }

    # standardise observed PGI
    data[, c(PGI1, PGI2) := lapply(.SD, scale), .SDcols=c(PGI1, PGI2)]
    
    # compute scaling factor
    if (within){
        data[, fam_n := .N, by=FID]
        data <- data[fam_n >= 2]
        data[, PGI1_dm := get(PGI1) - mean(get(PGI1)), by=FID]
        data[, PGI2_dm := get(PGI2) - mean(get(PGI2)), by=FID]
        R <- sqrt(data[, cor(PGI1_dm, PGI2_dm)])
    } else {
        R <- sqrt(data[, cor(get(PGI1), get(PGI2))])
    }

    # scale PGI
    data[, (PGI1) := get(PGI1) / R]
    data[, (PGI2) := get(PGI2) / R]
    cat("Scaling factor =", R, "\n")
    
    # stack the data
    N = nrow(data)
    data = rbind(data, data)
    data$rep = c(rep(0,N), rep(1,N))
    
    data[, PGI_MAIN := ifelse(rep==0, get(PGI1), get(PGI2))]
    data[, PGI_IV := ifelse(rep==0, get(PGI2), get(PGI1))]
    
    # Run estimation
    if (within){
        return(feols(as.formula(FM), data=data, se="twoway", cluster=c(FID, IID)))
    } else if (!is.null(FID)){
        return(feols(as.formula(FM), data=data, se="twoway", cluster=c(FID, IID)))
    } else {
        return(feols(as.formula(FM), data=data, se="cluster", cluster=IID))
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
