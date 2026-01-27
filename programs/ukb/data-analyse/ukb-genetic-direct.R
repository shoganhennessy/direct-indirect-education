#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> First-stage analysis of PGI causal design (done with Ed PGIs).
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# Library for equations in plots
library(latex2exp)
# Library to sumamrise models in plots
library(modelsummary)
library(gtsummary)
# The standard, linear, IV estimator package.
library(fixest)
setFixest_notes(FALSE)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.height <- 10
fig.width <- 1.5 * fig.height
presentation.width <- (3 / 2) * fig.width
presentation.height <- presentation.width
# List of 3 default colours.
colour.list <- c(
    "#1f77b4", # Blue
    "#2ca02c", # Green
    "#d62728") # Red
# Define folder path for the general data.
data.folder <- file.path("..", "..", "..", "data", "ukb-restricted", "cleaned")
# Graphics output folder
presentation.folder <- file.path("..", "..", "..", "text", "presentation-files")
figures.folder <- file.path("..", "..", "..", "text", "sections", "figures")
tables.folder <- file.path("..", "..", "..", "text", "sections", "tables")

# Define a function to get the F stat of a specific variable
fstat.get <- function(input.reg, variable.name){
    # Get the variable specific F stats.
    f.test <- car::linearHypothesis(input.reg, test = "F",
        c(paste0(variable.name, "=0")))
    # Extract the statistic.
    f.stat <- f.test["F"] %>%
        unlist() %>%
        nth(2) %>%
        as.numeric() %>%
        round()
    return(f.stat)
}


################################################################################
## Import relevant UKB data. 

# Load the pre-cleaned UKB panel data.
analysis.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv() %>%
    # Get the sibling imputed analysis sample.
    filter(analysis_sample == 1)
# Show the construction
print(analysis.data)
print(names(analysis.data))

# Load the full UKB, for the MLSA IV first-stage.
full.data <- data.folder %>%
    file.path("ukb-cleaned-pheno.csv") %>%
    read_csv() %>%
    filter(!is.na(edyears) & !is.na(soc_median_hourly)) %>%
    transmute(
        eid = eid,
        famid = famid,
        edyears = edyears,
        soc_median_hourly = soc_median_hourly,
        soc_median_annual = soc_median_annual,
        bandwidth = 11,
        cutoff = 1957 + (8 / 12),
        birthyearmonth = birthyear + ((birthmonth - 1) / 12),
        rosla_year = birthyearmonth - cutoff,
        rosla = as.integer(rosla_year >= 0),
        rosla_year = birthyearmonth - cutoff,
        rosla_year_above = rosla_year * rosla,
        rosla_year_below = rosla_year * (1 - rosla),
        dist = abs((birthyearmonth - cutoff) / bandwidth),
        kernel.wt = (1 - dist) * (dist <= 1) / bandwidth) %>%
    filter(kernel.wt > 0)

# Put these MLSA IV columns onto the analysis.data.
analysis.data <- full.data %>%
    select(eid,
        bandwidth,
        cutoff,
        birthyearmonth,
        rosla_year,
        rosla,
        rosla_year,
        rosla_year_above,
        rosla_year_below,
        dist,
        kernel.wt) %>%
    right_join(analysis.data, by = "eid") %>%
    mutate(kernel.wt = ifelse(is.na(kernel.wt), 0, kernel.wt))

# Standardise the controls.
control.formula <- paste0("sex_male + recruitedage + sibling_count + urban +",
    "adhd_pgi + asthma_pgi + bipolar_pgi + bmi_pgi + height_pgi +",
    "schizophrenia_pgi + t2diabetes_pgi")


################################################################################
## 1. OLS mediation results.

# Estimate mediation with OLS.
mediate_ols_edpgi.reg <- function(outcome, mediator, input.data, control.vars,
    indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        input.data <- input.data[indices, ]
    }
    # 1. Total effect (not with ORIV)
    totaleffect.reg <- lm(formula(paste0(outcome,
            "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental + ",
            control.vars)),
        data = input.data)
    #print(summary(totaleffect.reg))
    # 2. Mediation first-stage (controlling for parents).
    firststage.reg <- lm(formula(paste0(mediator,
            "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental + ",
            control.vars)),
        data = input.data)
    # 3. Mediation second-stage (controlling for parents).
    secondstage.reg <- lm(formula(paste0(outcome,
            "~ 1 + edpgi_all_imputed_parental + edpgi_all_imputed_self * ",
            mediator, " + ", control.vars)),
        data = input.data)
    # Extract the total effect.
    total.est <- coeftable(totaleffect.reg)["edpgi_all_imputed_self", "Estimate"]
    # Extract the direct effect.
    direct.effect <- coeftable(secondstage.reg)["edpgi_all_imputed_self", "Estimate"]
    interaction.effect <- coeftable(secondstage.reg)[
        paste0("edpgi_all_imputed_self:", mediator), "Estimate"]
    ade.est <- direct.effect + (interaction.effect * mean(input.data[[mediator]]))
    # Extract the indirect effect.
    firststage.effect <- coeftable(firststage.reg)["edpgi_all_imputed_self", "Estimate"]
    indirect.effect <- coeftable(secondstage.reg)[mediator, "Estimate"]
    aie.est <- firststage.effect * (indirect.effect +
        interaction.effect * mean(input.data[["edpgi_all_imputed_self"]]))
    # Return the estimates.
    output.list <- c(
        firststage.effect,
        total.est,
        ade.est,
        aie.est,
        aie.est / total.est)
    return(output.list)
}


################################################################################
## 2. Two sample IV mediation estimates.

# Estimate mediation with OLS.
mediate_iv_edpgi.reg <- function(outcome, mediator, input.data, control.vars,
    external.data = full.data, poly.count = 2, bandwidth = 11, cutoff = 1957.75,
    indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        input.data <- input.data[indices, ]
        # ALso boot index the full data.
        full.indicies <- sample(
            1:nrow(external.data), nrow(external.data), replace = TRUE)
        external.data <- external.data[full.indicies, ]
    }
    # MLSA IV First-stage, for IV estimated returns to education.
    dist <- abs((external.data$birthyearmonth - cutoff) / bandwidth)
    external.data$kernel.wt <- (1 - dist) * (dist <= 1) / bandwidth
    mlsa_education.reg <- feols(formula(paste0(outcome, " ~ 1 + ",
        "poly(rosla_year_below, poly.count) + ",
        "poly(rosla_year_above, poly.count)",
        " | ", mediator, " ~ rosla")),
        weights = external.data$kernel.wt,
        data = external.data)
    indirect.effect <- coeftable(mlsa_education.reg)[
        paste0("fit_", mediator), "Estimate"]
    # 1. Total effect (not with ORIV)
    totaleffect.reg <- lm(formula(paste0(outcome,
            "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental + ",
            control.vars)),
        weights = kernel.wt,
        data = input.data)
    # 2. Mediation first-stage (controlling for parents).
    firststage.reg <- lm(formula(paste0(mediator,
            "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental + ",
            control.vars)),
        weights = kernel.wt,
        data = input.data)
    # 3. Mediation second-stage (controlling for parents).
    secondstage.reg <- lm(formula(paste0(outcome,
        "~ 1 + edpgi_all_imputed_parental + edpgi_all_imputed_self + ",
        "edpgi_all_imputed_self * ", mediator, " + ", control.vars)),
        weights = kernel.wt,
        data = input.data)
    print(summary(secondstage.reg))
    # Extract the total effect.
    total.est <- coeftable(totaleffect.reg)["edpgi_all_imputed_self", "Estimate"]
    # Extract the direct effect.
    direct.effect <- coeftable(secondstage.reg)["edpgi_all_imputed_self", "Estimate"]
    interaction.effect <- coeftable(secondstage.reg)[
        paste0("edpgi_all_imputed_self:", mediator), "Estimate"]
    ade.est <- direct.effect + (
        interaction.effect * mean(input.data[[mediator]], na.rm = TRUE))
    # Extract the indirect effect.
    firststage.effect <- coeftable(firststage.reg)["edpgi_all_imputed_self", "Estimate"]
    aie.est <- firststage.effect * (indirect.effect + interaction.effect *
        mean(input.data[["edpgi_all_imputed_self"]], na.rm = TRUE))
    #ade.est <- total.est - aie.est
    # Return the estimates.
    output.list <- c(
        firststage.effect,
        total.est,
        ade.est,
        aie.est,
        aie.est / total.est)
    return(output.list)
}

#!TEST:
mediate_iv_edpgi.reg(
    outcome = "log(soc_median_hourly)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula)
    #indices = sample(1:nrow(analysis.data), nrow(analysis.data), replace = TRUE))


################################################################################
## 3. Sibling FE mediation estimates.

# Estimate mediation with sibling FEs.
mediate_fe_edpgi.reg <- function(outcome, mediator, input.data, control.vars,
    indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        input.data <- input.data[indices, ]
    }
    # 1. Total effect (not with ORIV)
    totaleffect.reg <- feols(formula(paste0(outcome,
            "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental + ",
            control.vars, "| birthyear + famid")),
        data = input.data)
    # 2. Mediation first-stage (controlling for parents).
    firststage.reg <- feols(formula(paste0(mediator,
            "~ 1 + edpgi_all_imputed_self + edpgi_all_imputed_parental + ",
            control.vars, "| birthyear + famid")),
        data = input.data)
    # 3. Mediation second-stage (controlling for parents).
    secondstage.reg <- feols(formula(paste0(outcome,
            "~ 1 + edpgi_all_imputed_parental + edpgi_all_imputed_self * ",
            mediator, " + ", control.vars, "| birthyear + famid")),
        data = input.data)
    # Extract the total effect.
    total.est <- coeftable(totaleffect.reg)["edpgi_all_imputed_self", "Estimate"]
    # Extract the direct effect.
    direct.effect <- coeftable(secondstage.reg)["edpgi_all_imputed_self", "Estimate"]
    interaction.effect <- coeftable(secondstage.reg)[
        paste0("edpgi_all_imputed_self:", mediator), "Estimate"]
    ade.est <- direct.effect + (interaction.effect * mean(input.data[[mediator]]))
    # Extract the indirect effect.
    firststage.effect <- coeftable(firststage.reg)["edpgi_all_imputed_self", "Estimate"]
    indirect.effect <- coeftable(secondstage.reg)[mediator, "Estimate"]
    aie.est <- firststage.effect * (indirect.effect +
        interaction.effect * mean(input.data[["edpgi_all_imputed_self"]]))
    # Return the estimates.
    output.list <- c(
        firststage.effect,
        total.est,
        ade.est,
        aie.est,
        aie.est / total.est)
    return(output.list)
}


################################################################################
## Wrapper for boostrapping the mediation estimates.

# Define a function to bootstrap.
mediate.bootstrap <- function(outcome, mediator, input.data, control.vars,
        external.data = full.data,
        type = c("OLS", "IV", "Sibling FE"),
        boot.reps = 10){
    # Define an empty data.frame.
    boot.data <- data.frame(matrix(ncol = 5, nrow = 0))
    names(boot.data) <- c(
        "First-stage", "ATE", "ADE", "AIE", "AIE / ATE")
    j <- 1
    for (i in 1:boot.reps){
        if ((boot.reps >= 100) & ((100 * i / boot.reps) %% 5 == 0)){
            cat(paste0(i, " out of ", boot.reps, ", ", 100 * (i / boot.reps),
                "% done.", "\n"))
        }
        boot.indices <- sample(
            1:nrow(input.data), nrow(input.data), replace = TRUE)
        if (type == "OLS"){
            point.est <- mediate_ols_edpgi.reg(outcome, mediator,
                input.data, control.vars, indices = boot.indices)
            }
        else if (type == "IV"){
            point.est <- mediate_iv_edpgi.reg(outcome, mediator,
                input.data, control.vars, indices = boot.indices)
            }
        else if (type == "Sibling FE"){
            point.est <- mediate_fe_edpgi.reg(outcome, mediator,
                input.data, control.vars, indices = boot.indices)
            }
        else {
            stop(paste0("The type option only takes values of ",
                'c("OLS", "IV", "Sibling FE").'))
        }
        boot.data[i, ] <- point.est
    }
    return(boot.data)
}

# Test it out.
mediate.boot <- mediate.bootstrap(
    outcome = "log(soc_median_hourly)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "IV",
    boot.reps = 10)
print(mediate.boot)

## Define a function to wrap around all the others.
mediate.model <- function(outcome, mediator, input.data, control.vars,
        type = c("OLS", "IV", "Sibling FE"), boot.reps = 10){
    # Calculate the point estimates.
    if (type == "OLS"){
        point.est <- mediate_ols_edpgi.reg(outcome, mediator,
            input.data, control.vars)
    }
    else if (type == "IV"){
        point.est <- mediate_iv_edpgi.reg(outcome, mediator,
            input.data, control.vars)
        }
    else if (type == "Sibling FE"){
        point.est <- mediate_fe_edpgi.reg(outcome, mediator,
            input.data, control.vars)
        }
    else {
        stop(paste0("The type option only takes values of ",
            'c("OLS", "IV", "Sibling FE").'))
    }
    # Calculate the SEs by a non-parametric bootstrap.
    if (!is.null(boot.reps)){
        if (boot.reps < 500){
            print(paste0("Attempting to bootstrap with fewer than 500 reps.",
                "  Are you sure?  This is likely not enough for convergence."))
        }
        point.boot <- mediate.bootstrap(
            outcome, mediator, input.data, control.vars,
            type = type, boot.reps = boot.reps) 
    }
    # Report output
    point.est <- as.matrix(point.est)
    point.se <- as.matrix(c(
        sd(point.boot$"First-stage"),
        sd(point.boot$"ATE"),
        sd(point.boot$"ADE"),
        sd(point.boot$"AIE"),
        sd(point.boot$"AIE / ATE")))
    tratio <- as.matrix(point.est / point.se)
    ptratio <- as.matrix(2 * pt(abs(tratio),
        df = nrow(input.data), lower.tail = FALSE))
    # Preapred object to putput.
    out <- list(
        coefficients = point.est,
        SE = point.se,
        tratio = tratio,
        ptratio = ptratio,
        type = type,
        variables = paste("Ed PGI", mediator, outcome, sep = ", "),
        boot.reps = boot.reps)
    rownames(out$coefficients) <-
        c("First-stage", "ATE", "ADE", "AIE", "Proportion, AIE / ATE")
    rownames(out$SE)      <-rownames(out$coefficients)
    rownames(out$tratio)  <-rownames(out$coefficients)
    rownames(out$ptratio) <-rownames(out$coefficients)
    class(out) <- "mediate.model"
    return(out)
}

# Print applied to the function.
print.mediate.model <- function(x, digits = 4, ...){
    cat("Treatment, Mediator, Outcome: \n")
    cat(x$variables)
    cat("\n")
    est <- cbind(x$coefficients, x$SE)
    colnames(est) <- c("Coefficients", "SE")
    cat(paste0("\n", x$type, " estimates, SEs from ", x$boot.reps, " bootstrap replications."))
    cat("\n\n")
    print.default(format(est, digits = digits), quote = FALSE)
}

# Apply the summary function, to get a presentable output.
summary.mediate.model <- function(object, ...){
    TAP <- cbind(
        Estimate = coef(object),
        SE = object$SE,
        ptratio = object$ptratio)
    colnames(TAP) <- c("Estimate", "SE", "P")
    res <- list(variables = object$variables, coefficients = TAP)  
    class(res) <- "summary.larf"
    return(res)
}

# Presentable summary, via printing.
print.summary.mediate.model <- function(x, digits = 4, ...){
    cat("Treatment, Mediator, Outcome: \n")
    cat(x$variables)
    cat("\n")
    print.default(round(x$coefficients, digits = digits), quote = FALSE)
}


################################################################################
## Estimate mediation on the data, to fill the mediation table.

# Decide how many bootstrap samples to use.
boot.n <- 10^2

# 1. OLS regression
mediate_ols.reg <- mediate.model(
    outcome = "log(soc_median_hourly)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "OLS",
    boot.reps = boot.n)
print(mediate_ols.reg)

# 2. IV regression
mediate_iv.reg <- mediate.model(
    outcome = "log(soc_median_hourly)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "IV",
    boot.reps = boot.n)
print(mediate_iv.reg)

# 3. Sibling FE regression
mediate_fe.reg <- mediate.model(
    outcome = "log(soc_median_hourly)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "Sibling FE",
    boot.reps = boot.n)
print(mediate_fe.reg)


################################################################################
## Estimate mediation on the data, to fill the mediation table.

# Decide how many bootstrap samples to use.
boot.n <- 10^2

# 1. OLS regression
mediate_ols.reg <- mediate.model(
    outcome = "log(soc_median_annual)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "OLS",
    boot.reps = boot.n)
print(mediate_ols.reg)

# 2. IV regression
mediate_iv.reg <- mediate.model(
    outcome = "log(soc_median_annual)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "IV",
    boot.reps = boot.n)
print(mediate_iv.reg)

# 3. Sibling FE regression
mediate_fe.reg <- mediate.model(
    outcome = "log(soc_median_annual)",
    mediator = "edyears",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "Sibling FE",
    boot.reps = boot.n)
print(mediate_fe.reg)

#! Test: for each level education.
mediate_ols.reg <- mediate.model(
    outcome = "log(soc_median_annual)",
    mediator = "edqual_gcses",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "OLS",
    boot.reps = boot.n)
print(mediate_ols.reg)

mediate_ols.reg <- mediate.model(
    outcome = "log(soc_median_annual)",
    mediator = "edqual_alevels",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "Sibling FE",
    boot.reps = boot.n)
print(mediate_ols.reg)
mediate_ols.reg <- mediate.model(
    outcome = "log(soc_median_annual)",
    mediator = "edqual_highered",
    input.data = analysis.data,
    control.vars = control.formula,
    type = "Sibling FE",
    boot.reps = boot.n)
print(mediate_ols.reg)
