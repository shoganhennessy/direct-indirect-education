#!/usr/bin/R
## Senan Hogan-Hennessy, 24 March 2025
## UKB data -> Genetic effects, Ed PGI random component -> Ed Years + Income.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# The standard, linear, IV estimator package.
library(ivreg)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.width <- 7.5
fig.height <- fig.width
presentation.width <- (5 / 3) * fig.width
presentation.height <- (2 / 3) * presentation.width
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


################################################################################
## Convenience functions.

# Define a function to automate Bin scatter.
binscatter.plot <- function(data, x, y, colour.name, option = "none"){
    library(binsreg)
    # Run the binscatter regression.
    binscatter.data <- binsreg(data[[y]], data[[x]], randcut = 1,
            line = c(1, 1), ci = c(1, 1), cb = c(1, 1), polyreg = 3)  %>%
        magrittr::use_series("data.plot") %>%
        magrittr::use_series("Group Full Sample")
    # Plot the binscatter
    binscatter.ggplot <- ggplot() +
        # Add the bin means
        geom_point(data = binscatter.data$data.dots, aes(x = x, y = fit),
            colour = colour.name, size = 2) +
        # Add the bin error bars
        #geom_errorbar(data = binscatter.data$data.ci,
        #    aes(x = x, ymin = ci.l, ymax = ci.r),
        #    colour = colour.name, size = 0.5, width = 0.02, linetype = "solid") +
        geom_ribbon(data = binscatter.data$data.cb,
            aes(x = x, ymin = cb.l, ymax = cb.r),
            fill = colour.name, alpha = 0.2) +
        # Add the line between bin means
        geom_smooth(data = binscatter.data$data.dots, aes(x = x, y = fit),
            se = FALSE, colour = colour.name, size = 0.5, linetype = "solid")
        if (option == "half-line"){
            binscatter.ggplot <- binscatter.ggplot +
                geom_smooth(data = binscatter.data$data.dots,
                aes(x = x, y = (6 + fit / 2)), se = FALSE,
                colour = "orange", size = 1, linetype = "solid")
        }
    # Return the ggplot of this.
    return(binscatter.ggplot)
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


################################################################################
## Plot comparison of causal med (1) OLS (2) OLS + controls (3) DML.

# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)

# Function to extract the coefficients from the mediate objects
effects.extract <- function(mediate.reg, model.name){
    # Compile the mediation regression results.
    reg.summary <- summary(mediate.reg)
    # Get the total effect estimates.
    total.est <- reg.summary$tau.coef
    total.ci.upper <- as.numeric(reg.summary$tau.ci[1])
    total.ci.lower <- as.numeric(reg.summary$tau.ci[2])
    # Get the direct effect estimates.
    direct.est <- reg.summary$z.avg
    direct.ci.upper <- as.numeric(reg.summary$z.avg.ci[1])
    direct.ci.lower <- as.numeric(reg.summary$z.avg.ci[2])
    # Get the indirect effect estimates.
    indirect.est <- reg.summary$d.avg
    indirect.ci.upper <- as.numeric(reg.summary$d.avg.ci[1])
    indirect.ci.lower <- as.numeric(reg.summary$d.avg.ci[2])
    # Get the percent mediated estimates.
    permediated.est <- reg.summary$n.avg
    permediated.ci.upper <- as.numeric(reg.summary$n.avg.ci[1])
    permediated.ci.lower <- as.numeric(reg.summary$n.avg.ci[2])
    # Put it all into a dataframe
    data.return <- data.frame(
        effect = c("Total", "Direct", "Indirect", "Percent Mediated"),
        pointest = c(total.est, direct.est, indirect.est, permediated.est),
        upperest = c(total.ci.upper, direct.ci.upper, indirect.ci.upper, permediated.ci.upper),
        lowerest = c(total.ci.lower, direct.ci.lower, indirect.ci.lower, permediated.ci.lower))
    # Label it with the model name
    data.return <- data.return %>% mutate(model = model.name)
    return(data.return)
}

# Put relevant data into matrix form
Z <- analysis.data %>% pull(edpgi_self)
D <- analysis.data %>% pull(indiv_edyears)
Y <- analysis.data %>% pull(indiv_earnings_real) %>% log(.)
X <- ukb.data %>%
    dplyr::select(c(
        # Regular demographics
        survey_year, indiv_agey, parent_edyears, gender_female,
        # Childhood SES
        family_poor, family_move, family_finhelp,
        father_missing, father_unemp, father_manualjob,
        parents_smoke, child_headinjury)) %>%
    as.matrix()
# Consider D = higher education (binary needed)
#Z <- as.integer(Z >= 0)
D <- as.integer(D >= 14)

# Define number of samples to bootstrap over.
boot.samples <- 10^4

# Define the simple OLS model
mediation_ols.reg <- 
    mediate(
        lm(D ~ (1 + Z)),
        lm(Y ~ (1 + Z * D)),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_ols.reg))

# Define the OLS model with linear controls
mediation_controls.reg <- 
    mediate(
        lm(D ~ (1 + Z) + X),
        lm(Y ~ (1 + Z * D) + X),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_controls.reg))

# Define the OLS model with non-linear controls
mediation_nonparametric.reg <-
    mediate(
        lm(D ~ (1 + Z) * X),
        lm(Y ~ (1 + Z * D) * X),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_nonparametric.reg))

# Get all the results together.
effects.data <- rbind(
    effects.extract(mediation_ols.reg, "OLS"),
    effects.extract(mediation_controls.reg, "+ Controls"),
    effects.extract(mediation_nonparametric.reg, "+ Non-Linear"))

# Plot in a bar chart.
mediation.plot <- effects.data %>%
    filter(effect != "Percent Mediated") %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("OLS", "+ Controls", "+ Non-Linear")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        #limits = c(0, 0.12),
        name = "",
        breaks = seq(0, 0.5, by = 0.025)) +
    ggtitle("Estimate, Percent on Earnings") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-plot.png"),
    plot = mediation.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Plot in a bar chart, just the total effect
mediation_total.plot <- effects.data %>%
    #filter(effect == "Total") %>%
    mutate(
        pointest = ifelse(effect == "Total", pointest, 0),
        upperest = ifelse(effect == "Total", upperest, 0),
        lowerest = ifelse(effect == "Total", lowerest, 0)) %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("OLS", "+ Controls", "+ Non-Linear")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        #limits = c(0, 0.12),
        name = "",
        breaks = seq(0, 0.5, by = 0.025)) +
    ggtitle("Estimate, Percent on Earnings") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-total-plot.png"),
    plot = mediation_total.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Plot in a bar chart, just the percent mediatied.
mediation_percent.plot <- effects.data %>%
    filter(effect == "Percent Mediated") %>%
    ggplot(aes(x = model)) +
    geom_bar(aes(y = pointest), fill = "orange",
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_x_discrete(
        name = "",
        limits = c("OLS", "+ Controls", "+ Non-Linear")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(0, 1),
        name = "",
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Estimate, Percent Mediated Through Higher Education") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-percent-plot.png"),
    plot = mediation_percent.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)


################################################################################
## Selection model estimates of causal mediation

# Point estimate with the Heckman parametric selection model.
library(sampleSelection)
library(boot)
boot.samples <- 10^4
# Put relevant data into matrix form
Z <- ukb.data %>% pull(edpgi_self)
D <- ukb.data %>% pull(indiv_edyears)
Y <- ukb.data %>% transmute(Y = log(indiv_earnings_real)) %>% pull(Y)
X <- ukb.data %>%
    dplyr::select(c(
        survey_year, indiv_agey, parent_edyears, gender_female)) %>%
    as.matrix()
# Consider D = higher education (binary needed)
#Z <- as.integer(Z >= 0)
D <- as.integer(D >= 14)

# Define a function that outputs a parametric selection model's estimates of
# (1) total effect (2) First-stage (3) direct (4) indirect
selection_mediation.reg <- function(data, type = "parametric"){
    # Get the relevant columns from the provided data.
    Z <- data[["Z"]]
    D <- data[["D"]]
    Y <- data[["Y"]]
    X <- as.matrix(data[4:dim(data)[2]])
    # Unified first-stage
    firststage.reg <- probit(D ~ (1 + Z) * X)
    # Unified total effect
    totaleffect.reg <- lm(Y ~ (1 + Z) * X)
    ## Get the prediction models by type
    if (type == "parametric"){
        # Get the second-stage by using a parametric control function.
        # Unobserved part estimated with a Heckman (switching) selection approach.
        selection.reg <- selection(
            # First-stage
            as.formula(D ~ (1 + Z) * X),
            # Second-stage
            list(as.formula(Y ~ (1 + Z)), as.formula(Y ~ (1 + Z) * X)),
                method = "2step")
        # Generate predictions
        EY_Zz_D0 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 0, X = X))[, 1]
        EY_Zz_D1 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 1, X = X))[, 2]
        EY_Z1_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z + 0.5, D = D, X = X))
        EY_Z1_Dd <- ifelse(D == 1, EY_Z1_Dd[, 2], EY_Z1_Dd[, 1])
        EY_Z0_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z - 0.5, D = D, X = X))
        EY_Z0_Dd <- ifelse(D == 1, EY_Z0_Dd[, 2], EY_Z0_Dd[, 1])
    }
    else if (type == "semi-parametric"){
        controlfun <- (D - predict(firststage.reg, type = "response"))
        selection.reg <- lm(Y ~ (1 + Z * D) + controlfun * D)
        # Generate predictions
        EY_Zz_D0 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 0, X = X,
                controlfun = controlfun), rankdeficient = "simple")
        EY_Zz_D1 <- predict(
            selection.reg, newdata = data.frame(Z = Z, D = 1, X = X,
                controlfun = controlfun), rankdeficient = "simple")
        EY_Z1_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z + 0.5, D = D, X = X,
                controlfun = controlfun), rankdeficient = "simple")
        EY_Z0_Dd <- predict(
            selection.reg, newdata = data.frame(Z = Z - 0.5, D = D, X = X,
                controlfun = controlfun), rankdeficient = "simple")
    }
    # Get the total effect by prediction
    totaleffect.est <- mean(
        (predict(totaleffect.reg, newdata = data.frame(Z = Z + 0.5, X = X)) -
            predict(totaleffect.reg, newdata = data.frame(Z = Z - 0.5, X = X))))
    # Get the first-stage by prediction
    complier.score <- (predict(firststage.reg, 
        newdata = data.frame(Z = Z + 0.5, X = X), type = "response") -
            predict(firststage.reg,
                newdata = data.frame(Z = Z - 0.5, X = X), type = "response"))
    firststage.est <- mean(complier.score)
    # Direct Effect, E[ Y_i(1, D) - Y_i(0, D)].
    direct.gains <- EY_Z1_Dd - EY_Z0_Dd
    direct.est <- mean(direct.gains)
    # Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
    treatment.gains <- EY_Zz_D1 - EY_Zz_D0
    indirect.est <- mean(complier.score * treatment.gains)
    treatment.est <- weighted.mean(treatment.gains, complier.score)
    # Percent mediated
    # note: probably does not follow a normal at boundary -> bootstrap inapprop.
    mediated.est <- indirect.est / totaleffect.est
    # Return the relevant statistics
    return(c(totaleffect.est, firststage.est,
        direct.est, indirect.est, mediated.est, treatment.est))
}

# Show how the function works
selection.example <- selection_mediation.reg(
    data.frame(Z = Z, D = D, Y = Y, X = X), type = "parametric")
print(selection.example)
# Or the semi-para est
selection.example <- selection_mediation.reg(
    data.frame(Z = Z, D = D, Y = Y, X = X), type = "semi-parametric")
print(selection.example)

# Define the bootstrap function.
selection_mediation.boot <- function(data, indicies){
    estimates.list <- selection_mediation.reg(
        data[indicies, ], type = "parametric")
    return(estimates.list)
}

# Define the bootstrap function.
selection_semiparametric.boot <- function(data, indicies){
    estimates.list <- selection_mediation.reg(
        data[indicies, ], type = "semi-parametric")
    return(estimates.list)
}

# Bootstrap the selection mediated estimates.
selection_mediation.bootest <-
    boot(data = data.frame(Z = Z, D = D, Y = Y, X = X),
        statistic = selection_mediation.boot, R = boot.samples)
selection_semiparametric.bootest <-
    boot(data = data.frame(Z = Z, D = D, Y = Y, X = X),
        statistic = selection_semiparametric.boot, R = boot.samples)

# Function to extract the coefficients from the bootstrapped estimates
bootstrap.extract <- function(model.bootest, model.name){
    # Extract the columns of the boot object by name
    total.boot       <- model.bootest$t[, 1]
    firststage.boot  <- model.bootest$t[, 2]
    direct.boot      <- model.bootest$t[, 3]
    indirect.boot    <- model.bootest$t[, 4]
    permediated.boot <- model.bootest$t[, 5]
    # Get the total effect estimates.
    total.est      <- mean(total.boot)
    total.se       <- sd(total.boot)
    total.ci.upper <- as.numeric(quantile(total.boot, probs = 0.975))
    total.ci.lower <- as.numeric(quantile(total.boot, probs = 0.025))
    # Get the direct effect estimates.
    direct.est      <- mean(direct.boot)
    direct.se       <- sd(direct.boot)
    direct.ci.upper <- as.numeric(quantile(direct.boot, probs = 0.975))
    direct.ci.lower <- as.numeric(quantile(direct.boot, probs = 0.025))
    # Get the indirect effect estimates.
    indirect.est      <- mean(indirect.boot)
    indirect.se       <- sd(indirect.boot)
    indirect.ci.upper <- as.numeric(quantile(indirect.boot, probs = 0.975))
    indirect.ci.lower <- as.numeric(quantile(indirect.boot, probs = 0.025))
    # Get the percent mediated estimates.
    permediated.est      <- mean(permediated.boot)
    permediated.se       <- sd(permediated.boot)
    permediated.ci.upper <- as.numeric(quantile(permediated.boot, probs = 0.975))
    permediated.ci.lower <- as.numeric(quantile(permediated.boot, probs = 0.025))
    # Put it all into a dataframe
    data.return <- data.frame(
        effect = c("Total", "Direct", "Indirect", "Percent Mediated"),
        pointest = c(total.est, direct.est, indirect.est, permediated.est),
        upperest = c(total.ci.upper, direct.ci.upper, indirect.ci.upper, permediated.ci.upper),
        lowerest = c(total.ci.lower, direct.ci.lower, indirect.ci.lower, permediated.ci.lower))
    # Label it with the model name
    data.return <- data.return %>% mutate(model = model.name)
    return(data.return)
}

# Get the regular mediation estimates of the binary model.
mediation_binary.reg <-
    mediate(
        lm(D ~ (1 + Z) * X),
        lm(Y ~ (1 + Z * D) * X),
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = boot.samples)
print(summary(mediation_binary.reg))

# Append the estimates to those extracted earlier
selection_est.data <- rbind(
    effects.extract(mediation_binary.reg, "OLS + Controls"),
    bootstrap.extract(selection_mediation.bootest,
        "Parametric Heckman\nSelection Model"),
    bootstrap.extract(selection_semiparametric.bootest,
        "Semi-Parametric\nSelection Model"))

# Plot in a bar chart.
mediation_selection.plot <- selection_est.data %>%
    filter(effect != "Percent Mediated") %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("OLS + Controls",
            "Parametric Heckman\nSelection Model",
            "Semi-Parametric\nSelection Model")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        #limits = c(-0.03, 0.175),
        name = "",
        breaks = seq(-0.5, 0.5, by = 0.025)) +
    ggtitle("Estimate, Percent on Earnings") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-selection-plot.png"),
    plot = mediation_selection.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Plot in a bar chart, just the percent mediatied.
mediation_selection_percent.plot <- selection_est.data %>%
    #mutate(
    #    pointest = ifelse(pointest > 1, 1, pointest),
    #    upperest = ifelse(upperest > 1, 1, upperest),
    #    lowerest = ifelse(lowerest > 1, 1, lowerest)) %>%
    filter(effect == "Percent Mediated") %>%
    ggplot(aes(x = model)) +
    geom_bar(aes(y = pointest), fill = "orange",
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_hline(yintercept = 1, alpha = 0.5, linetype = "dotted") +
    scale_x_discrete(
        name = "",
        limits = c("OLS + Controls",
            "Parametric Heckman\nSelection Model",
            "Semi-Parametric\nSelection Model")) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0),
        #limits = c(0, 1),
        name = "",
        breaks = seq(-2, 2, by = 0.1)) +
    ggtitle("Estimate, Percent Mediated Through Higher Education") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = c(0.775, 0.925),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-selection-percent-plot.png"),
    plot = mediation_selection_percent.plot,
    units = "cm",
    width = presentation.width * 1.125, height = presentation.height * 1.125)

# Question: What are the implied returns to education estimates?
