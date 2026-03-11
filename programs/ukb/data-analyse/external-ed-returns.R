#!/usr/bin/R
## Senan Hogan-Hennessy, 23 Feb 2026.
## Meta-analysis of ed returns.
set.seed(47)
print(Sys.time())

# Functions for data manipulation and visualisation
library(tidyverse)
# Functions for tables into TeX
library(xtable)
# Library for equations in plots
library(latex2exp)
# Library to sumamrise models in plots
library(gtsummary)
library(modelsummary)
# The standard, linear, IV estimator package.
library(fixest)
library(ivreg)
library(rdrobust)
# Define number of digits in tables and graphs
digits.no <- 3
# Size for figures
fig.height <- 8.5
fig.width <- 1.25 * fig.height
presentation.width <- (3 / 2) * fig.width
presentation.height <- presentation.width
# List of 3 default colours.
colour.list <- c(
    "#1f77b4", # Blue
    "#2ca02c", # Green
    "#d62728", # Red
    "Sea Green")      # Dark green (for similar to green variables).
# Define folder path for the general data.
data.folder <- file.path("..", "..", "..", "data", "ed-returns")
# Graphics output folder
figures.folder <- file.path("..", "..", "..", "text", "sections", "figures")
tables.folder <- file.path("..", "..", "..", "text", "sections", "tables")


################################################################################
## Import ed returns data.

# Load the spreadsheet, as provided.
edreturns.data <- data.folder %>%
    file.path("causal-returns-database.xlsx") %>%
    readxl::read_excel(sheet = 2)

# Show the construction
print(edreturns.data)
print(names(edreturns.data))


################################################################################
## Plot the UK distribution.

# Plot predicted edyears, before, after.
uk.plot <- edreturns.data %>%
    filter(Country == "UK") %>%
    ggplot(aes(x = IV)) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
    # Density plot.
    geom_density(colour = "black",, fill = colour.list[2], alpha = 0.75) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Percent Effect of an Extra Education Year on Income",
        limits = 100 * c(-0.04, 0.275),
        breaks = 100 * seq(0, 0.5, by = 0.05)) +
    scale_y_continuous(expand = c(0, 0),
        limits = (1 / 100) * c(-0.05, 15.5),
        breaks = (1 / 100) * seq(0, 15, by = 2.5),
        name = "") +
    ggtitle(TeX(r"(Density of British IV Estimates)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edreturns-uk.png"),
    plot = uk.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the UK distribution.

# Plot education returns, among the entire world.
world.plot <- edreturns.data %>%
    ggplot(aes(x = IV)) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
    # Density plot.
    geom_density(color = "black", fill = colour.list[4], alpha = 0.75) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Percent Effect of an Extra Education Year on Income",
        limits = 100 * c(-0.04, 0.275),
        breaks = 100 * seq(0, 0.5, by = 0.05)) +
    scale_y_continuous(expand = c(0, 0),
        limits = (1 / 100) * c(-0.01, 8.5),
        breaks = (1 / 100) * seq(0, 10, by = 1),
        name = "") +
    ggtitle(TeX(r"(Density of Worldwide IV Estimates)")) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1, 0), "mm"))
# Save this plot
ggsave(file.path(figures.folder, "edreturns-world.png"),
    plot = world.plot,
    units = "cm", width = fig.width, height = fig.height)
