#!~/anaconda3/bin/python3
## Senan Hogan-Hennessy, 09 December 2022
## Code to partialling identify an IV system, where
## The exclusion restriction is violated.
## Ths system has single dimension, where the treatment is one-dimensional

# Various packages for use.
import numpy as np
npSeed = np.random.default_rng(47)
import pandas as pd
import os
# data viewer for exploration
from gtabview import view
# Plotting
import matplotlib.pyplot as plt
figWidth = 5
figHeight = figWidth / 1.5
# Linear models
from scipy import stats
import statsmodels.api as sm
from stargazer.stargazer import Stargazer
# Instrument linear estimation
from linearmodels.iv import IV2SLS
from linearmodels.iv import compare


################################################################################
## Set up the 2SLS system, following the Conley+ (2012) specification, for
## violations to the exclusion restriction.

## Set up the parameters to begin with
# Length of simulation
nLength = 10 ** 4
# True parameters for the following specification,
# where the first stages are identified, but the second stage partially 
# First stage X:  X = φ_x + π Z + ε_x
# Second stage Y: Y = φ_y + β X + γ Z + ε_y
# Reduced form Y: Y = φ   + α Z + ε (where α = βπ + γ)
φ_x = 1
π = 1
φ_y = 1
β = 3
γ = 2
α = (β * π) + γ

# Generate the instrument, Z, and error terms.
Z = npSeed.binomial(1, 1 / 2, nLength)
ε_x = npSeed.normal(0, 1, nLength)
ε_y = npSeed.normal(0, 1, nLength)
# Generate data for X, with selection on Y (which confounds the OLS est for β).
X = [2 * φ_x + π * Z[i] + ε_x[i] if ε_y[i] >= 0
    else φ_x + π * Z[i] + ε_x[i]
    for i in range(0, nLength)]
# Generate the list of data for Y
Y = [φ_y + β * X[i] + γ * Z[i] + ε_y[i] for i in range(0, nLength)]

# Put the data, of observed data (X, Y, Z), to a dataframe
simData = pd.DataFrame([[X[i], Y[i], Z[i]] for i in range(0, nLength)],
    columns=["X", "Y", "Z"])


################################################################################
## Estimation of the system

# First-stage OLS estimate of π (this model is identified).
firststageEstimate = IV2SLS.from_formula("X ~ 1 + Z", data=simData)
firststageEstimate = firststageEstimate.fit()
firststageEstimate = firststageEstimate._params["Z"]

# Naive OLS estimate of β
olsEstimate = IV2SLS.from_formula("Y ~ 1 + X", data=simData)
olsEstimate = olsEstimate.fit()
olsEstimate = olsEstimate._params["X"]

# SHow how badly OLS performs, including both X + Z.
fullEstimate = IV2SLS.from_formula("Y ~ 1 + X + Z", data=simData).fit()
print(β, γ)
print(fullEstimate._params["X"], fullEstimate._params["Z"])

# Naive IV estimate of β, which corrsponds to the restricted version γ = 0
ivEstimate = IV2SLS.from_formula("Y ~ 1 + [X ~ Z]", data=simData)
ivEstimate = ivEstimate.fit()
ivEstimateSE = ivEstimate.std_errors["X"]
ivEstimate = ivEstimate._params["X"]

# Calculate the partially identified point estimates for (γ, β), β = (α - γ) / π
πHat = IV2SLS.from_formula("X ~ 1 + Z", data=simData).fit()._params["Z"]
αHat = IV2SLS.from_formula("Y ~ 1 + Z", data=simData).fit()._params["Z"]
γHatList = np.linspace(0, αHat, nLength, endpoint=True).tolist()
βHatList = [(αHat - γHatList[i]) / πHat for i in range(0, nLength)]
βUpperList = [1.96 * ivEstimateSE + (
    (αHat - γHatList[i]) / πHat) for i in range(0, nLength)]
βLowerList = [-1.96 * ivEstimateSE + (
    (αHat - γHatList[i]) / πHat) for i in range(0, nLength)]
print("α value and estimate", [α, αHat])
print("π value and estimate", [π, πHat])
print("β value and estimate", [β, olsEstimate])
print("γ value and estimate", [γ, 0])


################################################################################
## Plot the sensitivity analysis.

# Show the partially identified values of the estimates for (γ, β).
plt.figure(figsize=(figWidth, figHeight))
# Plot the Real value
plt.axhline(y = β, color = "r", linestyle = "--", alpha=0.5)
plt.text(β, β + 0.1, r"$\beta$", color = "r")
plt.axvline(x = γ, color = "r", linestyle = "--", alpha=0.5)
plt.text(γ + 0.1, ivEstimate - 0.1, r"$\gamma$", color = "r")
# Plot the naive OLS estimate
plt.axhline(y = olsEstimate, color = "orange", linestyle = "--", alpha=0.5)
plt.text(β - 0.2, olsEstimate - 0.2, r"$\widehat{\beta}_{OLS}$", color = "orange")
# Plot the naive IV estimate
plt.axhline(y = ivEstimate, color = "orange", linestyle = "--", alpha=0.5)
plt.text(β - 0.2, ivEstimate - 0.2, r"$\widehat{\beta}_{IV}$", color = "orange")
# Plot the partially identified estimates (and 95% confidence interval)
plt.plot(γHatList, βUpperList, linestyle="--", alpha=0.25, color="grey")
plt.plot(γHatList, βLowerList, linestyle="--", alpha=0.25, color="grey")
plt.fill_between(γHatList, βLowerList, βUpperList, alpha=0.25, color="grey")
plt.plot(γHatList, βHatList)
plt.xlabel(r"$\widehat{\gamma}$", fontsize=15)
plt.ylabel(r"$\widehat{\beta}$", fontsize=15, rotation=0, loc="top")
#plt.show()
plt.savefig("../../paper/figures/single-dimension.png",
    dpi=300, transparent=False, bbox_inches="tight")

#! Conclusions:
# This simulation shows how to use the Conely+ (2012) method 
# for sensistivity analysis.

# Note that if β, γ >= 0 is certain, then 0 <= γ <= α provide natural bounds.


#TODO: Next steps
#TODO: 2. Look into further theory on how to advance this model, and inform 
#TODO:    which applications are most appropriate.
#TODO:    (e.g the JMP on identifying IV compliers)
#TODO: 2. Look at the Roth (2022) sensitivity analysis, and see whether a similar
#TODO:    Approach is applicable here.
#TODO: 3. Write code to point identify (perhaps by Bayesian/MLE) the IV model.
#TODO:    where the priors are context dependent
#TODO:    -> consider which priors achieve identification, and thus which contexts
#TODO:       this model works for.
