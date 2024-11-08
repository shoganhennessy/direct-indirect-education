#!~/anaconda3/bin/python3
## Senan Hogan-Hennessy, 09 December 2022
## Code to partialling identify an IV system, where
## The exclusion restriction is violated.
## Ths system has two dimensions,
## 2 instruments which affect X in opposite directions, for upper + lower bounds

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
nLength = 10 ** 3
# True parameters for the following specification,
# where the first stages are identified, but the second stage partially 
# First stage X:  X = φ_x + π Z + ε_x
# Second stage Y: Y = φ_y + β X + γ Z + ε_y
# Reduced form Y: Y = φ   + α Z + ε (where α = βπ + γ)
β = 4
φ_x = 1
φ_y = 1
π_1 = 2
π_2 = 2
γ_1 = β * (1 / 2)
γ_2 = β * (- 5 / 3)
α_1 = (β * π_1) + γ_1
α_2 = (β * π_2) + γ_2

# Generate the instrument, Z, and error terms.
Z_1 = npSeed.binomial(1, 0.5, nLength)
Z_2 = npSeed.binomial(1, 0.5, nLength)
ε_x = npSeed.normal(0, 1, nLength)
ε_y = npSeed.normal(0, 1, nLength)
# Generate data for X, with selection on Y (which confounds the OLS est for β).
X = [3 * φ_x + π_1 * Z_1[i] + π_2 * Z_2[i] + ε_x[i] if ε_y[i] >= 0
    else φ_x + π_1 * Z_1[i] + π_2 * Z_2[i] + ε_x[i]
    for i in range(0, nLength)]
# Generate the list of data for Y
Y = [φ_y + β * X[i] + γ_1 * Z_1[i] + γ_2 * Z_2[i] + ε_y[i]
    for i in range(0, nLength)]

# Put the data, of observed data (X, Y, Z_1, Z_2), to a dataframe
simData = pd.DataFrame(
    [[X[i], Y[i], Z_1[i], Z_2[i]] for i in range(0, nLength)],
    columns=["X", "Y", "Z_1", "Z_2"])


################################################################################
## Estimation of the system

# Naive IV estimate of β from only Z_1, which corrsponds to γ_1 = 0
ivEstimate_1 = IV2SLS.from_formula("Y ~ 1 + [X ~ Z_1]", data=simData)
ivEstimate_1 = ivEstimate_1.fit()
ivEstimate_1SE = ivEstimate_1.std_errors["X"]
ivEstimate_1 = ivEstimate_1._params["X"]

# Naive IV estimate of β from only Z_2, which corrsponds to γ_2 = 0
ivEstimate_2 = IV2SLS.from_formula("Y ~ 1 + [X ~ Z_2]", data=simData)
ivEstimate_2 = ivEstimate_2.fit()
ivEstimate_2SE = ivEstimate_2.std_errors["X"]
ivEstimate_2 = ivEstimate_2._params["X"]

# Naive IV estimate of β from only Z_1 + Z_2, which corrsponds to γ_1 = γ_2 = 0
ivEstimate = IV2SLS.from_formula("Y ~ 1 + [X ~ Z_1 + Z_2]", data=simData)
ivEstimate = ivEstimate.fit()
ivEstimateSE = ivEstimate.std_errors["X"]
ivEstimate = ivEstimate._params["X"]

# Calculate the identified point estimates, for each instrument.
π_1Hat = IV2SLS.from_formula("X ~ 1 + Z_1 + Z_2", data=simData).fit()._params["Z_1"]
π_2Hat = IV2SLS.from_formula("X ~ 1 + Z_1 + Z_2", data=simData).fit()._params["Z_2"]
α_1Hat = IV2SLS.from_formula("Y ~ 1 + Z_1 + Z_2", data=simData).fit()._params["Z_1"]
α_2Hat = IV2SLS.from_formula("Y ~ 1 + Z_1 + Z_2", data=simData).fit()._params["Z_2"]
# Calculate the partially identified point estimates
# (γ_1, β_1), β_1 = (γ_1 - α_1) / π_1
γ_1HatList = np.linspace(-α_1Hat, 0, nLength, endpoint=True).tolist()
β_1HatList = [(α_1Hat - γ_1HatList[i]) / π_1Hat for i in range(0, nLength)]
# Calculate the partially identified point estimates
# (γ_2, β_2), β_2 = (γ_2 - α_2) / π_2
γ_2HatList = np.linspace(0, α_2Hat, nLength, endpoint=True).tolist()
β_2HatList = [(α_2Hat - γ_2HatList[i]) / π_2Hat for i in range(0, nLength)]
print("α_1 value and estimate", [α_1, α_1Hat])
print("α_2 value and estimate", [α_2, α_2Hat])
print("π_1 value and estimate", [π_1, π_1Hat])
print("π_2 value and estimate", [π_2, π_2Hat])
print("β value", β)
print("γ_1 =", γ_1,
    "so γ_1 <= 0 gives natural β lower bound is β_IV", ivEstimate_1)
print("γ_2 =", γ_2,
    "so γ_2 >= 0 gives natural β upper bound is β_IV", ivEstimate_2)


print("β estimate by Z_1 instrument " + str(ivEstimate_1))
print("β estimate by Z_1 + Z_2 instruments " + str(ivEstimate))
print("β estimate by Z_2 instrument " + str(ivEstimate_2))
print("Mid point of each bound " + str((ivEstimate_1 + ivEstimate_2) / 2))


################################################################################
## Plot the sensitivity analysis.

# Show the partially identified values of the estimates for (γ_1, β_1).
plt.figure(figsize=(figWidth, figHeight))
plt.xlabel(r"$\widehat{\gamma}$", fontsize=15)
plt.ylabel(r"$\widehat{\beta}$", fontsize=15, rotation=0, loc="top")
# Plot the Real value
plt.axhline(y=β, color="red", linestyle="--", alpha=0.5)
plt.axvline(x=γ_1, color="C0", linestyle="--", alpha=0.5)
plt.axvline(x=γ_2, color="C1", linestyle="--", alpha=0.5)
# Plot the multiple instrument estimate
plt.axhline(y=ivEstimate, color="C3", linestyle="--", alpha=0.5)
# Plot the partially identified estimates
plt.plot(γ_1HatList, β_1HatList)
# Plot the partially identified estimates
plt.plot(γ_2HatList, β_2HatList)
# Plot the bounded region
plt.annotate("", xy=(0, ivEstimate_1), xytext=(0, ivEstimate_2),
    #draws an arrow from one set of coordinates to the other
    arrowprops=dict(arrowstyle="<->"),
    annotation_clip=False)

## Save the figure
#plt.show()
plt.savefig("../../paper/figures/two-instruments.png",
    dpi=300, transparent=False, bbox_inches="tight")


#! Conclusion:
# Using both instruments gives the midpoint of the bounds, so does not add info.
