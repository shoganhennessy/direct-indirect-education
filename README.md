# direct-indirect-education

Senan Hogan-Hennessy, most recent update 15 November 2024

## Explanation

This folder contains all relevant files (not data files) for my academic project ``The Direct and Indirect Effects of Genetics and Education.''

# To-do list

Based on recent meetings with Professors

## Current to-do list in writing the paper:

1. Connect the OCC score wages in pheno data, provided by Kweon+ (2025).

- Update my UKB data dispense https://dnanexus.gitbook.io/uk-biobank-rap/getting-started/data-structure/updating-dispensed-data

- Update the phenotype extract variables, to use multiple instances (i0, i1, ...) of SOC and education qualifications data.
- Get separate columns for each qualification in EdQuals
- Get the following other PGIs in https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=302
Asthma                  https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26211
Bipolar Disorder        https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26215
Height                  https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26241
BMI                     https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26217
Schizophrenia           https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26276
Type 2 Diabetes         https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26286
Coronary Artery Disease https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=26228

Eventually:
-> calculate the ADHD PGI, which would be an important control in education.
-> Also use the Tan+ (2024) family-based GWAS Ed PGI, as a robustness check (possibly new measure) in the analysis (do total effect estimates change?).
-> Calculate a second Ed PGI (weights from a different sample), to use in an obviously related IV analysis.  This adjusts for measurement error in the PGI weights.

2. Data section

- Summary tables as follows:

      | UKB sibling sample | UKB entire sample
----------------------------------------------
Var 1 | Mean               | Mean             
      | (SD)               | (SD)             
...                                            
Var N | Mean               | Mean             
      | (SD)               | (SD)             
      |                    |                  
      | Obs                | Obs              

- Keep the HRS side-by-side figures, and update to new UKB data.

3. Write the randomisation design section, with figures for the new differences measure
- Connect to the formula instrument (Z - E{Z | parents} ~ Z_i) idea, showing ``first-stage'' results with estimates ~ 1.
- (summer) eventually to be rewritten with the one copy inheritance measure to replace the Z - mean{parents, siblings} measure.

Get a figure showing Ed PGI is correlated with other things (including other PGIs), but the differences measure is not.

4. Write the total effect estimates, Z -> Ed years, Income.

Table 2: total genetic effects.
- Raw OLS
- Raw OLS + controls
- OLS with differenced measure
- 2SLS with differences measure instrumenting the Ed PGI

Figures which show the OLS slope vs causal slope for Ed Years + Income.

5. Write the causal framework for mediation (draw heavily from my other paper), noting the non-identification of direct/indirect effects if selection-into-education.

Write one table of the first-stage IV, School leaving year change + education tuition fees for ed years (note: these are both year of birth-based).

6. Show the mediation results, side-by-side in a table:
- OLS
- OLS + controls
- Heckman selection model
- Semi-parametric selection model

2-3 paragraphs writing on what the mediation results actually mean.

### Newest Abstract Ideas.

I use Mendelian independent inheritance, where a child inherits some genes from parents and not other due to chance, to estimate causal effects.
This approach weights heavily children with greater differences from their parents, similar to a formula instrument approach.
RESULTS ON WAGES.
Next, I use a selection model approach


Main theme: putting labour economic insight into genetic effects
- Outlining where randomisation comes in, adjusting for choice in first- and second-stages
- Economic insight for what the Mendelian causal effects mean
- Implications for policy

Question for policy: if genetic effects are mostly direct, then they are immutable 

If mostly indirect, implies the education system has a big endowment-advantage, and are not equal opportunity.
ANyone can study; not just anyone can succeed, making this an equality of opportunity issue.


natural variation in the Ed PGI.
This approach is 

### Comments from Labour Seminar

1. From Tak: let the questioners finish speaking, and do not interrupt them.
He thought the ``deviation away from parents'' figure was a very compelling aha moment, and thus should be as soon in the presentation as possible.
The DAG shows no interaction effects, so he asked whether this had interaction effects; this has come up many time, so it may be important to write out the direct and indirect effects, and talk through the ``non-restrictions on interactions of same treatment effects.''
If I include a quote from somewhere else, then I (the presenter) must read it verbatim, and not allude to it + confuse the now-multitasking audience.
He thought ht efocus was backwards, and should be what genes are about intelligence?  Are they the ones in the education index, or for other things, etc.

2. From Zhuan Pei:
Needs a punch statement for why decompose.  Currently lacking one.
Careful about saying no previous natural experiment mediation (e.g., Kling Lietman 2007 Metrica).
``Natural experiment for what?''
Big fan of the formula instrument connection, and notes other labour economists will be familiar with that reasoning, and thus presentation in that may be useful.
Senan expansion: posing the formula instrument of Z_random := Z_child - 1 / 2 (Z_father + Z_mother) is a different approach to Cavhalo (2024), Young+ (2022), Biroli+ (2024) who instead do a linear regression Y = \beta Z_child + phi_f Z_father + phi_m Z_mother --- i.e., linear controls for parent values.
Give this some extra thought, as using this (and explicitly acknowledging the measurement error in each component) could be an avenue to distinguish from those previous papers.

3. From Evan Riehl.
He was a big fan of the figures for the random deviation away from parents, and encouraged focusing more on the total effect (i.e., before even going into the mediation).
This could be fruitful if it is plausibly different from Cavhalo's work; at the moment it seems too close.
Figure out if the formula instrument leads to meaningful differences; it could be the case that OLS with linear controls for mother and father relies heavily on the extrapolation argument, while the instrument is a purer design-based proposition.

4. From Senan:
If the parental deviation formula IV yields distinct results from Cavhalo (2024), then this paper should focus more on total effects, and a spun out version on the metrics of mediation becomes a different project, ``causal mediation under selection'' with a Roy model and control function solution.
It feels like the metrics work for mediation methods is shoe-horned into this paper.
Using the UK Biobank data in a section of that paper is still possible, too.
This could then accompany a paper singling out the Vietnam draft IV, too.

5. From Doug Miller, March 2025.
- External literature on returns to education show that OLS estimates are, in fact, not far off real (despite concerns for Roy-style ability selection).  A slide showing the CM results OLS + selection model, with a button for OLS vs selection-model education returns (to illustrate how I am over-coming the identification problem).
- Consider a whisker plot, summarising correlation of Ed PGI with other PGIs, then parent differenced one no longer correlated with Ed PGI.
- Recommends eventually controlling for other PGIs, in the analysis. 

## General notes from faculty meetings.

### Notes from Meeting with Ben Goldman $+$ Zhuan Pei

- Be very clear and precise about the research question: what am I asking, which has a clear answers?  What is the ideal experiment?
- What is important, in economic terms, about learning the decomposition of direct and indirect effects of education genetics?
How would this influence policy, understanding, our view of society, etc.?
- Be very precise on estimand of interest, and show in a presentation a simple version of this.
- A very open question on what role genetic analyses do when genes are only important up to their interactions with environment.
For example, a genetic disability in the 1950s may hinder education take-up (despite high returns to education) leading to high direct genetic effect estimates; this is levelled out after the passage of the ADA in 1970s, despite nothing changing about genetics.
In this sense, genetic effects are not distinguishable from effects of discrimination; consider this story in the case of Educ PGI negatively correlated with ADHD PGI.
- In the real world, $Y_i | Z_i = 0, D_i = 1$ where low Educ PGI people don't ever take-on education, might barely exist.
This limits the practical relevance of a mediation exercise.
- Variation in education must be sufficiently large for this to work; a one year of extra schooling from a new law might not cut this.
Especially if many of the gains are from higher education; maybe I should look into costs of education, too.
- Remember the story of the Earnings PGI being highly correlated with Educ PGI $\implies$ compelling in explaining that measurable genetic differences.
- See Florens Heckman (2008) on justifying the control function assumptions.

5. Writing tips for the introduction:
https://blogs.ubc.ca/khead/research/research-advice/formula
https://sites.google.com/site/amandayagan/writingadvice?authuser=0
