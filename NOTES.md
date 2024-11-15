# To-do list

Based on recent meetings with Professors

## To-DO:

2. Work with UK Biobank data, for causal PGI analysis
- Use the relatedness file to get parents, and use the imputation package https://github.com/AlexTISYoung/snipar
- Use the school law change birth cohorts as an instrument for education, and/or use another instrument (distance to college, etc.).

3. Write my results into the document as I go.

4. Take the code for parental imputation of EA PGI, and apply it to the 3 mental health scores, starting that project document as I go with the code from this one.

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
