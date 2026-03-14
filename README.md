# direct-indirect-education

Senan Hogan-Hennessy, 12 March 2026.

## The Direct and Indirect Effects of Genetics and Education.

Genes associated with educational attainment causally improve labour market income, but the economic mechanism behind this relationship remains poorly understood.
I use quasi-random variation in genetic inheritance to estimate the causal effect of the Education PolyGenic Index (Ed PGI) on education and labour market income, using data from the UK Biobank, controlling for imputed parental Ed PGI values constructed from sibling genetic data to isolate the random component of inheritance.
A one standard deviation increase in the Ed PGI increases completed education years by 0.5 years and raises later life income by around 5 per cent, replicating the main estimates in Carvalho (2025).
I decompose these total genetic effects into a direct genetic channel and an indirect channel operating through education years, combining correlational estimates of the returns to education with a sensitivity analysis benchmarked against the distribution of quasi-experimental estimates from the economics literature for Britain.
At correlational returns of around 6 per cent, roughly 65 to 75 per cent of the total genetic income effect operates through the education years channel; under higher values for returns to education from the economics literature, the direct genetic effect becomes indistinguishable from zero.
Education-linked genes earn their labour market return primarily by inducing more education years, not through direct biological pathways that bypass education.
The institutions governing access to education are therefore the primary mechanism through which genetic endowment shapes income inequality, and the primary policy lever through which that inequality could be addressed.

- [direct-indirect-education-2026.pdf](https://github.com/shoganhennessy/direct-indirect-education/blob/main/direct-indirect-education-2026.pdf)  is the latest version of the working paper.

## Replication

This folder is the replication package for the paper, with statistical analysis on UK Biobank data using the programming language R and associated packages. UK Biobank data access is available via application to the UK Biobank (Project #335854); the analysis code can be run on approved data following the instructions below.

### Data

UK Biobank data are not sharable online, so the data folder is not shared here.

After being granted access to Uk Biobank data, then the data folder contains code to extract data files to a tabular format.

### Programs

Folder "programs/" contains all analysis scripts in the R language.

- "run.sh" sequentially calls all relevant scripts for complete replication of the paper.
- "data-build/" contains R scripts that construct analysis data from raw files.
- "analysis/" contains R scripts that produce all tables and figures presented in the paper.

### Text

Folder "text/" contains all files for the final paper, including the LaTeX source and compiled output.
