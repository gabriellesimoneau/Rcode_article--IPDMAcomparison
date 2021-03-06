# Rcode_article--IPDMAcomparison
R code and support for article "A comparison of bivariate, multivariate random-effects, and Poisson correlated gamma-frailty models to meta-analyze individual patient data of ordinal scale diagnostic tests"

## Description

As the 9-item Patient Health Questionnaire (PHQ-9) individual patient dataset can not be publish, we provide an example dataset to reproduce the data analysis in the article. 

Data:
* `simulated_dataset.csv` is an example dataset that mimics the PHQ-9 IPD dataset.
* `simulation_scenarios.csv` contains information on the parameters that are varied across the 16 simulation scenarios.

R code:
* `tables_figures.R` reproduces Table 1-2-3-5 and Figure 1 using the example dataset, as well as Table 6 and Figure 2-3 in the article
* `simulation_bivariate.R` reproduces the simulations via the Bivariate Random Effects Model (BREM or Bivariate) as described in the article.
* `simulation_poisson.R` reproduces the simulations via the Poisson correlated gamma-frailty model (Poisson) as described in the article.
* `functions_general.R` contains useful functions to set up a dataset to apply the BREM or Poisson model, and a function used to simulate correlated diagnostic data as described in the article.
* `functions_poisson.R` contains useful functions to estimate all parameters of the Poisson model.

Results:
* `simulation_results_poisson.csv` contains the results from the 16 simulation scenarios analyzed with the Poisson model.
* `simulation_results_bivariate.csv` contains the results from the 16 simulation scenarios analyzed with the BREM.
* `parametric_boot_example.csv` contains the results from the parametric bootstrap used in the example data analysis. 
* `Results_example_dataset.pdf` presents the resulting tables and figures as obtained using the example dataset.
* `Results_example_dataset.Rnw` is the corresponding R Sweave file.

The code has been written using `R-3.0.2` (platform: x86_64-apple-darwin10.8.0, 64-bit) with 
package versions `mvtnorm_1.0-2`, `rje_1.9`, `copula_0.999-13`, `MASS_7.3-37` and `lme4_1.1-7`. 

## Data Analysis

Reproduce Table 1-2-3-5 and Figure 1 as shown in `Results_example_dataset.pdf` by running all R code in the Data Analysis section of `tables_figures.R`.

## Simulations

Reproduce Table 6 and Figure 2-3 in the manuscript by running all R code in the Simulation section of `tables_figures.R`.

## Support

For questions, comments or remarks about the code please contact [Gabrielle Simoneau] (mailto:gabrielle.simoneau@mail.mcgill.ca).