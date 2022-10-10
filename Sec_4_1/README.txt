This directory contains the scripts used to generate the results in section 4.1.
It contains these files:

- objects.RData: It contains a list with designs, a matrix with hyperparameters that define a prior scenario, and two sets of ground truths and initial points for adaptive designs.
- functions.R: Functions that are used.
- example.R: An example to run the main functions for EIBV approximations I, II and MCMC. A case to generate one adaptive path is included too.
- prescripted_part.R: Script used to obtain the results of the pre-scripted part.
- results.RData: Results for section 4.1.
- metrics.R: A script to generate the results of Table 1.
- adaptive_part.R: Script used to obtain the results of the adaptive part in Table 2.

We recommend to start with example.R.
To install the required R packages use these commands in R:

install.packages("mvtnorm")
install.packages("magrittr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tibble")
install.packages("parallel")

