# article-disease-pred-with-kcde

[Current PDF](https://github.com/reichlab/article-disease-pred-with-kcde/raw/master/inst/article/infectious-disease-prediction-with-kcde.pdf)

This repository contains code and data used in writing the paper "Infectious Disease Prediction with Kernel Conditional Density Estimation". It is organized as an R package, but formal installation as such is not necessary. Installation of 3 other R packages hosted on github is either helpful or required:

1. https://github.com/reichlab/kcde (required): A package implementing estimation of KCDE models, oriented towards time series
2. https://github.com/reichlab/pdtmvn (required): A package implementing partially discretized truncated multivariate normal distributions.
3. https://github.com/reichlab/mvtnorm-mod-kcde (optional): A version of the mvtnorm package for R, modified for speed improvements. Installing this instead of the mvtnorm package will result in approximately 30% reductions in run time for some kcde specifications.

This repository is organized as follows:
* data-raw/ contains "raw" data in the form of .csv files. These data were obtained from competition administrators (http://dengueforecasting.noaa.gov/) or the R package providing data on inidence of flu in the United States (https://cran.r-project.org/web/packages/cdcfluview/index.html, http://www.cdc.gov/flu/weekly/).
* inst/ contains everything else:
    * inst/code contains code used to estimate models and make predictions:
        * inst/code/estimation/ contains code for estimation:
            * inst/code/estimation/kcde-estimation-step.R does KCDE estimation for a particular combination of factors describing an estimation task (e.g., KCDE model specification, data set, prediction horizon, etc.)
            * inst/code/estimation/submit-cluster-job-kcde-estimation-step.R sets up jobs on our cluster to run kcde-estimation-step.R for all relevant combinations of factors describing an estimation task.
            * inst/code/estimation/copula-estimation-step.R does copula estimation given KCDE fits obtained by running kcde-estimation-step.R
            * inst/code/estimation/sarima-estimation.R estimates a SARIMA model.
            * inst/code/estimation/surveillance-estimation.R estimates an HHH4 model from the surveillance package.
        * inst/code/prediction/ contains code for prediction and model evaluation:
            * inst/code/prediction/kcde-evaluation-simstudy.R evaluates the quality of density estimates obtained in the simulation study
            * inst/code/prediction/kcde-prediction.R makes predictions for incidence in individual weeks from KCDE in the applications to dengue and influenza, and obtains log scores
            * inst/code/prediction/kcde-peak-prediction.R makes predictions for peak week timing and incidence from KCDE+copulas in the applications to dengue and influenza
            * inst/code/prediction/sarima-prediction.R makes predictions for incidence in individual weeks from SARIMA in the applications to dengue and influenza, and obtains log scores
            * inst/code/prediction/sarima-peak-prediction.R makes predictions for peak week timing and incidence from SARIMA in the applications to dengue and influenza
            * inst/code/prediction/surveillance-prediction.R makes predictions for incidence in individual weeks from an HHH4 model in the application to dengue, and obtains log scores
            * inst/code/prediction/surveillance-peak-prediction.R makes predictions for peak week timing and incidence from an HHH4 model in the application to dengue
        * inst/code/postprocessing/influenza-results-model/ contains a sketch of a start at a spline model to examine the results of the application to influenza from different models. This didn't make it into the paper.
        * inst/code/sim-densities-sim-study-discretized-Duong-Hazelton.R contains the code used in the simulation study to simulate from and evaluate the distributions data were simulated from.
    * inst/results/ contains intermediate results from the applications and the simulation study: model fits and data frames with summaries of the predictions that were made
    * inst/article/ contains the .Rnw and related files for the article.
    * inst/article-supplement/ contains the .Rnw and related files for the supplement.
