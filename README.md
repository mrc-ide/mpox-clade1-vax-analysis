# mpox-clade1-vax-analysis
This repository reproduces the analysis in the paper: "Epidemiological characteristics and vaccination impact scenario modelling of concurrent Clade I mpox outbreaks in the Democratic Republic of the Congo and Burundi". 

Preprint: https://www.medrxiv.org/content/10.64898/2026.02.24.26346883v1 (under review)
 
# Package installation
This repository uses the package mpoxseir to generate the model. In turn, this repository relies heavily upon odin (a high-level language for implementing mathematical models) and monty (a package used to fit odin models to data). These packages require a compiler to install dependencies for the package and to build any models with odin. Windows users should install Rtools. Be sure to select the “edit PATH” checkbox during installation or the tools will not be found.

After installation of odin, ensure you have the devtools package installed by running the following:

install.packages("devtools")

Then install the mpoxseir package (https://github.com/mrc-ide/mpoxseir) directly from GitHub by running:

devtools::install_github("mrc-ide/mpoxseir")

If you have any problems installing then please raise an issue on the mpoxseir GitHub.

If everything has installed correctly, you then need to load the package:

library(mpoxseir)

You should now be able to run the analysis in this repository.

## Orderly repository
This is an orderly repository. All analyses are organised in a series of interdependent tasks. Please see https://mrc-ide.github.io/orderly/ for more information. 
Install the orderly package by:
remotes::install_github("mrc-ide/orderly", upgrade = FALSE)

# Running the analysis 
Please use the scripts "run_scripts/run_pipeline.R" and "run_scripts/run_pipeline_burundi.R" to reproduce the analyses for DRC and Burundi, respectively. 
Note that the model-fitting task "pmcmc" is highly computationally expensive. We used a high-performance computing cluster to perform this task and would strongly recommend that no one runs this locally on their device. 


