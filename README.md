# **South Africa: Hepatitis B birth dose model**



*Data and code to accompany Cost-effectiveness of the hepatitis B birth dose vaccine in South Africa: a mathematical modelling analysis (under review)*



### Required Python Packages

Atomica: https://atomica.tools/

Sciris: https://docs.sciris.org/en/latest/



### Reproducing Analyses

* The complete analysis can be reproduced by running the ```sa\_hepbd\_model.py``` script 
* Provided are the required framework (outlining model transitions), databook (transition rates, population sizes, intervention impacts etc), and calibration adjustment (y\_) factors to reproduce the analysis. Code to run the calibration is not provided, but YAML files which outline the calibration process are included in the calibration folder.
* Modelled costs in 2024 ZAR are included in a separate csv within the data folder
* Model results are stored as pickle (.pkl) files for post-processing to produce figures and data for tables. 
* Results in output excel files tab2epi and tab2econ contain results consolidated into Table 2 (main text)
* **Outputs will save in a user defined directory**









