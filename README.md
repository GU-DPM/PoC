[12/04/2025]This repository contains the Python (v3.11.9) code for the paper titled "[**A Proof-of-Concept Clinical Trial Design for Evolutionary Guided Precision Medicine for Cancer**](https://www.medrxiv.org/content/10.1101/2025.05.23.25328210v1)" submitted to submitted to medRxiv.
<div align="center">
<img src="/Fig/Fig 1.png" width="400" height="400" title="Graphical depiction of phenotypic states in genetic space illustrating key 606 principles of DPM, along with a schematic of the mathematical model structure">
</div>

The **Code** folder contains the functions used in the work.\
The **Fig** folder contains the figures used in the paper.\
All original data required for running the code is provided in the .csv files within the **./Code/pnas** directory. These files store the original, non-misspecified parameter values.\
**Functions:**\
**(1).DPM_drugmis.py**\
This function is the main function used to run the analysis. It consists of five steps that generate all the figures and tables in the study.\
i.The first step is to generate the misspecified drug 2 sensitivies.\
ii.The second step is to run the simulations under the different misspecification scenarios.\
iii.The third step is to run the simulations without misspecification.\
iv.The fourth step is to process the simulation outputs, covering both misspecification and non-misspecificaiton scenarios.\
v.The fifth step is to analyze the results from the fourth step and genrate the raw figures (Fig.3,4,5,S2, Tab.S1) used in the paper.\
**(2).DPM_analysis.py**\
This script analyzes the simulation results from the model.\
**(3).DPM_assign_check.py**\
This script assigns and checks the arguments used in the model.\
**(4).DPM_cl.py**\
This script defines the argument format for command-line execution.\
**(5).DPM_constant.py**\
This script defines the model parameter names, input parameter names, and the default values and ranges used in the model.\
**(6).DPM_generate.py**\
This script generates the function arguments and parameters used in the model.\
**(7).DPM_lib.py**\
This script loads the packages required to run the work.\
**(8).DPM_miscellaneous.py**\
This script contains the miscellaneous functions used in the work.\
**(9).DPM_model.py**\
This script contains the model functions and performs the model simulation.\
**(10).DPM_plot.py**\
This script plots the simulation results of the model.\
**(11).DPM_print.py**\
This script defines the print functions used in the work.\
**(12).DPM_read_save.py**\
This script contains functions for reading data from files and saving data into files.\
**(13).DPM_run.py**\
This script contains functions for running the model simulations.\
**(14).DPM_run_drugmis.py**\
This script contains functions for analyzing the simulation results related to drug misspecification.\
**(15).DPM_strategy.py**\
This script contains various treatment strategies.





