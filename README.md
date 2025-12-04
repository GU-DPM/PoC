[12/04/2025]This repository contains the Python (v3.11.9) code for the paper titled "[**A Proof-of-Concept Clinical Trial Design for Evolutionary Guided Precision Medicine for Cancer**](https://www.medrxiv.org/content/10.1101/2025.05.23.25328210v1)" submitted to submitted to medRxiv.
<div align="center">
<img src="/Fig/Fig 1.png" width="400" height="400" title="Graphical depiction of phenotypic states in genetic space illustrating key 606 principles of DPM, along with a schematic of the mathematical model structure">
</div>

The **Code** folder contains the functions used in the work.\
The **Fig** folder contains the figures used in the paper.\
All original data required for running the code is provided in the .csv files within the **./Code/pnas** directory. These files store the original, non-misspecified parameter values.\
**Functions:**\
**(1).DPM_drugmis**\
This function is the main function used to run the analysis. It consists of five steps that generate all the figures and tables in the study.\
*The first step is to generate the misspecified drug 2 sensitivies.\
*The second step is to run the simulations under the different misspecification scenarios.\
*The third step is to run the simulations without misspecification.\
*The fourth step is to process the simulation outputs, covering both misspecification and non-misspecificaiton scenarios.\
*The fifth step is to analyze the results from the fourth step and genrate the raw figures (Fig.3,4,5,S2, Tab.S1) used in the paper.\


