# *Code for:* Accounting for Cellular-Level Variation in Lysis: Implications for Virus-Host Dynamics
*By Marian Dominguez-Mirazo, 2024*

## General Description

In this repository you can find all code needed to replicate the figures and analyses shown in the paper. Most code is written in Matlab R2021_a, except for the Parameter Inference framework, which is written in Julia v1.7.2. Information on the inference framework, including Julia versioning and dependencies, can be found in the README on folder *./ParameterInference. 

## Folder content description

### ./Data

Contains data retrieved from the literature and pre-run data to improve plotting speed. 

- **OneStep_data.csv :** Metadata table containing information on one-step growth curves retrieved from the literature. Required to plot Figure 1. A version of this table is available as Table S1. 

- **OneStep_datasets/ :** Contains the individual datasets for the one-step growth curves retrieved from the literature. Required to plot Figure 1. First column of the file is time, second column is free-virus measurement. The units for each column are shown in the metadata file. 

- **Simulated OneSteps :** One-step growth curve simulations for multiple coefficient of variation (CV) values for the three datasets used in parameter inference (see main text Table 3). The first column corresponds to the simulated CV, the second column to the time of first visible burst. Required to plot Figure 6b. **Origin script: ./FigureScripts/GenerateOneStepTimesforFig6b.m**

- **host_CI_data1_6.mat :** Contains multiple total host simulations using the parameters inferred (see ParameterInference folder) for dataset: data1 and data ID: 6. Required to plot Figure 6a. **Origin plot: ./FigureScripts/GenerateCIforFig6a.m**

- **vir_CI_data1_6.mat :** Contains multiple free-virus simulations using the parameters inferred (see ParameterInference folder) for dataset: data1 and data ID: 6. Required to plot Figure 6a. **Origin plot: ./FigureScripts/GenerateCIforFig6a.m**

- **WangDennehy.csv: **  Lysis time population and cellular-level measurements recovered from Wang 2006, and Dennehy et al., 2011 respectively. See supplementary material Table S2. 

### ./FigureScripts

The folder contains all scripts required to replicate main and supplementary text figures, as well as some pre-run data required for plotting. 

- **Figure1_OneSteps.m**

- **Figure2_IndividualVariationEffect.m**

- **Figure3_WangDennehy.m**

- **Figure4_identifiability.m**

- **Figure5_multicycle.m**

- **Figure6a_inferenceExample.m : ** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how to generate the data. 

- **Figure6b_inferenceAll.m :** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

- **FigureS1_EtoCV.m**

- **FigureS2_TurbidityvsOneStepErrorBars.m** 

- **FigureS3_insilicodata.m :** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

- **FigureS4_otherParams.m : ** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

- **FigureS5_diagnostic.m : ** The data required to plot this figure is obtained from the Parameter Inference framework. See the ParameterInference folder README for information on how the data was generated. 

### ./Figures
Figure storage. Output folder for figure producing scripts found in ./FigureScripts folder. 

### ./ParameterInference
The folder contains all code required to run a parameter inference framework that predicts host and viral-life history traits based on parameter fitting. The code is written in Julia v1.7.2 and contains its own README. 


