# From Spectra to Local Networks: Evaluating MS² Similarities

This repository contains the scripts used to perform the calculations for the study **“From Spectra to Local Networks: Evaluating MS² Similarities.”**

## Repository Structure

- **1_SpectraSimilarity.jl**  
  Contains the similarity calculations originally written in Python by Yuanyue Li  
  (https://doi.org/10.1038/s41592-021-01331-z) and translated into Julia.

- **2_imports.jl**  
  Lists all packages required for the calculations.

- **3_Functions.jl**  
  Includes functions for:
  - filtering the internal database (compiled from MassBank, MoNA, GNPS, and NIST20),
  - generating subdatasets and local networks,
  - computing performance metrics,
  - and other processing steps.

- **4_main.jl**  
  Contains the main analysis workflow, including all calculations and visualizations.


