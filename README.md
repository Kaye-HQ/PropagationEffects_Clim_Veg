# PropagationEffects_Clim_Veg

![R](https://img.shields.io/badge/R-4.x-blue.svg)
![Rcpp](https://img.shields.io/badge/Rcpp-supported-blue.svg)
![Hydrology](https://img.shields.io/badge/Domain-Hydrology-green.svg)
![Status](https://img.shields.io/badge/Status-Reproducible%20workflow-orange.svg)

This repository contains the code and data used for the study titled "Disentangling the impacts of changes in climate and vegetation on hydrological processes across 2,252 global catchments". The workflow combines an R-based modeling pipeline with a C++ implementation of the HBV-PML model.

> A streamlined utility implementation of this workflow is also available in the hbv.pml R package: https://github.com/Kaye-HQ/hbv.pml

## Overview

This project is designed to:

- calibrate a hydrological model against observed streamflow data,
- evaluate model performance under different calibration schemes,
- separate and quantify the contributions of climate and vegetation changes to hydrological behavior.

## Repository structure

- [code](code): R and C++ source code for the modeling workflow
  - [code/01_model_cal](code/01_model_cal): calibration and effect-separation scripts
    - [code/01_model_cal/functions.R](code/01_model_cal/functions.R): helper functions for data preparation and model execution
    - [code/01_model_cal/HBV_PML.cpp](code/01_model_cal/HBV_PML.cpp): C++ implementation of the model
    - [code/01_model_cal/run_seq.R](code/01_model_cal/run_seq.R): main entry script for the workflow
    - [code/01_model_cal/step1_cal_schemes.R](code/01_model_cal/step1_cal_schemes.R): calibration experiments
    - [code/01_model_cal/step2_effects_seperation.R](code/01_model_cal/step2_effects_seperation.R): effect-separation experiments
- [Data](Data): input data and metadata
  - [Data/info_2252.csv](Data/info_2252.csv): basin metadata as input for the 2,252 catchments
  - [Data/Merged_Par_Full_cal.csv](Data/Merged_Par_Full_cal.csv): merged parameter and calibration data for the 2,252 catchments
  - [Data/sample_data](Data/sample_data): sample meteorological and vegetation inputs
- [Result](Result): generated outputs and experiment results
  - [Result/HBV_PML](Result/HBV_PML): calibration results under different schemes
  - [Result/Result_Exps_detrend](Result/Result_Exps_detrend): experiment outputs for climate and vegetation scenarios

## Requirements

The workflow is implemented in R and requires the following packages:

- data.table
- GA
- Rcpp
- lubridate
- doParallel

A working C++ toolchain is also required for Rcpp compilation.

## Quick start

1. Open the R project file [code/Code.Rproj](code/Code.Rproj) in RStudio.
2. Install the required R packages if they are not already available.
3. Run the main workflow script:
   - [code/01_model_cal/run_seq.R](code/01_model_cal/run_seq.R)

You can also run the main steps separately by sourcing:

- [code/01_model_cal/step1_cal_schemes.R](code/01_model_cal/step1_cal_schemes.R)
- [code/01_model_cal/step2_effects_seperation.R](code/01_model_cal/step2_effects_seperation.R)

## Outputs

The workflow generates results in the following folders:

- [Result/HBV_PML/01_cal_daily](Result/HBV_PML/01_cal_daily): daily calibration results
- [Result/HBV_PML/02_cal_log](Result/HBV_PML/02_cal_log): log-based calibration results
- [Result/Result_Exps_detrend/Exp_data](Result/Result_Exps_detrend/Exp_data): experiment outputs for different forcing and vegetation scenarios

## Notes

- The analysis is computationally intensive because it performs calibration for many basins.
- File paths in the scripts may need to be adjusted when running on a different machine.
- This repository is intended for reproducible research and further exploration of hydrological responses to climate and vegetation change.

## Citation

If you use this workflow or code in your research, please cite the associated study and the repositories accordingly.

- Huang, Q., & Zhang, Y. (2026). PropagationEffects_Clim_Veg: Sample Data and HBV-PML model (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.21253454
- Huang, Q., & Zhang, Y. (2026). Disentangling the impacts of changes in climate and vegetation on hydrological processes across 2,252 global catchments. Water Resources Research, 62, e2025WR043326. https://doi.org/10.1029/2025WR043326

## Contact

Qi Huang  
Email: qihuang@ninhm.ac.cn
