# Simulated-Data-Generation-Methods
This repository contains the R source code and simulation framework for the manuscript: "A Computational Framework and Software Tool for Generating Complex Survival Data with Predefined Censoring Rates in Simulation Studies"

## Overview
A simulated data generation method that simultaneously integrates four key functionalities: modeling complex survival time distributions, supporting two types of right censoring, precisely calibrating to predefined censoring rates, and incorporating survival time-related covariates.

## Repository Structure
The repository is organized into two main R scripts:
* "00_core_functions.R": Contains all the core user-defined functions (described below).
* "01_Example_Run.R": Several typical examples to reproduce the article's results.

## Core Functions Description
The framework is driven by the following main functions, defined in "00_core_functions.R":
### 1. Key Parameter θ Solution Functions
* "thetacal()": The primary numerical inversion function. It solves for the exact underlying parameter (θ) of the specified Type III right censoring distribution (e.g., the rate parameter of an exponential distribution or the upper bound of a uniform distribution) required to achieve a predefined total censoring rate. (Note: When both a normally distributed covariate and a binary covariate are present, kernel density estimation is employed to address the covariate effects).
* "thetacal_a()": An advanced variant of the parameter solution function. (Note: When either a normal distribution covariate or a binary covariate is present, kernel density estimation is employed to address the covariate effect).

### 2. Data Generation & Validation
* "datagen()": The survival data generator with predefined censoring rate. Utilizing the parameter θ obtained from the functions above, it generates realistic survival datasets with right-censoring, alongside corresponding covariates and survival times.
* "censoringcal()":The main simulation loop function. It embeds 'datagen()' to perform 1,000 repeated simulations under specific conditions. It directly calculates and outputs the final empirical censoring rate for validation.

## Getting Started
To reproduce the methodology:
1. Clone or download this repository.
2. Open the R project or set your working directory to the cloned folder.
3. Open "01_Example_Run.R" in RStudio.
4. Run the script. It will automatically source the functions from "00_core_functions.R" and execute specific examples.

## Citation
If you utilize this code or the proposed methodology in your own research, please consider citing our manuscript:
> Chen, H., Yao, Q., Lin, Y. et al. A computational framework and software tool for generating complex survival data with predefined censoring rates in simulation studies. BMC Med Res Methodol (2026). https://doi.org/10.1186/s12874-026-02843-y
        
        
        
        
        
        
        
        
        
