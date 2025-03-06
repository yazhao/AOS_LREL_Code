---
output:
  html_document: default
  pdf_document: default
---
# Testing High-Dimensional Regression Coefficients in Linear Models

Assembled code for Zhao et al (2024) published in the Annals of Statistics

Paper citation: Alex Zhao. Changcheng Li. Runze Li. Zhe Zhang. "Testing high-dimensional regression coefficients in linear models." Ann. Statist. 52 (5) 2034 - 2058, October 2024. https://doi.org/10.1214/24-AOS2420 

Open access link: https://scholarsphere.psu.edu/resources/231ddcf5-36b8-40bf-b08a-d453b00aaf25

## 1. Path setup and data source
      Path setup: open the `AOS_LREL_Code.Rproj` file, which will make sure that you are in the correct relative working directory.

      Data source: The data for the real data analysis example comes from the transNOAH breast cancer trial (GEO series GSE50948), which can be found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50948

      NOTE: there is no need to download this file directly: the real data example code script uses an R library GEOQuery that automatically can pull down the file. If you would like to download the file to not repeat this process, it can be done, although the first run of the real data example script should also download and save the data object to a .RData file.

## 2. Instructions for simulations

      Note: each of these stages will take at least hours to fully run. For convenience a full set of simulation results are included.

      Step 1: open the `AOS_LREL_Code.Rproj` file, which will make sure that you are in the correct relative working directory.

      Step 2: Generate the null distribution results for our new method:
        Part A: Open `simulation_code\simulations_densities_lrel_new.R` in RStudio.
        Part B: Run the full script with cov_rho = 0.25
        Part C: Change cov_rho to 0.5, repeat part B, then do this again by setting cov_rho = 0.75

        Note: only the "knownnew_*" csv files are necessary even though an "unknownnew_*" set of csv files will be generated as well.

      Step 3: Generate the power results for our new method:
        Part A: Open `simulation_code\simulations_power_lrel_new.R` in RStudio.
        Part B: Run the full script with cov_rho = 0.25 and mdr = md25r.
        Part C: Repeat part B with the settings cov_rho = 0.5 and mdr = md5075r, and then cov_rho = 0.75 and mdr = md5075r, respectively.

        Note: there are four settings (normal errors and full signal, normal errors and random uniform signal, t errors and full signal, and t errors and random uniform signal). Each takes over 6 hours to run.

      Step 4: Generate the null distribution results for the Zhong and Chen method:
        Part A: Open `simulation_code\simulations_densities_zc.R` in RStudio.
        Part B: Run the full script with cov_rho = 0.5
        Part C: Repeat part B but change cov_rho to 0.25 and then 0.75

      Step 5: Generate the power results for the Zhong and Chen method:
        Part A: Open `simulation_code\simulations_power_zc.R` in RStudio.
        Part B: Run the full script with cov_rho = 0.25 and mdr = md25r.
        Part C: Repeat part B with the settings cov_rho = 0.5 and mdr = md5075r, and then cov_rho = 0.75 and mdr = md5075r, respectively.

        Note: there are four settings (normal errors and full signal, normal errors and random uniform signal, t errors and full signal, and t errors and random uniform signal). Each takes over 36 hours to run.

      Step 6: Generate the null distribution results for the Cui et al. method:
        Part A: Open `simulation_code\simulations_densities_cui.R` in RStudio.
        Part B: Run the full script with cov_rho = 0.5
        Part C: Repeat part B but change cov_rho to 0.25 and then 0.75

        Note: this script will generate a "*_null_cui_*" and "*_null_cuircv_*" .csv file for each error type and setting, only the "*_null_cuircv_*" csv files are used for comparison purposes.

      Step 7: Generate the power results for the Cui et al. method:
        Part A: Open `simulation_code\simulations_power_cui.R` in RStudio.
        Part B: Run the full script with cov_rho = 0.25 and mdr = md25r.
        Part C: Repeat part B with the settings cov_rho = 0.5 and mdr = md5075r, and then cov_rho = 0.75 and mdr = md5075r, respectively.

        Note: there are four settings (normal errors and full signal, normal errors and random uniform signal, t errors and full signal, and t errors and random uniform signal). Each takes over 15 hours to run.

      Step 8: Generate the sensitivity results for our new method:
        Part A: Open `simulation_code\simulations_sensitivity_lrel.R` in RStudio.
        Part B: Run the full script.

## 3. Real data analysis instructions

    Step 1: Open transNOAH_example.R in RStudio. Make sure that you have opened the `AOS_LREL_Code.Rproj` file so that the relative path is correct to find the `simulation_code` folder.

    Step 2: Ensure that the following R libraries are installed:
      - GEOquery
      - data.table
      - compiler
      - Rfast
      - SIS

      Note: while data.table, compiler, Rfast, and SIS can be installed from CRAN, GEOquery should be installed from BioConductor: https://bioconductor.org/packages/release/bioc/html/GEOquery.html

      Note 2: Rfast requires Rcpp, so it will also require Rtools on Windows: https://cran.r-project.org/bin/windows/Rtools/index.html

    Step 3: Once the path has been set correctly, you can run the entire script, or run the R code line by line.

## R code in this folder:

### `simulation_code` folder:

1. `simulations_functions.R` stores all of the functions that will be used for all of the simulations and real data analysis.

2. `simulations_densities_lrel_new.R` stores all the scripts to generate the null distribution case for our new method.

3. `simulations_densities_zc.R` stores all the scripts to generate the null distribution case for the Zhong and Chen method.

4. `simulations_densities_cui.R` stores all the scripts to generate the null distribution case for the Cui et al. method.

5. `simulations_power_lrel_new.R` stores all the scripts to generate the power results for our new method.

6. `simulations_power_zc.R` stores all the scripts to generate the power results for the Zhong and Chen method.

7. `simulations_power_cui.R` stores all the scripts to generate the power results for the Cui et al. method.

8. `simulations_sensitivity_lrel.R` stores the scripts to generate the sensitivity results for our method.

### other R code:

1. `transNOAH_example.R` is the R file for the real data analysis. This script will generate the file `transNOAH.rData` which stores the initial GEO file for the transNOAH breast cancer trial data.

## Other files:

### results folder:
  Contains the results output from simulation code, all stored in the distance subfolder.

AOS_LREL_Code.RProj: R project file to set relative paths

## git repository files
.git folder: folder for git repository purposes, can be ignored
.gitignore: file to note what files to ignore for GitHub
LICENSE: denotes the license type for git repository
README.md: markdown file (this readme file you are reading)
README.html: HTML version generated by the README.md file
README.pdf: PDF version generated by the README.md file
