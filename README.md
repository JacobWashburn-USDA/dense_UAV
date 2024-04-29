# High Temporal Resolution UAS Phenotyping Provides Unique Information Between Flight Dates #

If you use this data or scripts please cite: <MANUSCRIPT CITATION>

### Interactive Figures ###

* Figure 1. RmB vegetative index correlations and autocorrelations between: [A)](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig1A.html) flight dates throughout the season, and [B)](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig1B.html) auto-correlations for all pairs of flights and every possible distance between flights across the season in days.
* [Figure 2.](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig2.html) Correlations between standard manual phenotypes taken for each plot at a single time point in the season and the vegetative index GmR from flights throughout the season. The few days before and after the average flowering time are highlighted with a gray box.
* Figure 3. [A)](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig3A.html) Temporal plant height after Weibull fit modeling shown by tester. Gray box represents the mean flowering time plus and minus the standard deviation. PHK76, PHP02, and PHZ51 are the names of tester inbreds used in the study. [B)](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig3B.html) Densities of phenotype values by tester for manual phenotypes including days to anthesis (DTA), days to silking (DTS), anthesis silking interval (ASI), plant height (PHT), ear height (EHT), and grain yield as well as descriptive phenotypes derived from the Weibull fit temporal plant height data including asymptote, inflection point, and growth rate. [C)](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig3C.html) Percent variance explained by both manual and temporal plant height phenotypes. Variance components were derived from the BLUPs model described in section 2.4. Temporal plant height was derived from a similar model but with weibull fit height as the phenotype.
* [Figure 4.](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig4.html) Variance components and repeatability for the VARI and RmB vegetative indices for each flight date individually and all flight dates together (far right).
* [Figure 5.](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig5.html) Prediction ability results from predicting yield using genomic prediction and two different phenomic prediction models: one with all flight dates, the other with only flight dates before flowering. The phenomic models included all calculated vegetative indices as well as UAS derived Weibull fit plant height. All predictions were made using a CV1-like data splitting scenario by holding out (e.g. hidden or untested) a subset of genotypes.
* [Figure 6.](https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig6.html) Yield prediction models including: Genomic prediction (M1) compared to phenomic prediction from the first flight date (M2), the first and second flight date (M3), and subsequent flight dates following the same pattern of adding one new flight date to each model (M4-M44). M2-M44 include all vegetative indices and plant height. Displayed values are the means and the simultaneous Tukey’s HSD 95% confidence intervals for each model. Dotted lines are for ease in visualizing models that are significantly more accurate than the first flight date model (M2) and/or significantly less accurate than the last model. Average flowering time is highlighted by the box.





(https://htmlpreview.github.io/?https://github.com/JacobWashburn-USDA/dense_UAV/blob/main/Figures/Fig3A.html)

### Scripts for running analysis ###

To run this analysis you will need to clone the repository onto your computer and install the dependincies required by each script. Scripts are written in R and Python 3 using the jupyter notebooks interface. Scripts for the analysis can be found in the directories "r_scripts" and "Jupyter_notebooks.
