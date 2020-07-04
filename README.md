# Interactive rank test
This code is an accompaniment to the paper titled "Which Wilcoxon should we use? An interactive rank test and other alternatives", and produce all plots in the paper.

## Overview
The i-Wilcoxon test is a test for two-sample comparison with unpaired data in a randomized trial that allows human interaction.
To reproduce the results:
```
$Rscript batch-experiments.R
```
To generate the figures:
```
$Rscript plots.R
```

## Dependencies
The code was tested using R (version 3.6.0) and the following packages:
magrittr, reshape2, ggplot2, randomForest, ash, caTools, tibble, dplyr, tidyr, doParallel, foreach, energy, MASS
