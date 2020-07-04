# Interactive rank test
This code is an accompaniment to the paper titled "Which Wilcoxon should we use? An interactive rank test and other alternatives", and produce all plots in the paper.

## Overview
The i-Wilcoxon test is for two-sample comparison with unpaired data in a randomized trial that allows human interaction. We also investigate several non-interactive variants of the Wilcoxon signed-rank tests. 

To reproduce the experiments:
```
$Rscript reproduce.R
```
To generate the figures:
```
$Rscript plots.R
```

## Functions and files
tests.R: the i-Wilcoxon test; variants of the Wilcoxon test; linear-CATE-test.

sim_dat.R: simulate outcomes of different types.

result/.: p-values of each method in each experiment.

## Dependencies
The code was tested using R (version 3.6.0) and the following packages:
magrittr, reshape2, ggplot2, randomForest, ash, caTools, tibble, dplyr, tidyr, doParallel, foreach, energy, MASS
