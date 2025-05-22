# Conformal Prediction for Hierarchical Data
This repository is the official implementation of Conformal Prediction for Hierarchical Data. 

## Requirements

To install the packages:

```setup
Rscript install_packages.R
```

>📋  The experiments are run under [R version 3.6.1](https://cran-archive.r-project.org/bin/windows/base/old/3.6.1/).

## Running Experiments

For the sake of clarity, we distinguish the code that runs a single instance (or job) from the one which runs the entire experiment presented in the article.

### Single run

To perform one instance of the experiment, run this command:

```
Rscript R/study.R --simulation_indice 0 --config "3 12 1000000 1000 1" 
```

>📋 This code can be run on a personnal computer. The "config" parameter can be adapted to correspond to one of the six configurations considered in the article. The arguments of the "config" parameter in the above example are depth = 3, n = 12, T = $10^6$, N= $10^3$ (not used here) and 1, the configuration index. The complete experiment can be run using a "for loop" on the simulation indices in ./R/study.R (but we prefer the parallelized setup described in the following section).

### Complete Experiment

To perform the complete experiment, run this command:

```
sbatch launcher.sh
```

>📋 This code is based on SLURM, thus it require adequate computational ressources (typically a Linux high performance computing environment). The files launcher.sh and run.sh shall be modified to corresponds to your own computational ressources.

## Evaluation

To evaluate the performances of the complete experiment for the different approaches, run this command:

```
sbatch eval.sh
```

>📋 The outputs are, among other things, the figures and two .txt files containing the LaTeX code used to reconstruct the tables.


## Results

We copy here the main experimental results from the article:

### SCP for joint coverage based on ellipsoidal sets

|     Matrix A      | Config. |     Alg. (1)   |     Alg. (2)       |
|-------------------|---------|----------------|--------------------|
| $\widehat{\Sigma}^{-1}$ |    1    | 17.4 ± 0.6     | **16.5 ± 0.55**    |
| $\widehat{\Sigma}^{-1}$ |    2    | 3.36 ± 0.12    | **3.19 ± 0.11**    |

### Component-wise SCP

| Config. | Direct         | OLS            | WLS              | Combi            | MinT             |
|---------|----------------|----------------|------------------|------------------|------------------|
| 1       | 876 ± 254      | 787 ± 226      | 322 ± 131        | 364 ± 101        | **216 ± 47**     |
| 2       | 871 ± 253      | 753 ± 216      | 308 ± 116        | 361 ± 92         | **246 ± 51**     |
| 3       | 3032 ± 467     | 2954 ± 455     | 1869 ± 377       | 1758 ± 395       | **1502 ± 578**   |
| 4       | 3036 ± 479     | 2901 ± 458     | 1581 ± 340       | 1604 ± 349       | **1404 ± 571**   |
| 5       | 10424 ± 885    | 10358 ± 880    | **8861 ± 853**   | 9664 ± 850       | 10613 ± 918      |
| 6       | 10621 ± 889    | 10460 ± 875    | **7673 ± 785**   | 9068 ± 806       | 10503 ± 905      |
