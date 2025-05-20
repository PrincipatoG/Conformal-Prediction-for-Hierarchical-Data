# Conformal-Prediction-for-Hierarchical-Data
This repository is the official implementation of [Conformal Prediction for Hierarchical Data](https://arxiv.org/abs/2411.13479). 

## Requirements

To install the packages:

```setup
Rscript install_packages.R
```

>📋  The experiments are run under R version 3.6.1 

## Run of the Experiments

For the sake of clarity, we distinguish the code that runs a single instance of the code from the one which runs the entire experiment presented in the article.

### Single run

### Complete Experiment
To train the model(s) in the paper, run this command:

```train
python train.py --input-data <path_to_data> --alpha 10 --beta 20
```

>📋  Describe how to train the models, with example commands on how to train the models in your paper, including the full training procedure and appropriate hyperparameters.

## Evaluation

To evaluate my model on ImageNet, run:

```eval
python eval.py --model-file mymodel.pth --benchmark imagenet
```

>📋  Describe how to evaluate the trained models on benchmarks reported in the paper, give commands that produce the results (section below).

## Results

We copy here the results from the article:

### SCP for joint coverage based on ellipsoidal sets

|     Matrix A      | Config. |     Alg. (1)   |     Alg. (2)       |
|-------------------|---------|----------------|--------------------|
| $\hat\Sigma^{-1}$ |    1    | 17.4 ± 0.6     | **16.5 ± 0.55**    |
| $\hat\Sigma^{-1}$ |    2    | 3.36 ± 0.12    | **3.19 ± 0.11**    |

### Component-wise SCP

| Config. | Direct         | OLS            | WLS              | Combi            | MinT             |
|---------|----------------|----------------|------------------|------------------|------------------|
| 1       | 876 ± 254      | 787 ± 226      | 322 ± 131        | 364 ± 101        | **216 ± 47**     |
| 2       | 871 ± 253      | 753 ± 216      | 308 ± 116        | 361 ± 92         | **246 ± 51**     |
| 3       | 3032 ± 467     | 2954 ± 455     | 1869 ± 377       | 1758 ± 395       | **1502 ± 578**   |
| 4       | 3036 ± 479     | 2901 ± 458     | 1581 ± 340       | 1604 ± 349       | **1404 ± 571**   |
| 5       | 10424 ± 885    | 10358 ± 880    | **8861 ± 853**   | 9664 ± 850       | 10613 ± 918      |
| 6       | 10621 ± 889    | 10460 ± 875    | **7673 ± 785**   | 9068 ± 806       | 10503 ± 905      |
