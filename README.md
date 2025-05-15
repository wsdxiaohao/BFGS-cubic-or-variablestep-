# Quasi-Newton Methods without Line Search: Reproducibility

This repository contains the code used to reproduce the experiments in our NeurIPS 2025 submission on non-asymptotic global convergence of BFGS methods without line search.

## Overview

Two experiments (log-sum-exp problem and logistic regression) are included.

We provide implementations of the following optimization algorithms:
- Cubic BFGS Quasi-Newton method
- Variable Step BFGS Quasi-Newton method
- Baselines: Gradient Descent (GD), Cubic Newton, Gradient-Regularized Newton, and Heavy Ball Method (HBF)

## Reproducing

To reproduce the experiment, please open and run the `test.ipynb` in the corresponding folder.
## Requirements

Python 3.8.8
NumPy 1.22.4

## License

To be specified after acceptance. Usage permitted for reproducibility and review purposes only.
