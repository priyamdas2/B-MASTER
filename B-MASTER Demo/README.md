# B-MASTER Demo Code

This repository provides demo MATLAB code (`BMASTER_Demo.m`) for running B-MASTER on paired multivariate datasets `X.csv` and `Y.csv`. The script reads the data, performs a random train/test split, standardizes the training data, runs the MCMC sampler, estimates a sparse coefficient matrix, and evaluates prediction accuracy on the test set using RMSE and MAD.

To run the demo, place `X.csv`, `Y.csv`, and all required MATLAB function files in the same folder, then run the main demo script in MATLAB. For a quick test, `NumIter = 100` can be used; for more stable results, `NumIter = 1000` or larger is recommended.
