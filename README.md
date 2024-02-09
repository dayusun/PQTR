# Partial Quantile Tensor Regression (PQTR)
This repository contains the MATLAB implementation of PQTR and simulation example codes.
Partial Quantile Tensor Regression is an innovative method for quantile regression with tensor predictors based partial-least-sqaure-type (PLS) algorithm.
The functions and scripts require Tensor Toolbox for MATLAB (>=v3.6).

## Installation

1. Download all the function files 
2. Add the whole folder to the path in MATLAB

## Introduction and usage

`pqtr` is the main function of PQTR. Please use `help pqtr` for more detials of usage.
`pqtr` fits the quantile regression model
$$
    Y(\tau)  = \alpha(\tau)+\beta(\tau)^{T}Z+<\mathcal{B}(\tau),\mathcal{X}>
$$
with PLS-type algorithm.
Here, $Y$ is the response, $Z$ is the low-dimensional vector predictor, $\mathcal{X}$ is the tensor predictor and $\tau$ is the quantile of $Y$ to be fitted.
The PLS-type algorithm requires the dimension of reduced tensor predictor.
`pqtr` provides three methods for specifying the dimesion:

* Fix: Users provide a fixed dimension 
    
        [alpha, gamma, B_tau, ~] = pqtr(X, Y, Z, tau, 'd', [1,1],'dselect',"Fix");

* ER: the eigenvalue ration method

        [alpha, gamma, B_tau, u] = pqtr(X, Y, Z, tau, 'dselect',"ER");

* CV: the cross-validation method

        [alpha, gamma, B_tau, u] = pqtr(X, Y, Z, tau, 'dselect',"CV");

When the dimension `d` is not provided, the algorithm-selected dimesnion would be in the output as `u`.

## Code Example

`simulation.m` provides a comprehesive example on how to simulate the data in various cases and fit the model with `pqtr`. Here is a simple example.


        Y = normrnd[0,1,[100,1]];
        Z = normrnd[0,1,[100,1]];
        X = tensor(normrnd[0,1,[100,16,16,16]]);
        [alpha, gamma, B_tau, u] = pqtr(X, Y, Z, tau);

## Reference

* Brett W. Bader, Tamara G. Kolda and others, Tensor Toolbox for MATLAB, Version 3.6, [www.tensortoolbox.org](www.tensortoolbox.org) 
