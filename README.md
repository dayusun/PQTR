# PQTR: Partial Quantile Tensor Regression

MATLAB implementation of **partial quantile tensor regression (PQTR)**, with the
simulation study from the paper.

## What problem it solves

Quantile regression when one predictor is an **array rather than a vector**: a
3-D brain image, a subject-by-region-by-time block, any tensor-valued covariate,
alongside ordinary scalar covariates.

Fitting such a model directly is infeasible, because an unconstrained tensor
coefficient has more free parameters than any study has subjects. PQTR reduces
the tensor predictor first, with a partial-least-squares-type algorithm that
extracts directions relevant to the **conditional quantile** of the outcome, and
then fits the quantile regression in the reduced space. Because the reduction is
quantile-specific, the method recovers coefficient structure that varies across
the distribution rather than only at the mean.

The fitted model is

    Y(tau) = alpha(tau) + beta(tau)' Z + <B(tau), X>

where `X` is the tensor predictor, `Z` the scalar covariates, and `B(tau)` the
tensor coefficient at quantile `tau`.

## Paper

Sun, D., Qiu, Z., Peng, L., Guo, Y., and Manatunga, A. (2024). Partial Quantile
Tensor Regression. *Journal of the American Statistical Association*,
**120**(551), 1724-1735. https://doi.org/10.1080/01621459.2024.2422129

Free full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC12448065/

## Requirements

MATLAB with [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org) v3.6 or
later.

## Installation

Download the function files and add the folder to the MATLAB path.

## Usage

`pqtr` is the main function. The partial-least-squares reduction needs a rank
for the reduced tensor predictor, chosen by one of three methods:

| `method` | Selection rule |
|----------|----------------|
| `'fix'`  | rank supplied by the user |
| `'ER'`   | eigenvalue ratio |
| `'CV'`   | cross-validation |

`simulation.m` reproduces the simulation study. Set `scena`, `src`, `casenum`,
`varcase` and `errdist` at the top of the script to choose the coefficient shape
(tri-square, cross, bi-square, bi-circle, frame), homogeneous or heterogeneous
quantile effects, the covariance model (envelope, compound symmetric) and the
error distribution (t, chi-square, normal).

## R implementation

`pqtr()` is also available in [tensory](https://github.com/dayusun/tensory), an R
package for tensor algebra and tensor regression.

## Author

Dayu Sun, Department of Biostatistics and Health Data Science, Indiana
University School of Medicine. https://www.sundayu.me/
