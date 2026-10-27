# robustICP

**Robust Individualized Conformal Prediction for Multilevel Functional Data**

<!-- Add badges after confirming the corresponding links and status. -->

## Overview

`robustICP` provides methods for robust individualized conformal prediction in multilevel functional data settings, particularly when observations are repeated within subjects and exchangeability assumptions may be questionable.

The package implements two complementary approaches:

- **Subject-specific jackknife conformal prediction:** constructs individualized prediction regions using within-subject information.
- **Pooled conformal prediction:** combines information across subjects and uses distributional-similarity weighting to account for differences between observations.

Both approaches use robust functional principal component analysis (robust FPCA) to represent functional observations in a lower-dimensional space. The package also provides functionality for generating simulated multilevel functional data.

Potential applications include individualized uncertainty quantification, prediction-region construction, and anomaly detection in repeated functional measurements.

## Installation

### Install from GitHub

Install the package directly from the GitHub repository:

```r
install.packages("remotes")
remotes::install_github("HolgerCevallos/robustICP")
```

Load the package:

```r
library(robustICP)
```

### Install from CRAN

Once the package is accepted and published on CRAN, it will also be installable using:

```r
install.packages("robustICP")
```

## Main methods

### 1. Subject-specific jackknife conformal prediction

The jackknife approach uses observations from the target subject to construct an individualized prediction region.

Main functions:

- `jackknife.calib()`: performs the calibration procedure.
- `jackknife.pred()`: obtains prediction results using the calibrated procedure.

### 2. Pooled conformal prediction

The pooled approach borrows information across subjects. Distributional-similarity weights are used to account for differences between the target subject and the available calibration data.

Main functions:

- `pooled.calib()`: performs pooled calibration.
- `pooled.pred()`: obtains prediction results using the pooled procedure.

## Additional functions

- `calibration.data()`: generates simulated multilevel functional data.
- `test.data()`: provides additional data-related functionality.

Consult the function documentation for the exact arguments, return values, and supported data structures.

```r
?jackknife.calib
?jackknife.pred
?pooled.calib
?pooled.pred
?calibration.data
?test.data
```

## Basic workflow

A typical analysis follows these steps:

1. Prepare the functional observations and identify the target subject.
2. Select the jackknife or pooled approach according to the data structure and modeling objective.
3. Calibrate the conformal procedure using the appropriate function.
4. Obtain prediction results using the corresponding prediction function.
5. Evaluate the resulting prediction regions and their empirical coverage on suitable validation data.

The precise input format and argument names depend on the selected function. See the function documentation and package examples for executable code.

## Validation and reproducibility

Prediction performance should be assessed using a validation design appropriate for the intended application. In multilevel functional data, repeated curves from the same subject may be dependent; therefore, the evaluation should reflect the intended prediction target and account for the data hierarchy.

Empirical coverage on a finite test set is not, by itself, a guarantee of nominal coverage or a proof of a theoretical coverage property.

## Documentation

The function help pages are available after installing the package:

```r
help(package = "robustICP")
```

For the latest development version, visit the [GitHub repository](https://github.com/HolgerCevallos/robustICP).

## Citation

If you use `robustICP` in academic work, please cite the package and the associated methodological publication when the citation information is available.

**Citation information:** To be completed with the final publication or a verified package citation.

## License

This package is distributed under the GNU General Public License (GPL), version 3 or later. See the `LICENSE` file for details.

## Development status

`robustICP` is under development. Please consult the repository for the current version, documentation, and updates.
