# robustICP

**Robust Individualized Conformal Prediction for Multilevel Functional Data**

<!-- Add badges after confirming the corresponding links and status. -->

## Overview

`robustICP` provides methods for robust individualized conformal prediction for repeated 
    and non-exchangeable data. The package implements two conformal methods, both 
    using robust functional principal component analysis to obtain functional 
    representations: a subject-specific jackknife conformal procedure that relies 
    exclusively on within-individual information, and a pooled conformal procedure 
    that borrows information across individuals while accounting for departures from 
    exchangeability through distributional-similarity weighting. The methods construct 
    individualized prediction regions with user-specified coverage levels and can be used 
    for anomaly detection in multilevel functional data. The package also provides a function to generate simulated multilevel
    functional data based on the data-generating framework described in the
    relevant research literature.

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

