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

## Main methods

### 1. Subject-specific jackknife conformal prediction

The jackknife approach uses observations from the target subject to construct an individualized prediction region.

Main functions:

- `jackknife.calib()`: Computes the calibration for a jackknife conformal prediction method for detecting potentially anomalous trajectories
at the individual level. This approach is suitable for multilevel,
temporally structured (functional) data, such as those
commonly encountered in free-living digital health studies.
- `jackknife.pred()`: Performs subject-specific conformal prediction on one or more new trajectories from a target individual
using a previously calibrated jackknife conformal model.

### 2. Pooled conformal prediction

The pooled approach borrows information across subjects. Distributional-similarity weights are used to account for differences between the target subject and the available calibration data.

Main functions:

- `pooled.calib()`: Computes the calibration for a pooled conformal prediction
method that accounts for departures from exchangeability
in multilevel functional data.
- `pooled.pred()`: Performs subject-specific conformal prediction on
one or more new trajectories
from a target individual using a previously calibrated pooled conformal model.

## Additional functions

- `calibration.data()`: generates simulated multilevel functional data.
- `test.data()`: generates test curves for a new subject for prediction.

Consult the function documentation for the exact arguments, return values, and supported data structures.

```r
?jackknife.calib
?jackknife.pred
?pooled.calib
?pooled.pred
?calibration.data
?test.data
```

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

