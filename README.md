# Dependence Analysis between Investment Funds and Macroeconomic Variables

## About the Project 

This repository contains the code used in the dissertation "Assessment of the dependency structure between Investment Funds and macroeconomic variables: An application of copula theory" developed at UERJ. The project analyzes the relationship between Brazilian investment funds and macroeconomic variables through three main approaches:

1. Cointegration Analysis (ARDL)
2. Causality Analysis (Toda-Yamamoto)
3. Copula Theory

## Theoretical Background

This project combines three main theoretical approaches to analyze the relationship between investment funds and macroeconomic variables:

1. **Copula Theory**: Used to model complex dependency structures between variables, allowing for non-linear relationships and tail dependencies.

2. **Cointegration Analysis**: Employs ARDL bounds testing approach to investigate long-run relationships between variables that may be integrated of different orders.

3. **Causality Testing**: Utilizes the Toda-Yamamoto procedure to examine causal relationships while avoiding pre-test bias.

For detailed theoretical explanations and mathematical foundations, please refer to the documentation in the `docs/` directory.

## Project Structure

```
.
├── R/
│   ├── copula.R             # Functions for copula analysis
│   ├── cointegration.R      # Functions for ARDL cointegration tests
│   └── causality.R          # Functions for causality tests
├── data/
│   └── example_data/        # Example datasets
├── docs/
│   └── dissertation.pdf     # Full dissertation
├── tests/
│   └── testthat/           # Unit tests
├── LICENSE
└── README.md
```

## Main Features

### 1. Copula Analysis (`copula.R`)
- Estimation of bivariate copulas
- Goodness-of-fit tests
- Data simulation
- Calculation of dependence measures
- Visualization of dependence structures

### 2. Cointegration Analysis (`cointegration.R`)
- Unit root tests (ADF)
- ARDL model estimation
- Bounds Testing
- Error correction analysis

### 3. Causality Analysis (`causality.R`)
- Implementation of Toda-Yamamoto test
- Modified Granger causality analysis
- Robustness tests

## Requirements

```r
# Required packages
install.packages(c(
  "tidyverse",
  "VineCopula",
  "fitdistrplus",
  "kdecopula",
  "urca",
  "vars",
  "collapse",
  "aod"
))
```

## How to Use

### Copula Analysis Example

```r
source("R/copula.R")

# Load and prepare data
data <- read_data("path/to/data.csv")

# Find best distribution
dist_info <- find_best_distribution(data, distributions)

# Estimate copulas
copulas <- estimate_bivariate_copula(data, family)

# Simulate data
simulated_data <- simulate_bivariate_copula(copulas, orig_data, dist_info, 1000)
```

### Cointegration Analysis Example

```r
source("R/cointegration.R")

# ADF test
adf_results <- adf_test(data)

# ARDL estimation
ardl_model <- estimate_ardl(data_list, max_lag = 3)

# Bounds test
bound_test_results <- run_bound_test(ardl_model)
```

## Documentation and Theoretical Background

The theoretical foundation of this project is documented in two main files in the `docs/` directory:

1. `dissertation.pdf` - Contains detailed theoretical background on:
   - Copula theory and its applications
   - Different copula families
   - Estimation methods
   - Goodness-of-fit tests
   - Simulation procedures

2. `impact_funds.pdf` - Provides comprehensive information about:
   - ARDL cointegration methodology
   - Toda-Yamamoto causality testing
   - Application to investment funds
   - Empirical results and interpretation

These documents provide the theoretical framework and mathematical foundations for the implementations in this repository.

## Technical Details

### ARDL Cointegration Testing
The Autoregressive Distributed Lag (ARDL) approach to cointegration is used to test for long-run relationships between variables. This method has several advantages:
- Can be applied when variables are integrated of different orders
- Suitable for small sample sizes
- Allows for different lag lengths in the model

The bounds testing procedure involves:
1. Testing for unit roots to ensure no I(2) variables
2. Estimating the ARDL model
3. Computing the F-statistic for the bounds test
4. Testing for serial correlation, normality, and heteroscedasticity

### Toda-Yamamoto Causality Testing
The Toda-Yamamoto approach to testing for Granger causality:
- Robust to different orders of integration
- Does not require pre-testing for cointegration
- Minimizes risks of pre-test bias

The procedure involves:
1. Determining the maximum order of integration (dmax)
2. Selecting the optimal lag length (k)
3. Estimating a VAR model in levels with k + dmax lags
4. Performing modified Wald tests

## Contributing

Contributions are welcome! Please read CONTRIBUTING.md for details on our code of conduct and the process for submitting pull requests.

## Citation

If you use this code in your research, please cite:

```bibtex
@mastersthesis{arruda2024assessment,
  title={Assessment of the dependency structure between Investment Funds and macroeconomic variables: An application of copula theory},
  author={Arruda, Gabriel de Almeida},
  school={Rio de Janeiro State University},
  year={2024}
}
```

