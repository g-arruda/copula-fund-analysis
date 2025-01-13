# Cointegration and Causality Analysis Functions
# Author: Gabriel de Almeida Arruda
# Description: Functions for performing cointegration and causality tests on time series data

library(tidyverse)
library(urca)
library(vars)
library(collapse)
library(aod)

#' Perform Augmented Dickey-Fuller (ADF) Test
#'
#' This function performs the ADF test on multiple variables in a dataset to check
#' for stationarity. If non-stationary variables are found, it performs the test
#' again on their first differences.
#'
#' @param data A dataframe containing the variables to be tested
#' @return A list containing the test results for both levels and first differences
#' @examples
#' # Example with a dataframe containing variables x and y
#' data <- data.frame(x = rnorm(100), y = cumsum(rnorm(100)))
#' adf_test(data)
#'
adf_test <- function(data) {
  # Perform ADF test for all numeric variables
  adf_test <- data |>
    select(where(is.numeric)) |>
    map(~ tseries::adf.test(.x, k = 1) |> suppressWarnings())

  # Initialize empty dataframe for results
  df <- data.frame(
    variable = character(),
    adf = numeric(),
    p_value = numeric()
  )

  # Fill results dataframe
  for (i in seq_along(adf_test)) {
    df <- rbind(df, data.frame(
      variable = names(adf_test)[i],
      adf = round(as.numeric(adf_test[[i]]$statistic), digits = 3),
      p_value = adf_test[[i]]$p.value
    ))
  }

  # Add result column based on p-value
  df$result <- case_when(
    df$p_value >= 0.05 ~ "non-stationary",
    df$p_value < 0.05 ~ "stationary"
  )

  # Test first differences if non-stationary variables exist
  if (any(df$p_value >= 0.05)) {
    non_stationary_vars <- df |>
      filter(result == "non-stationary") |>
      pull(var = "variable")

    # Perform ADF test on first differences
    adf_test2 <- data |>
      select(any_of(non_stationary_vars)) |>
      map(~ diff(.x) |>
        tseries::adf.test(k = 1) |>
        suppressWarnings())

    df2 <- data.frame(
      variable = character(),
      adf = numeric(),
      p_value = numeric()
    )

    for (i in seq_along(adf_test2)) {
      df2 <- rbind(df2, data.frame(
        variable = names(adf_test2)[i],
        adf = round(as.numeric(adf_test[[i]]$statistic), digits = 3),
        p_value = adf_test2[[i]]$p.value
      ))
    }

    df2$result <- case_when(
      df2$p_value >= 0.05 ~ "non-stationary",
      df2$p_value < 0.05 ~ "stationary"
    )

    return(list(
      "Level_Test" = df,
      "First_Difference_Test" = df2
    ))
  }

  return(df)
}

#' Estimate ARDL Model
#'
#' This function estimates an ARDL model with automatic lag selection based on
#' information criteria.
#'
#' @param data_list A list of dataframes containing the variables
#' @param max_lag Maximum number of lags to consider
#' @return A list of estimated ARDL models
#' @examples
#' # Example with a list of dataframes containing time series data
#' data_list <- list(df1 = data.frame(y = rnorm(100), x = rnorm(100)))
#' estimate_ardl(data_list, max_lag = 3)
#'
estimate_ardl <- function(data_list, max_lag = 3) {
  # Remove date column if present
  data <- data_list |> map(~ .x[, -1])

  # Get names of independent variables
  independent_vars <- colnames(data[[1]])[-1]

  # Create formula for ARDL estimation
  formula <- as.formula(paste0(
    colnames(data[[1]])[1],
    " ~ ",
    paste0(independent_vars, collapse = "+")
  ))

  # Estimate ARDL model for each dataset using BIC
  best_ardl_bic <- data |>
    purrr::map(
      ~ ARDL::auto_ardl(formula,
        data = .x,
        max_order = rep(max_lag, ncol(.x))
      )$best_order,
      selection = "BIC"
    )

  # Initialize list for ARDL models
  ardl_list <- list()

  # Estimate ARDL model for each dataset
  for (i in 1:length(data)) {
    order <- as.numeric(best_ardl_bic[[i]])

    ardl_list[[i]] <- ARDL::ardl(formula,
      data = data[[i]],
      order = order
    )
  }

  names(ardl_list) <- names(data)
  return(ardl_list)
}

#' Perform Bounds Test
#'
#' This function performs the bounds test for cointegration on ARDL models and
#' filters results based on the Breusch-Godfrey test.
#'
#' @param ardl List of ARDL models
#' @return List of bounds test results for models that pass diagnostics
#' @examples
#' # Assuming ardl_models is a list of estimated ARDL models
#' run_bound_test(ardl_models)
#'
run_bound_test <- function(ardl) {
  # Initialize list for bounds test results
  bound_test <- list()

  # Perform bounds test for each ARDL model
  for (i in 1:length(ardl)) {
    bound_test[[i]] <- ARDL::bounds_f_test(ardl[[i]], case = 2)
  }

  names(bound_test) <- names(ardl)

  # Perform Breusch-Godfrey test
  bg_test <- set_names(
    ardl |> map(~ lmtest::bgtest(.x)),
    names(ardl)
  )

  df_bg_test <- data.frame(id = character(), pvalue = numeric())

  # Filter models based on Breusch-Godfrey test
  for (i in seq_along(ardl)) {
    pvalue <- bg_test[[i]]$p.value

    if (pvalue >= 0.05) {
      df_bg_test <- rbind(
        df_bg_test,
        data.frame(
          id = names(ardl)[i],
          pvalue = pvalue
        )
      )
    }
  }

  id <- df_bg_test$id
  filtered_results <- bound_test[id]

  return(filtered_results)
}
