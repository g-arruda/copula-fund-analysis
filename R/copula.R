# Copula Analysis Functions
# Author: Gabriel de Almeida Arruda
# Description: A collection of functions for copula-based dependence analysis

library(tidyverse)
library(VineCopula)
library(fitdistrplus)
library(kdecopula)

#' Find Best Distribution for Each Variable
#'
#' Tests multiple distributions against each numeric column in the dataset
#' using Kolmogorov-Smirnov test and returns the best fit.
#'
#' @param df Dataframe containing the variables to analyze
#' @param distributions Vector of distribution names to test
#' @return Dataframe with best distribution for each variable
#' @examples
#' df <- data.frame(x = rnorm(100), y = rexp(100))
#' distributions <- c("norm", "exp", "t", "logis")
#' find_best_distribution(df, distributions)
#'
find_best_distribution <- function(df, distributions) {
  if (!is.data.frame(df)) {
    stop("df must be a dataframe")
  }
  if (length(df) == 0) {
    stop("df is empty")
  }
  if (!is.character(distributions) || length(distributions) == 0) {
    stop("distributions must be a non-empty character vector")
  }

  results_list <- list()

  for (colname in colnames(df)) {
    if (!is.numeric(df[[colname]])) {
      next
    }

    best_p_value <- 0
    best_stats <- list()

    for (dist in distributions) {
      if (colname == "Return" && dist == "cauchy") {
        next
      }

      tryCatch(
        {
          fit <- if (dist == "t") {
            fitdistrplus::fitdist(df[[colname]], dist, start = list(df = 5))
          } else {
            fitdistrplus::fitdist(df[[colname]], dist)
          }

          params <- fit$estimate
          ks_test <- if (length(params) == 1) {
            ks.test(df[[colname]], paste0("p", dist), params[[1]])
          } else if (length(params) == 2) {
            ks.test(df[[colname]], paste0("p", dist), params[[1]], params[[2]])
          }

          if (ks_test$p.value > best_p_value) {
            best_p_value <- ks_test$p.value
            best_stats <- list(
              variable = colname,
              distribution = dist,
              ks_statistic = ks_test$statistic,
              p_value = ks_test$p.value,
              param1 = params[[1]],
              param2 = if (length(params) > 1) params[[2]] else NA
            )
          }
        },
        error = function(e) {
          warning(paste("Error fitting", dist, "to", colname, ":", e$message))
        }
      )
    }

    if (length(best_stats) > 0) {
      results_list[[colname]] <- best_stats
    }
  }

  # Convert list to tibble
  results_df <- do.call(rbind, results_list) |>
    as.data.frame() |>
    tibble::as_tibble()

  return(results_df)
}

#' Create Cumulative Distribution Functions
#'
#' Creates CDFs for each variable based on their best-fit distributions.
#'
#' @param df Dataframe with original data
#' @param dist Dataframe with distribution information from find_best_distribution
#' @return Dataframe of CDF values
#' @examples
#' df <- data.frame(x = rnorm(100))
#' dist_info <- find_best_distribution(df, c("norm", "t"))
#' create_cdf(df, dist_info)
#'
create_cdf <- function(df, dist) {
  if (ncol(df) != nrow(dist)) {
    stop("Number of variables in df differs from number of rows in dist")
  }

  cdf_list <- vector("list", length = nrow(dist))

  for (i in seq_len(nrow(dist))) {
    cdf_name <- paste0("p", dist[i, 2])
    cdf <- get(cdf_name)

    serie <- as.numeric(df[[i]])
    param_1 <- as.numeric(dist[i, 5])
    param_2 <- as.numeric(dist[i, 6])

    if (is.na(param_1)) {
      stop(paste("First parameter for variable", i, "is not numeric"))
    }

    cdf_list[[i]] <- if (!is.na(param_2) && param_2 != 0) {
      cdf(serie, param_1, param_2)
    } else {
      cdf(serie, param_1)
    }
  }

  cdf_list <- purrr::set_names(cdf_list, dist[[1]])
  return(do.call(cbind, cdf_list) |> tibble::as_tibble())
}

#' Find Best Copula Family
#'
#' Identifies the best copula family for the data using specified criterion.
#'
#' @param data Dataframe of CDFs
#' @param B Number of bootstrap replications
#' @param criterion Selection criterion ("logLik", "AIC", or "BIC")
#' @return Dataframe with best copula family for each pair
#' @examples
#' data <- data.frame(x = runif(100), y = runif(100))
#' find_best_family(data, 100, "AIC")
#'
find_best_family <- function(data, B, criterion) {
  set.seed(999)
  copula_families <- c(1, 2, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)
  num_col <- ncol(data) - 1
  results <- vector(mode = "list", length = num_col)
  u <- VineCopula::pobs(data)

  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (!is.numeric(B) || B <= 0 || B != round(B)) {
    stop("B must be a positive integer")
  }
  if (!criterion %in% c("logLik", "AIC", "BIC")) {
    stop("criterion must be 'logLik', 'AIC' or 'BIC'")
  }

  for (i in seq_len(num_col)) {
    copula_selection <- VineCopula::BiCopSelect(
      u[, 1], u[, i + 1],
      familyset = copula_families,
      rotations = FALSE,
      selectioncrit = criterion
    )

    param1 <- copula_selection$par
    param2 <- ifelse(copula_selection$par2 == 0, 0, copula_selection$par2)

    gof_test <- VineCopula::BiCopGofTest(
      u[, 1], u[, i + 1],
      family = copula_selection$family,
      par = param1,
      par2 = param2,
      method = ifelse(copula_selection$family == 2, "white", "kendall"),
      B = B
    )

    max_pvalue <- max(gof_test$p.value, gof_test$p.value.CvM)

    selec_number <- case_when(
      criterion == "logLik" ~ 11,
      criterion == "AIC" ~ 12,
      criterion == "BIC" ~ 13
    )

    result <- tibble::tibble(
      var = colnames(u)[i + 1],
      family_number = copula_selection$family,
      family_name = copula_selection$familyname,
      value = copula_selection[[selec_number]],
      gof_pvalue = max_pvalue
    ) |>
      dplyr::rename(!!criterion := value)

    results[[i]] <- result
    cat("\rFinalizada a copula", i, "de", num_col, "\n")
  }

  return(do.call(rbind, results))
}

#' Estimate Bivariate Copula
#'
#' Estimates bivariate copula parameters for given data and family.
#'
#' @param data CDF values
#' @param family Copula family information from find_best_family
#' @return List of estimated copula objects
#' @examples
#' data <- data.frame(x = runif(100), y = runif(100))
#' family <- find_best_family(data, 100, "AIC")
#' estimate_bivariate_copula(data, family)
#'
estimate_bivariate_copula <- function(data, family) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (!is.data.frame(family)) {
    stop("family must be a data frame")
  }

  num_col <- ncol(data) - 1
  result <- setNames(vector(mode = "list", length = num_col), names(data)[-1])

  for (i in seq_len(num_col)) {
    family_index <- family$family_number[i]
    if (!is.numeric(family_index)) {
      stop("Family indices must be numeric")
    }

    estimate_param <- VineCopula::BiCopEst(
      data[[1]], data[[i + 1]],
      family = family_index,
      max.df = 200
    )

    copula <- VineCopula::BiCop(
      family = family_index,
      par = estimate_param$par,
      par2 = estimate_param$par2
    )

    result[[i]] <- copula
  }

  return(result)
}

#' Extract Copula Parameters
#'
#' Extracts and formats parameters from estimated copulas.
#'
#' @param data List of BiCop objects
#' @return Tibble with copula parameters
#' @examples
#' # Assuming copulas is a list of estimated BiCop objects
#' get_copula_parameters(copulas)
#'
get_copula_parameters <- function(data) {
  if (!is.list(data)) {
    stop("data must be a list of BiCop objects")
  }

  result <- list()
  num_copulas <- length(data)

  for (i in seq_len(num_copulas)) {
    copula <- data[[i]]
    if (!inherits(copula, "BiCop")) {
      stop("All elements must be BiCop objects")
    }

    result[[i]] <- tibble::tibble(
      variable = names(data[i]),
      param_1 = copula[[2]],
      param_2 = copula[[3]],
      tau = copula[[6]]
    ) |>
      dplyr::mutate(
        dplyr::across(
          dplyr::where(is.numeric),
          ~ scales::number(., accuracy = 0.0001, big.mark = ".", decimal.mark = ",")
        )
      )
  }

  return(do.call(rbind, result))
}

#' Simulate Data from Bivariate Copula
#'
#' Simulates data from estimated copulas and transforms back to original scale.
#'
#' @param copulas List of estimated copulas
#' @param orig_data Original data
#' @param best_dist Best distribution information
#' @param simulations Number of simulations
#' @param seed Random seed
#' @return Dataframe of simulated values
#' @examples
#' # Assuming we have estimated copulas, original data, and distribution info
#' simulate_bivariate_copula(copulas, orig_data, dist_info, 1000, 123)
#'
simulate_bivariate_copula <- function(copulas, orig_data, best_dist,
                                      simulations, seed) {
  set.seed(seed)
  df_list <- list()

  for (j in seq_along(copulas)) {
    set.seed(seed)
    sim_test <- VineCopula::BiCopSim(N = simulations, obj = copulas[[j]]) |>
      as.data.frame()

    names(sim_test) <- c(names(orig_data)[1], names(orig_data)[j + 1])
    list_transform <- list()

    for (i in seq_len(ncol(sim_test))) {
      distr_function <- get(paste0("q", best_dist[i, 2]))
      list_transform[[i]] <- distr_function(
        sim_test[[i]],
        best_dist[[5]][i],
        best_dist[[6]][i]
      )
    }

    df_transform <- purrr::set_names(list_transform, names(sim_test)) |>
      dplyr::bind_cols()

    df_list[[j]] <- df_transform
  }

  result_list <- list()
  for (i in seq_along(df_list)) {
    result_list[[i]] <- if (i == 1) {
      df_list[[i]]
    } else {
      df_list[[i]][, 2]
    }
  }

  return(dplyr::bind_cols(result_list))
}



#' Check Copula Model Fit
#'
#' Compares simulated and observed data to assess the fit of copula models.
#' Calculates various statistics to evaluate how well the model captures the
#' original data's characteristics.
#'
#' @param copulas List of estimated copulas
#' @param orig_data Original data
#' @param best_dist Best distribution information
#' @param simulations Number of simulations
#' @param seed Random seed for reproducibility
#' @return Tibble with comparison statistics
#' @examples
#' # Assuming we have estimated copulas, original data, and distribution info
#' check_copula_fit(copulas, orig_data, best_dist, 1000, 123)
#'
check_copula_fit <- function(copulas, orig_data, best_dist, simulations, seed) {
  set.seed(seed)

  list_transform <- list()
  list_analysis <- list()
  list_result <- list()
  result_list <- list()

  for (j in seq_along(copulas)) {
    set.seed(seed)
    # Simulate from copula
    sim_test <- VineCopula::BiCopSim(N = simulations, obj = copulas[[j]]) |>
      as.data.frame()

    names(sim_test) <- c(names(orig_data)[1], names(orig_data)[j + 1])

    # Transform simulated data back to original scale
    for (i in seq_len(ncol(sim_test))) {
      distr_function <- get(paste0("q", best_dist[i, 2]))
      list_transform[[i]] <- distr_function(
        sim_test[[i]],
        best_dist[[5]][i],
        best_dist[[6]][i]
      )
    }

    # Combine transformed data
    df_transform <- purrr::set_names(list_transform, names(sim_test)) |>
      dplyr::bind_cols()

    # Calculate comparison statistics
    for (i in seq_len(ncol(df_transform))) {
      list_analysis[[i]] <- tibble::tibble(
        name = names(sim_test)[i],
        max_original = orig_data[[names(sim_test)[i]]] |> max(),
        max_sim = df_transform[[i]] |> max(),
        mean_original = orig_data[[names(sim_test)[i]]] |> mean(),
        mean_sim = df_transform[[i]] |> mean(),
        min_original = orig_data[[names(sim_test)[i]]] |> min(),
        min_sim = df_transform[[i]] |> min(),
        median_original = orig_data[[names(sim_test)[i]]] |> median(),
        median_sim = df_transform[[i]] |> median(),
        skewness_original = moments::skewness(orig_data[[names(sim_test)[i]]]),
        skewness_sim = moments::skewness(df_transform[[i]]),
        kurtosis_original = moments::kurtosis(orig_data[[names(sim_test)[i]]]),
        kurtosis_sim = moments::kurtosis(df_transform[[i]])
      )
    }

    # Format results
    list_result[[j]] <- dplyr::bind_rows(list_analysis) |>
      tidyr::pivot_longer(cols = -name, names_to = "statistic", values_to = "value") |>
      tidyr::pivot_wider(names_from = name, values_from = value)
  }

  # Combine results
  for (i in seq_along(list_result)) {
    result_list[[i]] <- if (i == 1) {
      list_result[[i]]
    } else {
      list_result[[i]][, 3]
    }
  }

  return(dplyr::bind_cols(result_list))
}

#' Calculate Conditional Probabilities
#'
#' Calculates probabilities of certain events based on simulated data from copula models.
#' Specifically, it calculates probabilities of returns being positive or negative
#' given that a variable is within certain intervals.
#'
#' @param df_sim List of dataframes with simulated data
#' @return List of dataframes with probability calculations
#' @examples
#' # Assuming df_sim is a list of dataframes with simulated data
#' calculate_probabilities(df_sim)
#'
calculate_probabilities <- function(df_sim) {
  list_result <- list()

  for (j in seq_along(df_sim)) {
    # Calculate intervals based on quantiles and mean
    interval <- c(
      quantile(df_sim[[j]][[2]]),
      mean(df_sim[[j]][[2]])
    ) |>
      sort()

    # Find median index
    median_index <- which(interval == quantile(df_sim[[j]][[2]], 0.5))
    interval <- interval[-median_index]

    # Initialize results dataframe
    results <- data.frame(
      lower = numeric(),
      upper = numeric(),
      probability = numeric(),
      mean_return = numeric(),
      return_type = character()
    )

    # Calculate probabilities for each interval
    for (i in 1:(length(interval) - 1)) {
      # Find indices within current interval
      indices <- which(
        df_sim[[j]][[2]] >= interval[i] &
          df_sim[[j]][[2]] < interval[i + 1]
      )

      # Calculate probabilities and mean returns
      prob_negative <- sum(df_sim[[j]][[1]][indices] < 0) /
        length(df_sim[[j]][[1]])
      prob_positive <- sum(df_sim[[j]][[1]][indices] > 0) /
        length(df_sim[[j]][[1]])

      # Calculate mean return for the interval
      mean_return <- if (length(indices) > 0) {
        mean(df_sim[[j]][[1]][indices])
      } else {
        0
      }

      # Add results to dataframe
      results <- rbind(
        results,
        data.frame(
          lower = interval[i] * 100,
          upper = interval[i + 1] * 100,
          probability = prob_negative * 100,
          mean_return = mean_return * 100,
          return_type = "negative"
        ) |>
          dplyr::mutate(
            dplyr::across(
              dplyr::where(is.numeric),
              ~ scales::number(., accuracy = 0.001, big.mark = ".", decimal.mark = ",")
            )
          ),
        data.frame(
          lower = interval[i] * 100,
          upper = interval[i + 1] * 100,
          probability = prob_positive * 100,
          mean_return = mean_return * 100,
          return_type = "positive"
        ) |>
          dplyr::mutate(
            dplyr::across(
              dplyr::where(is.numeric),
              ~ scales::number(., accuracy = 0.001, big.mark = ".", decimal.mark = ",")
            )
          )
      )
    }

    # Order results by return type
    results <- results[order(results$return_type, decreasing = TRUE), ]

    # Add descriptive row names
    rownames(results) <- c(
      "minimum to 1st quartile",
      "1st quartile to mean",
      "mean to 3rd quartile",
      "3rd quartile to maximum",
      "minimum to 1st quartile",
      "1st quartile to mean",
      "mean to 3rd quartile",
      "3rd quartile to maximum"
    )

    list_result[[j]] <- results
  }

  # Name the list elements
  names(list_result) <- names(df_sim)

  return(list_result)
}
