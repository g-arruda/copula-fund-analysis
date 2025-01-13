#' Perform Toda-Yamamoto Causality Test
#'
#' This function implements the Toda-Yamamoto approach to Granger causality testing.
#'
#' @param funds_list List of dataframes containing time series
#' @param ic Information criterion to use ("aic", "bic", or "hq")
#' @return List of causality test results
#' @examples
#' # Example with a list of dataframes containing time series data
#' data_list <- list(df1 = data.frame(y = rnorm(100), x = rnorm(100)))
#' toda_yamamoto(data_list, "aic")
#'
toda_yamamoto <- function(funds_list, ic) {
    require(tidyverse)
    require(vars)

    funds <- funds_list |> map(~ .x[, -1])
    result <- list()

    for (i in seq_along(funds)) {
        df_result <- data.frame(
            cause = character(),
            effect = character(),
            chi_square = numeric(),
            p_value = numeric()
        )

        for (j in 2:(ncol(funds[[i]]))) {
            var_selec <- funds[[i]][, c(1, j)] |>
                VARselect(lag.max = 15, type = "none")

            var_ic <- var_selec$selection[str_detect(
                names(var_selec$selection),
                toupper(ic)
            )] |>
                as.numeric()

            var_model <- funds[[i]][, c(1, j)] |>
                VAR(p = var_ic + 1, type = "none")

            # Wald test for causality
            test_result_forward <- wald.test(
                Sigma = vcov(var_model$varresult[[1]]),
                b = coef(var_model$varresult[[1]]),
                Terms = seq(from = 2, to = (var_model$p * 2), by = 2)
            )

            test_result_backward <- wald.test(
                Sigma = vcov(var_model$varresult[[2]]),
                b = coef(var_model$varresult[[2]]),
                Terms = seq(1, (var_model$p * 2), 2)
            )

            df_result <- rbind(
                df_result,
                data.frame(
                    cause = names(var_model$varresult)[[2]],
                    effect = names(var_model$varresult)[[1]],
                    chi_square = test_result_forward$result$chi2["chi2"] |>
                        as.numeric() |>
                        round(3),
                    p_value = test_result_forward$result$chi2["P"] |>
                        as.numeric() |>
                        round(3)
                ),
                data.frame(
                    cause = names(var_model$varresult)[[1]],
                    effect = names(var_model$varresult)[[2]],
                    chi_square = test_result_backward$result$chi2["chi2"] |>
                        as.numeric() |>
                        round(3),
                    p_value = test_result_backward$result$chi2["P"] |>
                        as.numeric() |>
                        round(3)
                )
            )
        }

        result[[i]] <- df_result
    }

    result <- set_names(result, names(funds))

    return(result)
}
