#!/usr/bin/env Rscript

# POM component shared by Supplementary Tables S2 and S3.
# Usage (normally called by the Python entry scripts):
#   Rscript code/white_wine_POM.R 7 /temporary/path/pom_repeat_scores.csv
#   Rscript code/white_wine_POM.R 5 /temporary/path/pom_repeat_scores.csv
# External R dependency: ordinal.

LOG_PROBABILITY_FLOOR <- 1e-12
PROBABILITY_TOLERANCE <- 1e-8
PREDICTORS <- c(
  "fixed_acidity",
  "volatile_acidity",
  "citric_acid",
  "residual_sugar",
  "chlorides",
  "free_sulfur_dioxide",
  "total_sulfur_dioxide",
  "density",
  "ph",
  "sulphates",
  "alcohol"
)

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

parse_arguments <- function() {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 2L) {
    stop(
      "Usage: Rscript code/white_wine_POM.R {5|7} OUTPUT_CSV",
      call. = FALSE
    )
  }
  n_classes <- suppressWarnings(as.integer(arguments[[1L]]))
  if (is.na(n_classes) || !n_classes %in% c(5L, 7L)) {
    stop("The number of categories must be 5 or 7", call. = FALSE)
  }
  list(n_classes = n_classes, output_path = arguments[[2L]])
}

load_inputs <- function(n_classes) {
  data_path <- file.path("data", "wine_white_processed.csv")
  split_path <- file.path("data", "white_repeated_stratified_5x5.csv")
  if (!file.exists(data_path) || !file.exists(split_path)) {
    stop(
      paste0(
        "Run from the repository root and place wine_white_processed.csv ",
        "and white_repeated_stratified_5x5.csv directly under data/."
      ),
      call. = FALSE
    )
  }

  data <- utils::read.csv(data_path, check.names = FALSE)
  splits <- utils::read.csv(split_path, check.names = FALSE)
  missing_data <- setdiff(c("quality", PREDICTORS), names(data))
  if (length(missing_data) > 0L) {
    stopf("The data file is missing: %s", paste(missing_data, collapse = ", "))
  }
  required_split <- c("row_index", "repeat_id", "fold_id", "split")
  missing_split <- setdiff(required_split, names(splits))
  if (length(missing_split) > 0L) {
    stopf("The split file is missing: %s", paste(missing_split, collapse = ", "))
  }

  if (!"row_index" %in% names(data)) {
    data$row_index <- seq_len(nrow(data)) - 1L
  }
  data$row_index <- as.integer(data$row_index)
  data$quality <- as.integer(data$quality)
  if (anyNA(data$row_index) || anyDuplicated(data$row_index)) {
    stop("row_index is missing or duplicated", call. = FALSE)
  }
  predictor_matrix <- as.matrix(data[PREDICTORS])
  storage.mode(predictor_matrix) <- "double"
  if (anyNA(predictor_matrix) || any(!is.finite(predictor_matrix))) {
    stop("A predictor is missing or non-finite", call. = FALSE)
  }

  if (n_classes == 7L) {
    response_map <- c(`3` = 1L, `4` = 2L, `5` = 3L, `6` = 4L,
                      `7` = 5L, `8` = 6L, `9` = 7L)
  } else {
    response_map <- c(`3` = 1L, `4` = 1L, `5` = 2L, `6` = 3L,
                      `7` = 4L, `8` = 5L, `9` = 5L)
  }
  data$y_order <- as.integer(unname(response_map[as.character(data$quality)]))
  data$y_index <- data$y_order - 1L
  if (anyNA(data$y_order) ||
      !identical(sort(unique(data$y_order)), seq_len(n_classes))) {
    stop("The response mapping failed", call. = FALSE)
  }

  splits$row_index <- as.integer(splits$row_index)
  splits$repeat_id <- as.integer(splits$repeat_id)
  splits$fold_id <- as.integer(splits$fold_id)
  splits$split <- as.character(splits$split)
  if (!identical(sort(unique(splits$repeat_id)), 1:5) ||
      !identical(sort(unique(splits$fold_id)), 1:5)) {
    stop("The split file must contain repeats 1,...,5 and folds 1,...,5", call. = FALSE)
  }
  if (!setequal(unique(splits$split), c("inner_fit", "validation", "test"))) {
    stop("Unexpected split labels", call. = FALSE)
  }
  for (repeat_id in 1:5) {
    tested <- splits$row_index[
      splits$repeat_id == repeat_id & splits$split == "test"
    ]
    if (length(tested) != nrow(data) || length(unique(tested)) != nrow(data)) {
      stopf("Repeat %d does not test every row exactly once", repeat_id)
    }
  }
  list(data = data, splits = splits)
}

fit_scaler <- function(data) {
  x <- as.matrix(data[PREDICTORS])
  storage.mode(x) <- "double"
  center <- colMeans(x)
  scale <- sqrt(colMeans(sweep(x, 2L, center, FUN = "-")^2))
  if (any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0)) {
    stop("The outer-training scaler is invalid", call. = FALSE)
  }
  list(center = center, scale = scale)
}

apply_scaler <- function(data, scaler) {
  x <- as.matrix(data[PREDICTORS])
  storage.mode(x) <- "double"
  x <- sweep(x, 2L, scaler$center, FUN = "-")
  x <- sweep(x, 2L, scaler$scale, FUN = "/")
  output <- as.data.frame(x, check.names = FALSE)
  names(output) <- PREDICTORS
  output
}

score_probabilities <- function(probabilities, y_index) {
  n <- nrow(probabilities)
  n_classes <- ncol(probabilities)
  observed_probability <- probabilities[cbind(seq_len(n), y_index + 1L)]
  log_score <- -log(pmax(observed_probability, LOG_PROBABILITY_FLOOR))
  predicted_cdf <- t(apply(probabilities, 1L, cumsum))[
    , seq_len(n_classes - 1L), drop = FALSE
  ]
  observed_cdf <- vapply(
    0:(n_classes - 2L),
    function(cutpoint) as.numeric(y_index <= cutpoint),
    numeric(n)
  )
  # Ordinary unnormalized RPS: no division by n_classes - 1.
  rps <- rowSums((predicted_cdf - observed_cdf)^2)
  if (any(!is.finite(log_score)) || any(!is.finite(rps))) {
    stop("A held-out score is non-finite", call. = FALSE)
  }
  list(LogS = log_score, RPS = rps)
}

run_fold <- function(data, splits, repeat_id, fold_id, n_classes) {
  assignment <- splits[
    splits$repeat_id == repeat_id & splits$fold_id == fold_id,
    c("row_index", "split"),
    drop = FALSE
  ]
  if (nrow(assignment) != nrow(data) || anyDuplicated(assignment$row_index)) {
    stop("A repeat/fold assignment must contain every row once", call. = FALSE)
  }
  fold_data <- merge(
    assignment, data, by = "row_index", all.x = TRUE, sort = FALSE
  )
  outer_training <- fold_data[fold_data$split != "test", , drop = FALSE]
  outer_test <- fold_data[fold_data$split == "test", , drop = FALSE]
  outer_training <- outer_training[order(outer_training$row_index), , drop = FALSE]
  outer_test <- outer_test[order(outer_test$row_index), , drop = FALSE]
  if (nrow(outer_training) == 0L || nrow(outer_test) == 0L) {
    stop("Outer training or test is empty", call. = FALSE)
  }
  if (!identical(sort(unique(outer_training$y_order)), seq_len(n_classes))) {
    stop("An outer training set is missing an outcome category", call. = FALSE)
  }

  scaler <- fit_scaler(outer_training)
  training_scaled <- apply_scaler(outer_training, scaler)
  test_scaled <- apply_scaler(outer_test, scaler)
  training_scaled$response <- ordered(
    outer_training$y_order,
    levels = seq_len(n_classes)
  )
  formula <- stats::reformulate(PREDICTORS, response = "response")
  fit <- ordinal::clm(
    formula = formula,
    data = training_scaled,
    link = "logit",
    threshold = "flexible",
    na.action = stats::na.fail,
    control = ordinal::clm.control(
      method = "Newton",
      maxIter = 100L,
      gradTol = 1e-6,
      convergence = "stop"
    )
  )
  convergence_code <- if (is.list(fit$convergence)) {
    suppressWarnings(as.integer(fit$convergence$code))
  } else {
    suppressWarnings(as.integer(fit$convergence))
  }
  if (length(convergence_code) != 1L || is.na(convergence_code) ||
      convergence_code != 0L) {
    stop("The POM did not converge regularly", call. = FALSE)
  }

  prediction <- stats::predict(
    fit,
    newdata = test_scaled[PREDICTORS],
    type = "prob"
  )
  probabilities <- as.matrix(prediction$fit)
  if (!identical(dim(probabilities), c(nrow(outer_test), n_classes))) {
    stop("The POM probability array has the wrong shape", call. = FALSE)
  }
  expected_names <- as.character(seq_len(n_classes))
  if (!is.null(colnames(probabilities))) {
    if (!setequal(colnames(probabilities), expected_names)) {
      stop("The POM probability columns have unexpected names", call. = FALSE)
    }
    probabilities <- probabilities[, expected_names, drop = FALSE]
  }
  if (any(!is.finite(probabilities)) ||
      min(probabilities) < -PROBABILITY_TOLERANCE ||
      max(abs(rowSums(probabilities) - 1)) > PROBABILITY_TOLERANCE) {
    stop("The POM probabilities failed validation", call. = FALSE)
  }
  scores <- score_probabilities(probabilities, outer_test$y_index)
  data.frame(
    row_index = outer_test$row_index,
    LogS = scores$LogS,
    RPS = scores$RPS
  )
}

main <- function() {
  arguments <- parse_arguments()
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop(
      "R package 'ordinal' is required. Install it with install.packages('ordinal').",
      call. = FALSE
    )
  }
  inputs <- load_inputs(arguments$n_classes)
  cat("R", as.character(getRversion()), "; ordinal",
      as.character(utils::packageVersion("ordinal")), "\n")
  cat("POM categories:", arguments$n_classes, "\n")
  cat("RPS: ordinary unnormalized sum over cumulative cutpoints\n")

  repeat_results <- vector("list", 5L)
  for (repeat_id in 1:5) {
    fold_results <- vector("list", 5L)
    for (fold_id in 1:5) {
      cat(sprintf("  POM repeat %d, fold %d\n", repeat_id, fold_id))
      fold_results[[fold_id]] <- run_fold(
        inputs$data,
        inputs$splits,
        repeat_id,
        fold_id,
        arguments$n_classes
      )
    }
    held_out <- do.call(rbind, fold_results)
    if (nrow(held_out) != nrow(inputs$data) ||
        anyDuplicated(held_out$row_index)) {
      stopf("POM repeat %d does not contain one score per row", repeat_id)
    }
    repeat_results[[repeat_id]] <- data.frame(
      repeat_id = repeat_id,
      LogS = mean(held_out$LogS),
      RPS = mean(held_out$RPS)
    )
  }
  repeat_scores <- do.call(rbind, repeat_results)
  output_directory <- dirname(arguments$output_path)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  utils::write.csv(repeat_scores, arguments$output_path, row.names = FALSE)
  print(repeat_scores, row.names = FALSE, digits = 8)
}

main()
