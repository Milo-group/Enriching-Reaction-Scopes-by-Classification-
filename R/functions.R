#' @title Package Load
#' @description Sets default options upon package loading.
#' @export
.onLoad <- function(libname, pkgname) {
  default::default(data.frame) <- list(check.names = FALSE)
  glm.control(maxit = 100)
}

#' Stratified Sampling Function
#'
#' This function returns a stratified sample from a dataframe `df` based on a grouping variable `group`
#' and the desired sample size `size`. The function accommodates different sample size scenarios,
#' including varying sample sizes per group, proportions, and handling of groups with fewer observations
#' than the requested sample size.
#'
#' @param df A dataframe to sample from.
#' @param group The grouping variable(s) by which the stratified sampling will be performed.
#' @param size Desired sample size, which can be:
#'             - A fixed number (>=1) of samples per group.
#'             - A proportion (<1) to sample from each group.
#'             - A vector of specific sample sizes for each group.
#' @return A dataframe containing the stratified sample.
#'
#' @details
#' The function handles different sampling scenarios:
#' - If `size` is a vector, it must match the number of groups or be named to match group levels.
#' - If `size` is less than 1, it is treated as a proportion of each group's size.
#' - If `size` is a fixed number, the function checks if each group has enough samples; if not, it returns all samples from smaller groups.
#'
#' @examples
#' # Example 1: Stratified sample with 10 samples per group
#' # result <- stratified(df = my_data, group = "my_group", size = 10)
#'
#' # Example 2: Stratified sample where 20% of each group is sampled
#' # result <- stratified(df = my_data, group = "my_group", size = 0.2)
#'
#' # Example 3: Stratified sample with varying sizes per group
#' # result <- stratified(df = my_data, group = "my_group", size = c(A = 5, B = 10, C = 3))
#' @export

stratified <- function(df, group, size) {
  # Create a factor by interacting the grouping variable(s) and dropping unused levels.
  df.interaction <- interaction(df[group], drop = TRUE)

  # Create a table of the counts of each group level.
  df.table <- table(df.interaction)

  # Split the dataframe into a list, where each element is a subset of the dataframe for a group level.
  df.split <- split(df, df.interaction)

  # Handle case where `size` is a vector of sizes for each group.
  if (length(size) > 1) {
    # Ensure the number of sizes matches the number of group levels.
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))

    # If the size vector is unnamed, assign names based on group levels.
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
    } else {
      # Ensure the named size vector matches the group levels.
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)], # Reorder size vector to match group levels.
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  }
  # Handle case where `size` is a proportion (less than 1) to sample from each group.
  else if (size < 1) {
    n <- round(df.table * size, digits = 0)  # Calculate the number of samples as a proportion of the group size.
  }
  # Handle case where `size` is a fixed number of samples per group.
  else if (size >= 1) {
    # Ensure all groups have enough samples or replace samples if necessary.
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))  # Assign the fixed size to each group.
    } else {
      # Provide a message if some groups have fewer samples than the desired size.
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")

      # Return all samples from groups that have fewer observations than the requested size.
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }

  # For each group, sample `n[x]` observations without replacement.
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = F), ])

  # Combine all the sampled subsets into a single dataframe.
  set1 <- do.call("rbind", temp)

  return(set1)  # Return the stratified sample.
}


#' K-Fold Logistic Cross-Validation
#'
#' This function performs K-fold cross-validation (CV) for ordinal or multinomial logistic regression
#' models. It allows stratification of data, calculates performance metrics, and returns accuracy,
#' the confusion matrix, and other relevant statistics.
#'
#' @param formula A formula specifying the model structure.
#' @param data A dataframe containing the data to fit the model.
#' @param ordinal A boolean indicating whether to fit an ordinal logistic regression (`TRUE`)
#'                or a multinomial logistic regression (`FALSE`). Defaults to `TRUE`.
#' @param folds Number of folds for cross-validation. If `NULL`, stratified sampling can be used.
#' @param outcome.column The column number of the outcome variable. Defaults to the "class" column.
#' @param stratify A boolean indicating whether to stratify the data for sampling. Defaults to `FALSE`.
#' @param sample.vector A vector specifying the number of samples to be drawn from each class during stratified sampling.
#'                      Defaults to the ratio of the smallest class size.
#'
#' @return A list containing:
#'         - `acc`: Accuracy of the model across folds.
#'         - `J.small`: J statistic (Youden's Index) for the smallest class.
#'         - `ct`: Confusion matrix.
#'         - `probs`: Dataframe of predicted probabilities and class predictions.
#'
#' @details
#' - If `folds` is specified, the function performs K-fold cross-validation.
#' - If `stratify` is `TRUE`, the data is stratified based on `sample.vector` before performing CV.
#' - Models can either be ordinal (`MASS::polr`) or multinomial (`nnet::multinom`).
#' - Performance metrics like accuracy and J statistic (for each class) are calculated.
#' @export
k.fold.log.cv <- function(formula, data, ordinal = T,
                          folds = NULL,
                          outcome.column = which(colnames(data) == 'class'),
                          stratify = F,
                          sample.vector = floor(round(summary(data$class)/min(summary(data$class))))) {
  # Initialize lists to store models, probabilities, accuracy, and class predictions.
  models <- list()
  probalities <- list()
  acc <- list()
  class.pred <- list()

  # Stratify the data if requested and no folds are provided.
  if (stratify == TRUE && !is.null(sample.vector) && is.null(folds)) {
    sets <- list()  # List to store stratified sets.

    # Create stratified folds based on sample.vector.
    for (i in 1:ceiling(dim(data)[1] / sum(sample.vector))) {
      # Create a pool of data to sample from.
      if (length(sets) == 0) {
        pool.data <- data
      } else {
        pool.data <- dplyr::setdiff(data, do.call(rbind, sets))
      }

      # If the pool contains fewer observations than the desired sample size, add them all to the set.
      if (dim(pool.data)[1] < sum(sample.vector) & dim(pool.data)[1] > 0) {
        sets[[i]] <- pool.data
        pool.data <- data.frame()  # Clear pool after adding remaining data.
      }

      # Handle cases where there are fewer observations than desired in any group.
      if (nrow(pool.data) > 0 & any(summary(pool.data$class) < sample.vector)) {
        for (j in 1:nrow(pool.data)) {
          # Distribute remaining data across sets.
          ifelse(j <= length(sets),
                 sets[[j]] <- rbind(sets[[j]], pool.data[j, ]),
                 sets[[(j - length(sets))]] <- rbind(sets[[(j - length(sets))]], pool.data[j, ]))
        }
        is_empty <- function(x) (nrow(x) == 0 || is.null(x))
        sets <- sets[sapply(sets, is_empty) == FALSE]  # Remove empty sets.
      } else {
        # Stratify the pool data if there are enough observations.
        if (nrow(pool.data) > 0) suppressWarnings(sets[[i]] <- stratified(pool.data, 'class', sample.vector))
      }

      # Check if all data has been stratified as desired.
      if (i == ceiling(dim(data)[1] / sum(sample.vector))) {
        check <- dim(do.call(rbind, sets)) == dim(data)
        if (!all(check)) {
          print("Wasn't able to stratify as wanted. Check with another sample vector.")
        }
      }
    }

    # Assign each set to a split and arrange by a custom flag.
    for (i in 1:length(sets)) {
      sets[[i]] <- dplyr::mutate(sets[[i]], split.assign = i)
    }
    new_dat <- do.call(rbind, sets)
    new_dat <- dplyr::arrange(new_dat, flag)

    # Train models and compute probabilities and class predictions for each split.
    for (i in 1:length(sort(unique(new_dat$split.assign)))) {
      train <- new_dat[new_dat$split.assign != i,]
      test <- new_dat[new_dat$split.assign == i,]

      # Define starting values for the ordinal logistic regression.
      num.of.vars <- stringi::stri_count(formula, fixed = '+')
      start <- c(rep(0, (num.of.vars + 2)), 1)

      # Fit the ordinal logistic regression model.
      models[[match(i,sort(unique(new_dat$split.assign)))]] <- if (ordinal == T) {
        MASS::polr(formula, data = train, Hess = T, start = start, control = list(maxit = 1000))
      } else {
        nnet::multinom(formula,data = train, maxit = 2000, trace = FALSE)
      }
      # Predict probabilities and classes.
      probalities[[match(i,sort(unique(new_dat$split.assign)))]] <- data.frame(predict(models[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                                       newdata = test, "probs")*100)
      class.pred[[match(i,sort(unique(new_dat$split.assign)))]] <- data.frame(predict(models[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                                      newdata = test, "class"))

      # Combine probabilities, predictions, and flags.
      probalities[[match(i,sort(unique(new_dat$split.assign)))]] <- cbind(probalities[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                          class.pred[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                          test$flag)
      names(probalities[[match(i,sort(unique(new_dat$split.assign)))]])[dim(probalities[[1]])[2] - 1] <- 'prediction'
      names(probalities[[match(i,sort(unique(new_dat$split.assign)))]])[dim(probalities[[1]])[2]] <- 'flag'
    }
  } else {
    # Handle case when folds are specified or folds equal the number of rows (Leave-One-Out Cross Validation).
    if (folds == nrow(data)) {
      split.assign <- sample(1:folds, nrow(data), replace = FALSE)  # LOO case.
    } else {
      split.assign <- caret::createFolds(1:dim(data)[1], folds, list = FALSE)  # K-fold CV.
    }
    new_dat <- cbind(data, split.assign)

    # Train models and compute probabilities for each fold.
    for (i in 1:folds) {
      train <- new_dat[split.assign != i,]
      test <- new_dat[split.assign == i,]

      num.of.vars <- stringi::stri_count(formula, fixed = '+')
      start <- c(rep(0, (num.of.vars + 2)), 1)

      models[[match(i,1:folds)]] <- if (ordinal == T) {
        MASS::polr(formula, data = train, Hess = T, start = start, control = list(maxit = 1000))
      } else {
        nnet::multinom(formula, data = train, maxit = 2000, trace = FALSE)
      }

      # Predict probabilities and classes.
      if (folds == nrow(data)) {
        probalities[[match(i,1:folds)]] <- data.table::transpose(data.frame(predict(models[[match(i,1:folds)]],
                                                                                    newdata = test, "probs")*100))
      } else {
        probalities[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]],
                                                              newdata = test, "probs")*100)
      }

      class.pred[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], newdata = test, "class"))

      # Combine probabilities, predictions, and flags.
      probalities[[match(i,1:folds)]] <- cbind(probalities[[match(i,1:folds)]],
                                               class.pred[[match(i,1:folds)]], test$flag)
      names(probalities[[match(i,1:folds)]])[dim(probalities[[1]])[2] - 1] <- 'prediction'
      names(probalities[[match(i,1:folds)]])[dim(probalities[[1]])[2]] <- 'flag'
    }
  }

  # Aggregate probabilities, predictions, and flags across all folds.
  probs <- data.frame(do.call(rbind, probalities))
  probs <- probs[order(probs$flag),]
  probs[,1:(dim(probalities[[1]])[2] - 2)] <- round(probs[,1:(dim(probalities[[1]])[2] - 2)], digits = 0)

  # Extract predictions and actual outcomes.
  pred <- probs[,(dim(probalities[[1]])[2] - 1)]
  actual <- data[[outcome.column]]

  # Calculate confusion table and accuracy.
  ct <- table(actual, pred)
  acc <- round((sum(diag(ct)) / sum(ct)) * 100, 2)

  # Calculate performance metrics (J statistic) for each class.
  ct.df <- data.frame(ct)
  TP <- ct.df$Freq[ct.df$actual == ct.df$pred]
  TP_FN <- list()
  TN <- list()
  TN_FP <- list()
  J <- list()

  for (i in levels(ct.df$actual)) {
    TP_FN[[i]] <- sum(ct.df$Freq[ct.df$actual == i])
    TP_FN[[i]] <- TP[which(levels(ct.df$actual) == i)] / TP_FN[[i]]
    TN[[i]] <- sum(ct.df$Freq[dplyr::intersect(which(ct.df$pred != i), which(ct.df$actual != i))])
    TN_FP[[i]] <- sum(ct.df$Freq[(ct.df$actual != i)])
    TN_FP[[i]] <- TN[[i]] / TN_FP[[i]]
    J[[i]] <- TP_FN[[i]] + TN_FP[[i]] - 1
  }

  # Return the overall accuracy and the J statistic for the smallest class.
  J.small <- J[which.min(summary(data$class))]
  return(list(acc, J.small, ct, probs))
}

#' Iterative K-Fold Logistic Cross-Validation
#'
#' This function performs multiple iterations of K-fold cross-validation (CV) for ordinal or multinomial
#' logistic regression models, allowing for repeated assessment of accuracy across multiple random splits.
#' The results include overall accuracy, best/worst iteration accuracy, and confusion tables for both best
#' and worst iterations.
#'
#' @param formula A formula specifying the model structure.
#' @param data A dataframe containing the data to fit the model.
#' @param ordinal A boolean indicating whether to fit an ordinal logistic regression (`TRUE`) or a multinomial
#'                logistic regression (`FALSE`). Defaults to `TRUE`.
#' @param folds Number of folds for cross-validation. If `NULL`, stratified sampling can be used.
#' @param out.col The column number of the outcome variable. Defaults to the "class" column.
#' @param stratify A boolean indicating whether to stratify the data for sampling. Defaults to `FALSE`.
#' @param sample.vector A vector specifying the number of samples to be drawn from each class during stratified sampling.
#'                      Defaults to the ratio of the smallest class size.
#' @param iterations Number of iterations to repeat the cross-validation. Each iteration results in a different random split.
#' @param verbose A boolean to control whether to print the classification tables and accuracy results.
#'                Defaults to `FALSE`.
#'
#' @return A list containing:
#'         - `over.all.accuracy`: Overall average accuracy across all iterations.
#'         - `best`: The best accuracy obtained across iterations.
#'         - `worst`: The worst accuracy obtained across iterations.
#'         - `cts`: A formatted table showing the best and worst confusion matrices.
#'
#' @details
#' - This function repeatedly calls the `k.fold.log.cv` function to perform cross-validation for the specified
#'   number of iterations.
#' - The function computes and prints the average accuracy, best and worst accuracies, and corresponding
#'   confusion matrices if `verbose` is set to `TRUE`.
#' @export
k.fold.log.iter <- function(formula, data, ordinal = T, folds = NULL,
                            out.col = which(colnames(data) == 'class'),
                            stratify = FALSE,
                            sample.vector = floor(round(summary(data$class) / min(summary(data$class)))),
                            iterations, verbose = FALSE) {

  # Initialize lists to store accuracy and classification table results
  iter.list <- list()  # Stores the accuracy for each iteration
  ct.list <- list()    # Stores the classification tables for each iteration

  # Loop over the specified number of iterations
  for (i in 1:iterations) {

    # Call the k-fold ordinal logistic regression function for each iteration
    mod <- k.fold.log.cv(formula, data, ordinal, folds, out.col, stratify, sample.vector)

    # Store the accuracy result from mod[[1]] and the classification table from mod[[3]]
    iter.list[[match(i, 1:iterations)]] <- mod[[1]]
    ct.list[[match(i, 1:iterations)]] <- mod[[3]]
  }

  # Calculate overall average accuracy across all iterations
  over.all.accuracy <- round(Reduce(`+`, iter.list) / iterations, digits = 2)

  # Identify the best and worst accuracy
  best <- iter.list[which.max(iter.list)]  # Best accuracy
  worst <- iter.list[which.min(iter.list)] # Worst accuracy

  # Create a table showing overall accuracy, best accuracy, and worst accuracy
  Accuracies <- knitr::kable(cbind(over.all.accuracy, best, worst))

  # Create a combined classification table for the best and worst accuracies
  tab <- suppressMessages(cbind(
    reshape2::dcast(data.frame(ct.list[which.max(iter.list)]), actual ~ pred),  # Best classification table
    rep("***", length(unique(data[, out.col]))),                                # Separator
    reshape2::dcast(data.frame(ct.list[which.min(iter.list)]), actual ~ pred)   # Worst classification table
  ))

  # Rename the column separating the two classification tables
  names(tab)[length(unique(data[, out.col])) + 2] <- ''

  # Format the combined classification table as a knitr table
  cts <- knitr::kable(tab)

  # If verbose is TRUE, print classification tables and accuracy results
  if (verbose == TRUE) {
    print(cts)
    print(Accuracies)
  }

  # Return the overall accuracy (invisibly)
  invisible(over.all.accuracy)
}

#' Submodel Logistic Regression (Ordinal/Multinomial)
#'
#' This function performs logistic regression by generating and fitting various combinations of predictor
#' variables. It can handle both ordinal logistic regression (via `MASS::polr`) and multinomial logistic
#' regression (via `nnet::multinom`). The models are evaluated using McFadden's R-squared value to rank
#' model performance.
#'
#' @param data A dataframe containing the data.
#' @param out.col The column number of the outcome variable. Defaults to the "class" column.
#' @param min The minimum number of predictor variables to include in combinations. Defaults to 1.
#' @param max The maximum number of predictor variables to include in combinations. Defaults to one-fifth of
#'            the number of rows in the dataset.
#' @param ordinal A boolean indicating whether to fit an ordinal logistic regression (`TRUE`) or a multinomial
#'                logistic regression (`FALSE`). Defaults to `TRUE`.
#'
#' @return A dataframe containing the top 10 models sorted by McFadden's R-squared value. The output includes:
#'         - `formula`: The model formula.
#'         - `McFadden R2`: McFadden's R-squared value, indicating model fit.
#'
#' @details
#' - This function generates combinations of predictor variables and fits logistic regression models for
#'   each combination.
#' - McFadden's R-squared value is used to evaluate model fit, with higher values indicating better fit.
#' - For ordinal regression, the `MASS::polr` function is used, while multinomial regression uses
#'   `nnet::multinom`.
#' - The function handles special characters in variable names by wrapping them in backticks.
#' @export
sub_model_log <- function(data, out.col = which(colnames(data) == 'class'),
                          min = 1, max = floor(dim(data)[1] / 5),
                          ordinal = T) {

  # Prepare the outcome variable name
  output <- stringr::str_c('`', names(data[out.col]), '`')

  # Get predictor variable names, excluding 'flag' (if present)
  vars <- names(data[, -out.col])
  vars <- vars[vars != 'flag']

  # Wrap predictor variable names with backticks to handle special characters
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c('`', vars[i], '`')
  }

  # Initialize lists to store combinations of variables and R-squared values
  comb.list <- list()
  R2.list <- list()

  # Generate combinations of predictors from 'min' to 'max' variables
  for (i in min:max) {
    comb.list[[i]] <- data.frame(aperm(combn(vars, i)), stringsAsFactors = FALSE)

    # Concatenate variable names to create model formulas
    comb.list[[i]][, dim(comb.list[[i]])[2] + 1] <- do.call(paste,
                                                            c(comb.list[[i]][names(comb.list[[i]])], sep = " + "))

    # Name the new column as 'formula'
    names(comb.list[[i]])[dim(comb.list[[i]])[2]] <- 'formula'

    # Remove individual variable columns, keeping only the formula
    for (co in names(comb.list[[i]])[1:length(names(comb.list[[i]])) - 1]) {
      comb.list[[i]][co] <- NULL
    }
  }

  # Remove empty combinations (if any) and combine them into a single data frame
  comb.list <- plyr::compact(comb.list)
  forms <- do.call(rbind, comb.list)
  names(forms) <- 'formula'

  # Define a helper function to fit an ordinal logistic regression model
  fit_polr_model <- function(formula, data, start_vector) {
    tryCatch({
      # Try fitting the model using MASS::polr
      model <- MASS::polr(formula, data = data, Hess = TRUE, start = start_vector, control = list(maxit = 100))
      return(list(success = TRUE, model = model))
    }, error = function(e) {
      return(list(success = FALSE, error = e))
    })
  }

  # Fit models for each formula
  if (ordinal == T) {
    for (i in 1:dim(forms)[1]) {

      # Construct the full formula with the outcome variable
      forms$formula[i] <- stringr::str_c(output, ' ~ ', forms$formula[i])

      # Count the number of variables in the formula (for initializing start vector)
      num.of.vars <- stringi::stri_count(forms$formula[i], fixed = '+')

      # Initialize start vector for model fitting (coefficients)
      start <- c(rep(0, num.of.vars + 2), 1)

      # Flag to control the loop in case the initial model fitting fails
      success <- FALSE

      # Keep trying to fit the model until it succeeds
      while (!success) {
        result <- fit_polr_model(forms$formula[i], data, start)

        if (result$success) {
          success <- TRUE
          test.1 <- result$model
        } else {
          # If model fitting fails, extend the start vector and retry
          start <- c(0, start)
        }
      }

      # Fit a null model (class ~ 1) for comparison
      test.0 <- MASS::polr(class ~ 1, data = data, Hess = TRUE)

      # Calculate log-likelihoods for the null and fitted models
      loglik.0 <- as.numeric(nnet:::logLik.multinom(test.0))
      loglik.1 <- as.numeric(nnet:::logLik.multinom(test.1))

      # Calculate McFadden's R-squared value and store it
      R2.list[i] <- round((1 - (loglik.1 / loglik.0)), digits = 3)
    }
  } else {
    for (i in 1:dim(forms)[1]) {
      # Construct the full formula with the outcome variable
      forms$formula[i] <- stringr::str_c(output,' ~ ',forms$formula[i])

      test.1 <- nnet::multinom(forms$formula[i],
                               data = data,
                               maxit = 2000, trace = F)
      test.0 <- nnet::multinom(class ~ 1, data = data, maxit = 2000, trace = F)
      loglik.0 <- as.numeric(nnet:::logLik.multinom(test.0))
      loglik.1 <- as.numeric(nnet:::logLik.multinom(test.1))
      R2.list[i] <- round((1 - (loglik.1/loglik.0)),digits = 3)
    }
  }


  # Add McFadden's R2 values to the forms data frame
  forms[, 2] <- do.call(rbind, R2.list)
  names(forms)[2] <- 'McFadden R2'

  # Return the top 10 models sorted by McFadden's R2
  out.models <- head(dplyr::arrange(forms, desc(forms[, 2])), 10)
  return(out.models)
}


#' Model Information and Performance Evaluation
#'
#' This function computes and displays key performance metrics for a trained classification model,
#' including accuracy and McFadden's pseudo R-squared. It also generates a confusion matrix and
#' displays the model coefficients.
#'
#' @param model A trained classification model. It could be an object of class `nnet::multinom`
#'   or `MASS::polr`.
#' @param data A data frame containing the predictors and the true class labels in a column
#'   named `class`.
#' @param verbose Logical; if `TRUE` (default), prints out the model's accuracy, pseudo R-squared,
#'   and coefficients.
#' @param save.acc Logical; if `TRUE`, saves the model's accuracy as a global variable `Acc.print`
#'   for further use. Default is `FALSE`.
#'
#' @return A list containing:
#' \item{class.table}{The confusion matrix of actual vs. predicted classes.}
#' \item{Accuracy}{The computed accuracy of the model as a percentage.}
#' \item{McFadden_R2}{The McFadden pseudo R-squared value of the model.}
#'
#' @details
#' The function generates predictions for the entire dataset and compares them with the true class
#' labels. It calculates the confusion matrix, overall accuracy, and McFadden's pseudo R-squared
#' for evaluating the model's goodness of fit. For multinomial regression models (`nnet::multinom`
#' or `MASS::polr`), the log-likelihood of the null model (intercept-only) is compared to that of
#' the fitted model to compute the pseudo R-squared.
#'
#' If `verbose` is set to `TRUE`, the function outputs:
#' \itemize{
#'   \item The confusion matrix.
#'   \item The accuracy of the model.
#'   \item McFadden's pseudo R-squared value.
#'   \item The model coefficients.
#' }
#' @export
mod.info <- function(model, data, verbose = TRUE, save.acc = FALSE) {

  # Generate predictions using the fitted model on the full dataset
  pred <- predict(model, newdata = data, 'class')

  # Extract the actual class labels from the dataset
  actual <- data$class

  # Create a confusion matrix (classification table) comparing actual vs. predicted classes
  class.table <<- table(actual, pred)

  # Optionally save accuracy as a global variable if save.acc is TRUE
  if (save.acc == TRUE) {
    Acc.print <<- round((sum(diag(class.table)) / sum(class.table)) * 100, 2)
  }

  # Compute the overall accuracy based on the classification table
  Accuracy <- paste(round((sum(diag(class.table)) / sum(class.table)) * 100, 2), "%", sep = '')

  # Null model (intercept-only model) for comparison to calculate McFadden's R-squared
  if (length(grep('MASS::polr', test$call)) > 0) {
    test.0 <- MASS::polr(class ~ 1, data = data, Hess = TRUE, control = list(maxit = 100))
  } else {
    test.0 <- nnet::multinom(class ~ 1, data = data, maxit = 2000, trace = F)
  }

  # Log-likelihood of the fitted model and null model
  loglik.0 <- as.numeric(nnet:::logLik.multinom(test.0))  # Log-likelihood of the null model
  loglik.1 <- as.numeric(nnet:::logLik.multinom(model))   # Log-likelihood of the fitted model

  # Compute McFadden's pseudo R-squared
  McFadden_R2 <- round((1 - (loglik.1 / loglik.0)), digits = 3)

  # Create a table of accuracy and McFadden's R-squared
  st <- cbind(Accuracy, McFadden_R2)
  names(st) <- c("Accuracy", "McFadden_R2")

  # Display the model statistics (Accuracy and Pseudo-R2) as a knitr table
  stats <- knitr::kable(st, caption = "\n\nFull Model Stats - Overall Accuracy and Pseudo-R2")

  # Display the model coefficients as a knitr table
  ce <- knitr::kable(coef(model), caption = "\n\nModel Coefficients")

  # If verbose is TRUE, print the model statistics and coefficients
  if (verbose == TRUE) {
    print(stats)
    print(ce)
  }
}



#' Confusion Matrix Heatmap Plot
#'
#' This function generates a heatmap visualization of a confusion matrix using `ggplot2`.
#' It calculates proportions for True Positives, False Positives, Total Size, Precision,
#' and Accuracy, and presents them as part of the plot. The plot includes both counts
#' and percentages of correct/incorrect classifications, along with other summary metrics.
#'
#' @param class.table A confusion matrix in table or matrix format.
#' @param plot.title The title to be displayed on the plot.
#' @param conformation A subtitle or caption to be displayed on the plot (e.g., model information).
#'
#' @return A heatmap plot visualizing the confusion matrix with counts, percentages, and classification metrics.
#'
#' @details
#' - The function takes a confusion matrix and computes the total row and column sums.
#' - It calculates proportions for recall, precision, accuracy, and sizes, and uses these metrics to create
#'   a heatmap. True classifications are marked as 'True' and incorrect ones as 'False'.
#' - Additional metrics like overall accuracy, class size, and precision are calculated and displayed.
#' - The heatmap tiles are colored based on the classification outcome, and the alpha level (transparency)
#'   of the tiles represents the computed proportions for better visual distinction.
#' - The plot uses `ggplot2` with customized aesthetics, and the plot title and caption can be personalized.
#' - Class-wise recall values are also computed and saved as a global variable `Low.Recall` for further use.
#' @export
#' @import ggplot2
#' @import ggrepel
ct.plot <- function(class.table, plot.title, conformation) {

  # Convert the confusion matrix to a matrix format
  ct <- as.matrix(class.table)
  classes <- nrow(ct)  # Number of classes

  # Initialize a vector for row totals and append to the confusion matrix
  total <- rep(NA, nrow(ct))
  ct <- cbind(ct, total)

  # Calculate total counts for each row
  for (i in 1:dim(ct)[1]) {
    ct[i, (dim(ct)[1] + 1)] <- sum(ct[i, 1:(dim(ct)[2] - 1)])
  }

  # Initialize a vector for column totals and append to the confusion matrix
  total <- rep(NA, ncol(ct))
  ct <- rbind(ct, total)

  # Calculate total counts for each column
  for (i in 1:dim(ct)[2]) {
    ct[dim(ct)[2], i] <- sum(ct[1:(dim(ct)[1] - 1), i])
  }

  # Reshape the confusion matrix for ggplot
  ct <- reshape2::melt(ct)
  names(ct) <- c('Exp', 'Pred', 'Freq')  # Rename columns for clarity
  ct$Exp <- as.factor(ct$Exp)  # Convert to factors for ordered plotting
  ct$Pred <- as.factor(ct$Pred)

  # Initialize columns for legend and proportion
  ct <- dplyr::mutate(ct, Legend = rep(NA, nrow(ct)))
  ct <- dplyr::mutate(ct, prop = rep(NA, nrow(ct)))

  # Calculate proportions and legends for each class
  for (i in as.numeric(row.names(ct[ct$Exp != 'total' & ct$Pred != 'total', ]))) {
    for (j in 1:classes) {
      if (ct$Exp[[i]] == as.character(j)) {
        ct$prop[[i]] <- ct$Freq[[i]] / sum(ct$Freq[ct$Exp == as.character(j) & ct$Pred != 'total'])
      }
    }
    ct$Legend[i] <- ifelse(ct$Exp[i] == ct$Pred[i], 'True', 'False')  # Classify as True or False
  }

  # Round and adjust proportions for display
  ct$prop <- round(ct$prop, digits = 3) * 100

  # Handle total predictions for size and precision
  for (i in as.numeric(row.names(ct[ct$Exp == 'total' | ct$Pred == 'total', ]))) {
    if (ct$Pred[i] == 'total' & ct$Exp[i] != 'total') {
      ct$prop[i] <- (ct$Freq[i] / ct$Freq[ct$Exp == 'total' & ct$Pred == 'total']) * 100
      ct$Legend[i] <- 'Size'
    }
    if (ct$Pred[i] != 'total' & ct$Exp[i] == 'total') {
      ct$prop[[i]] <- (ct$Freq[ct$Exp == ct$Pred[i] & ct$Pred == ct$Pred[i]] / ct$Freq[i]) * 100
      ct$Legend[i] <- 'Precision'
    }
    if (ct$Pred[i] == 'total' & ct$Exp[i] == 'total') {
      ct$prop[[i]] <- (sum(diag(class.table)) / ct$Freq[i]) * 100
      ct$Legend[i] <- 'Accuracy'
    }
  }

  ct$prop <- round(ct$prop, digits = 1)  # Final rounding of proportions

  # Create a column for plotting
  ct <- dplyr::mutate(ct, value = ct$prop)

  # Set proportions for non-true/false entries to a constant value for visualization
  for (i in 1:nrow(ct)) {
    if (ct$Legend[i] != 'True' & ct$Legend[i] != 'False') {
      ct$prop[i] <- 45  # Constant value for other metrics
    }
  }

  ct[is.na(ct)] <- 0  # Replace NA values with 0

  # Create the heatmap using ggplot2
  base <- ggplot(data = ct,
                 mapping = aes(x = ordered(Pred, levels = sort(unique(Pred))),
                               y = ordered(Exp, levels = rev(sort(unique(Exp)))),
                               fill = Legend,
                               alpha = prop)) +
    geom_tile(color = 'black', size = 1.5) +  # Tile borders
    coord_fixed() +
    geom_text(aes(label = paste(Freq, "\n", '(', value, '%', ')', sep = '')),
              size = 4, vjust = .5, fontface = "bold", alpha = 1) +  # Text labels on tiles
    scale_fill_manual(values = c(True = '#6CAE75', False = '#8A0526',
                                 Size = '#FCDEBE', Precision = '#DCAB6B',
                                 Accuracy = '#247BA0')) +  # Color mapping
    scale_alpha(guide = "none") +  # Hide prop legend
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 9, face = 'bold'),
          axis.text.y = element_text(size = 9, face = 'bold'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom")  # Move legend to bottom

  # Extract recall for true classifications
  true.recall <- data.frame(ct$prop[ct$Exp == ct$Pred & ct$Exp != 'total'])
  row.names(true.recall) <- unique(as.character(ct$Exp[ct$Exp != 'total']))
  names(true.recall) <- 'Class Recall'
  Low.Recall <<- true.recall  # Save true recall as a global variable

  # Add titles and captions to the plot
  base + labs(title = plot.title,
              subtitle = 'Confusion Matrix',
              caption = conformation)
}

#' Probability Heatmap Plot
#'
#' This function generates a heatmap to visualize the predicted probabilities of class membership
#' for a classification model using `ggplot2`. It highlights correct and incorrect predictions
#' with different colors and displays the predicted probabilities for each class as part of the plot.
#'
#' @param model A fitted classification model used to make predictions.
#' @param data A dataset with the same structure as the training data, including the true class labels.
#' @param plot.title The title to be displayed on the plot.
#' @param conformation A subtitle or caption to be displayed on the plot (e.g., model version or information).
#'
#' @return A heatmap plot visualizing the predicted probabilities for each class, along with
#' indicators for correct and incorrect predictions.
#'
#' @details
#' - **Predictions and Probabilities**:
#'   - The function predicts the class membership for each instance in the dataset and also retrieves
#'     the predicted class probabilities from the model.
#'   - These probabilities are converted into percentages and visualized in the heatmap.
#'   - Each row represents an observation, and the columns represent the predicted probabilities for each class.
#' - **Color Coding**:
#'   - Correct predictions are colored green.
#'   - Misclassifications where the second-highest probability corresponds to the true class are colored orange.
#'   - Other misclassifications are colored red.
#' - **Heatmap**:
#'   - The heatmap shows the predicted probabilities using a gradient from white to blue.
#'   - Text labels show the exact probability values in each tile.
#'   - For each observation, the true class is indicated by a separate tile and labeled with its class number.
#' - **Customization**:
#'   - The title and caption of the plot can be customized via the function parameters.
#' @export
#' @import ggplot2
prob.heatmap <- function(model, data, plot.title, conformation) {
  # Determine number of classes
  classes <- length(model$lev)

  # Generate predictions and probabilities
  pred <- predict(model, newdata = data, 'class')
  r.w <- pred == data$class
  probs <- predict(model, newdata = data, 'probs') * 100  # Convert to percentage

  # Adjust probabilities if binary classification
  if (dim(data.frame(probs))[2] == 1) {
    probs <- cbind(100 - data.frame(probs), data.frame(probs))
  }

  # Create data frame with expected and predicted classes, and correct/wrong classification
  verif <- data.frame(cbind(data$class, pred, r.w, probs, rep(NA, nrow(probs))))
  row.names(verif) <- row.names(probs)
  colnames(verif)[c(1, dim(verif)[2])] <- c('Exp', 'color')  # Rename columns

  # Assign colors based on correct classification, second-highest probability, or incorrect prediction
  for (i in 1:dim(verif)[1]) {
    second <- sort(as.numeric(verif[i, 4:(4 + classes)]), decreasing = T)[2]  # Get 2nd highest probability
    where.is.sec <- which(as.numeric(verif[i, 4:(4 + (classes - 1))]) == second)
    if (verif$r.w[i] == 1) {
      verif$color[i] <- "#66a180"  # Green for correct predictions
    } else {
      if (as.numeric(verif$Exp[i]) == where.is.sec) {
        verif$color[i] <- 'tan1'  # Tan if second-highest probability matches actual class
      } else {
        verif$color[i] <- '#d1505b'  # Red for incorrect prediction
      }
    }
  }

  # Prepare data for plotting
  pro.df <- data.frame(probs)
  pro.df <- tibble::rownames_to_column(pro.df)
  pro.df[, 1] <- factor(pro.df[, 1], levels = pro.df[, 1])  # Order row names as factors
  pro.df[, (dim(pro.df)[2] + 1)] <- as.numeric(data$class)
  colnames(pro.df) <- c('Name', model$lev, 'Exp')  # Rename columns for plotting
  row.names(pro.df) <- row.names(probs)

  # Convert data to long format for ggplot
  long <- reshape2::melt(pro.df, id.vars = 'Name')
  long[,3] <- round(long[,3], digits = 0)  # Round probabilities to integers for display
  long <- dplyr::mutate(long, exp_shape = rep(NA, nrow(long)))  # Add column for expected shapes

  # Assign expected shapes based on class
  for (i in 1:nrow(long)) {
    for (j in 1:classes) {
      if (long$variable[i] == 'Exp' & long$value[i] == j) {
        long$exp_shape[i] <- j
      }
    }
  }

  # Prepare color vector for labels
  col_vec <- vector(mode = 'numeric')
  coloring <- c("darkgoldenrod4", "slateblue", 'darksalmon',
                'darkblue', 'navajowhite4',
                'darkcyan', 'chocolate4',
                "coral3","cornsilk4",'darkslateblue')  # Define colors for each class

  for (i in 1:length(long[long$variable == 'Exp',4])) {
    col_vec[i] <- coloring[long[long$variable == 'Exp',4][i]]
  }

  # Prepare labels for text display
  label.vec <- as.character(long[long$variable == 'Exp',4])

  # Generate the probability heatmap
  prob.heatmap <- ggplot2::ggplot(mapping = ggplot2::aes(x = variable,
                                                         y = ordered(Name,
                                                                     levels = rev(factor(pro.df$Name,
                                                                                         levels = pro.df$Name))))) +
    ggplot2::geom_tile(data = long[long$variable != 'Exp',],
                       color = 'black', ggplot2::aes(fill = value)) +  # Add tiles for probabilities
    ggplot2::coord_fixed(ratio = 0.5) +  # Adjust cell height with ratio
    ggplot2::geom_text(data = long[long$variable != 'Exp',],
                       ggplot2::aes(label = value)) +  # Display probability values in each cell
    ggplot2::scale_fill_gradient(name = "% probability",
                                 low = "#FFFFFF",
                                 high = "dodgerblue2",
                                 guide = ggplot2::guide_colorbar(frame.colour = "black",
                                                                 ticks.colour = "black")) +  # Set fill gradient for probabilities
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, face = 'bold'),
                   axis.text.y = ggplot2::element_text(size = 10, face = 'bold',
                                                       colour = rev(verif$color)),  # Color y-axis labels by verification result
                   axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "bottom") +  # Place legend at bottom of plot
    ggplot2::scale_x_discrete(position = "top", limits = levels(long$variable)) +  # Display x-axis labels at top
    ggplot2::geom_tile(data = long[long$variable == 'Exp',],
                       alpha = 0, inherit.aes = F,
                       ggplot2::aes(x = rev(variable),
                                    y = ordered(Name, levels = rev(factor(pro.df$Name,
                                                                          levels = pro.df$Name))))) +  # Add empty tiles for 'Exp' rows
    ggplot2::geom_text(data = long[long$variable == 'Exp',], label = label.vec,
                       size = 4, color = col_vec, fontface = 'bold')  # Display 'Exp' labels in bold

  # Add titles and caption to the plot
  prob.heatmap + labs(title = plot.title,
                      subtitle = 'Probability Heatmap',
                      caption = conformation)
}


#' Similarity-Based Sampling Function
#'
#' This function samples a subset of observations from a dataset based on their similarity
#' to class-specific centroids (average class feature vectors). It calculates the similarity
#' of each observation to its own class as well as to other classes and returns the indices of
#' the sampled subset. It can optionally plot the similarity before and after sampling.
#'
#' @param data A data frame containing features and class labels. The class label column should be named 'class'.
#' @param class The class from which to sample observations based on similarity.
#' @param compare.with The index of the class similarity to compare (default is 0, which refers to the "same class" similarity).
#' @param plot Logical. If `TRUE`, plots of similarity before and after group truncation will be displayed. Default is `FALSE`.
#' @param sample.size The number of observations to sample from the specified class. Default is the size of the smallest class in the data.
#'
#' @return A vector of indices corresponding to the rows in the dataset that have been sampled.
#'
#' @details
#' - **Data Scaling**:
#'   - The feature columns are scaled to have mean zero and unit variance before similarity calculations.
#'   - The columns `flag` and `tag` (if present) are excluded from the analysis.
#' - **Similarity Calculation**:
#'   - For each class, a centroid vector is computed by averaging the scaled feature vectors of observations in that class.
#'   - The similarity between each observation and its class's centroid vector is calculated using the cosine similarity formula.
#' - **Sampling**:
#'   - The function samples observations from the specified class by evenly selecting points across the similarity range.
#'   - The sample size is specified by `sample.size`, and the sampled points are those closest to evenly spaced intervals within the range of similarities.
#' - **Plotting**:
#'   - If `plot = TRUE`, the function will generate two plots using `ggplot2`:
#'     - **Before Truncation**: A plot showing the similarity of each instance to its class centroid before sampling.
#'     - **After Truncation**: A plot showing the similarity of the sampled instances after truncation.
#'   - Both plots will use colors to represent different classes.
#' @export
#' @import ggplot2
simi.sampler <- function(data, class,
                         compare.with = 0,
                         plot = FALSE,
                         sample.size = min(summary(as.factor(data$class)))) {

  # Identify the output column and relevant variables
  out.col <- which(colnames(data) == 'class')
  vars <- names(data[,-out.col])
  vars <- vars[vars != 'flag' & vars != 'tag']  # Exclude 'flag' and 'tag' columns

  # Prepare data for sampling by scaling numeric values
  sampler.data <- data[vars]
  sampler.data <- data.frame(apply(sampler.data, 2, as.numeric))
  sampler.data <- data.frame(scale(sampler.data, TRUE, TRUE))  # Scale data

  classes <- length(unique(data$class))  # Count unique classes

  # Compute mean values and magnitudes for each class
  for (i in 1:classes) {
    assign(paste('class.', i, '.vector', sep = ''),
           apply(sampler.data[data$class == i, , drop = FALSE], 2, mean))
    assign(paste('class.', i, '.mag', sep = ''),
           sqrt(sum(apply(sampler.data[data$class == i, , drop = FALSE], 2, mean)^2)))
  }

  # Compute similarity of each instance with its class
  new.col <- ncol(sampler.data) + 1
  for (r in 1:nrow(sampler.data)) {
    for (i in 1:classes) {
      if (data$class[r] == i) {
        vec <- get(ls()[grepl(paste('.', i, '.vector', sep = ''), ls())])
        mag <- get(ls()[grepl(paste('.', i, '.mag', sep = ''), ls())])
        sampler.data[r, new.col] <- sum(vec * sampler.data[r, 1:(new.col - 1)]) /
          (mag * sqrt(sum(sampler.data[r, 1:(new.col - 1)]^2)))
      }
    }
  }

  # Compute similarity between instances and all classes
  simi.df <- data.frame(matrix(ncol = classes, nrow = nrow(data)))
  for (i in 1:nrow(simi.df)) {
    for (j in 1:classes) {
      vec <- get(ls()[grepl(paste('.', j, '.vector', sep = ''), ls())])
      mag <- get(ls()[grepl(paste('.', j, '.mag', sep = ''), ls())])
      simi.df[i, j] <- sum(vec * sampler.data[i, 1:(new.col - 1)]) /
        (mag * sqrt(sum(sampler.data[i, 1:(new.col - 1)]^2)))
    }
  }
  names(simi.df) <- as.character(seq(1, classes))  # Name columns by class numbers

  # Prepare similarity table
  names(sampler.data)[new.col] <- 'class.similarity'
  simi.table <- cbind(sampler.data[new.col], simi.df)
  simi.table <- dplyr::mutate(simi.table, class = data$class)
  simi.table <- dplyr::mutate(simi.table, Name = row.names(data))
  simi.table <- dplyr::mutate(simi.table, flag = data$flag)
  colnames(simi.table)[c(1, (2 + classes):ncol(simi.table))] <- c('same.class', 'class', 'Name', 'flag')

  # Plot similarity before truncation
  plot.sim.before <- ggplot(simi.table, aes(simi.table[, compare.with + 1], y = class, label = Name)) +
    geom_point(aes(color = class)) +
    geom_text_repel(aes(label = Name), max.overlaps = 25) +
    theme(axis.line = element_line(size = 1, colour = "black"),
          axis.text.x = element_text(colour = "black", size = 12, face = 'bold'),
          axis.text.y = element_text(colour = "black", size = 12, face = 'bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          legend.position = c(2, 2)) +
    xlab(names(simi.table)[compare.with + 1]) +
    ggtitle('Similarity Before Group Truncation')

  # Sampling based on similarity
  simi.class <- simi.table[simi.table$class == class, compare.with + 1]
  steps <- sort(seq(min(simi.class), max(simi.class), (max(simi.class) - min(simi.class)) / (sample.size - 1)))

  dis.mat <- matrix(ncol = length(steps), nrow = length(sort(simi.class)))
  for (i in 1:nrow(dis.mat)) {
    dis.mat[i, ] <- abs(simi.class[i] - steps)
  }

  pts <- vector()
  row.names(dis.mat) <- as.character(simi.table$flag[simi.table$class == class])

  # Determine which points to keep
  if (length(steps) < length(simi.class)) {
    for (i in 1:ncol(dis.mat)) {
      drop <- which.min(dis.mat[, i])
      pts[i] <- simi.table[as.numeric(names(drop)), 1]
      dis.mat <- dis.mat[-drop, ]
    }
  } else {
    pts <- simi.table[as.numeric(row.names(dis.mat)), 1]
  }

  keep <- as.numeric(simi.table$flag[simi.table$same.class %in% pts])

  # Determine rows to truncate
  class.rows <- simi.table$flag[simi.table$class == class]
  if (min(class.rows) != 1) {
    truncated <- unique(c(1:(min(class.rows) + 1), keep, (max(class.rows) + 1):nrow(data)))
  } else {
    truncated <- unique(c(keep, (max(class.rows) + 1):nrow(data)))
  }

  # Prepare data for after truncation plot
  simi.plot.data <- simi.table[truncated, ][complete.cases(simi.table[truncated, ]), ]
  plot.sim.after <- ggplot(simi.plot.data,
                           aes(x = simi.plot.data[, compare.with + 1], y = class, label = Name)) +
    geom_point(aes(color = class)) +
    geom_text_repel(aes(label = Name), max.overlaps = 25) +
    theme(axis.line = element_line(size = 1, colour = "black"),
          axis.text.x = element_text(colour = "black", size = 12, face = 'bold'),
          axis.text.y = element_text(colour = "black", size = 12, face = 'bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          legend.position = c(2, 2)) +
    xlab(names(simi.table)[compare.with + 1]) +
    ggtitle('Similarity After Group Truncation')

  # Display plots if requested
  if (isTRUE(plot)) gridExtra::grid.arrange(plot.sim.before, plot.sim.after, ncol = 2)

  return(keep)  # Return the indices of samples to keep
}

#' Fit a Proportional Odds Logistic Regression Model
#'
#' @description
#' `fit_polr_model` fits a proportional odds logistic regression model using the `MASS::polr` function.
#' The function attempts to fit the model with a provided formula, data, and starting coefficient values,
#' and includes error handling to ensure robustness.
#'
#' @param formula A formula representing the model structure. This should be in the standard R formula notation.
#' @param data A data frame containing the variables in the model.
#' @param start_vector A numeric vector of starting values for the model coefficients.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{success}{Logical, indicating whether the model fit was successful (TRUE) or not (FALSE).}
#'   \item{model}{If `success` is TRUE, the fitted `polr` model object. If `success` is FALSE, this will be `NULL`.}
#'   \item{error}{If `success` is FALSE, the error message produced during the fitting process.}
#' }
#'
#' @details
#' The function internally calls `MASS::polr` to fit the model. It includes error handling using `tryCatch`
#' to capture and return any errors that occur during model fitting.
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   data(iris)
#'   iris$Species <- as.factor(iris$Species)
#'   formula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
#'   start_vector <- c(0.1, 0.1, 0.1, 0.1)
#'   fit <- fit_polr_model(formula, iris, start_vector)
#'   if (fit$success) {
#'     print(summary(fit$model))
#'   } else {
#'     print(fit$error)
#'   }
#' }
#'
#' @seealso
#' \code{\link[MASS]{polr}}
#'
#' @export
fit_polr_model <- function(formula, data, start_vector) {
  tryCatch({
    model <- MASS::polr(formula, data = data, Hess = TRUE, start = start_vector, control = list(maxit = 100))
    return(list(success = TRUE, model = model))
  }, error = function(e) {
    return(list(success = FALSE, error = e))
  })
}

#' Fit a Proportional Odds Logistic Regression (POLR) Model
#'
#' This function iteratively fits a proportional odds logistic regression (POLR) model using the
#' provided formula and dataset. If the initial model fails to converge, the function adjusts the
#' starting coefficient vector and retries fitting the model until it succeeds.
#'
#' @param formula A formula specifying the model structure (e.g., `class ~ var1 + var2 + ...`).
#' @param data A data frame containing the variables in the formula and the response variable.
#' @return A fitted `polr` model object.
#'
#' @details
#' This function calls \code{fit_polr_model()} to fit a \code{MASS::polr()} model. If the model
#' fails to converge, it expands the starting coefficient vector and retries fitting the model
#' in a while loop until successful convergence is achieved.
#'
#' @examples
#' # Example usage:
#' \dontrun{
#' data(example_data)
#' fit_polr(class ~ var1 + var2, data = example_data)
#' }
#'
#' @seealso \code{\link[MASS]{polr}}, \code{\link{fit_polr_model}}
#' @export
fit_polr <- function(formula, data) {
  num.of.vars <- stringi::stri_count(formula, fixed = '+')
  start <- c(rep(0, num.of.vars + 2), 1)

  # Flag to control the loop
  success <- FALSE

  # While loop to keep trying until the model fits successfully
  while (!success) {
    result <- fit_polr_model(formula, data, start)

    if (result$success) {
      success <- TRUE
      test <- result$model
    } else {
      # Modify the start vector by adding another 0
      start <- c(0, start)
    }
  }
  return(test)
}

#' Clean Correlated Features Based on Multiple Criteria
#'
#' This function cleans up highly correlated features in a dataset using a user-defined threshold and a selection criterion.
#' Users can choose to retain the feature with the higher correlation to the outcome, higher mutual information,
#' higher feature importance from a random forest, or higher variance.
#'
#' @param data A data frame containing the features and outcome.
#' @param outcome_col A character string specifying the name of the outcome column.
#' @param corr_threshold A numeric value specifying the threshold for identifying highly correlated features. Defaults to 0.9.
#' @param method A character string specifying the method to choose between correlated features. Options include:
#'   \itemize{
#'     \item \code{"correlation"}: Select based on higher absolute correlation with the outcome.
#'     \item \code{"mutual_information"}: Select based on higher mutual information with the outcome.
#'     \item \code{"feature_importance"}: Select based on higher feature importance from a random forest model.
#'     \item \code{"variance"}: Select based on higher variance in the feature values.
#'   }
#'   Defaults to \code{"correlation"}.
#' @param outcome_type A character string specifying the type of the outcome. Options are:
#'   \itemize{
#'     \item \code{"numeric"}: Use if the outcome is continuous.
#'     \item \code{"categorical"}: Use if the outcome is categorical.
#'     \item \code{"auto"}: Automatically detect the type of the outcome. Defaults to \code{"auto"}.
#'   }
#' @return A cleaned data frame where correlated features have been removed based on the chosen criterion.
#'
#' @details The function computes the correlation matrix for the input features and identifies pairs of features
#' that exceed the specified correlation threshold. For each pair of correlated features, the function selects
#' one feature to retain based on the selected \code{method}.
#'
#' - When using \code{"correlation"}, the function retains the feature with the higher absolute correlation to the outcome.
#' - When using \code{"mutual_information"}, mutual information is computed to select the feature that provides more information about the outcome.
#' - The \code{"feature_importance"} method fits a random forest model to calculate feature importance and selects the more important feature.
#' - The \code{"variance"} method retains the feature with the higher variance.
#'
#' The function works with both numeric and categorical outcome variables. It automatically adjusts calculations
#' based on the type of the outcome, which can be inferred by the function or specified by the user.
#'
#' @examples
#' # Example 1: Using correlation to clean features (numeric outcome)
#' data(iris)
#' cleaned_data <- clean_correlated_features(data = iris, outcome_col = "Sepal.Length",
#'                                           corr_threshold = 0.8, method = "correlation")
#'
#' # Example 2: Using mutual information to clean features (categorical outcome)
#' iris$Species <- as.factor(iris$Species)
#' cleaned_data <- clean_correlated_features(data = iris, outcome_col = "Species",
#'                                           corr_threshold = 0.8, method = "mutual_information", outcome_type = "categorical")
#'
#' # Example 3: Using feature importance from a random forest (numeric outcome)
#' cleaned_data <- clean_correlated_features(data = iris, outcome_col = "Sepal.Length",
#'                                           corr_threshold = 0.8, method = "feature_importance")
#'
#' @export
#' @importFrom randomForest randomForest importance
#' @importFrom infotheo mutinformation discretize
#' @importFrom stats cor var cor.test
clean_correlated_features <- function(data, outcome_col = dim(data)[2], corr_threshold = 0.9, method = "correlation", outcome_type = "auto") {

  # Extract the outcome column
  outcome <- data[[outcome_col]]

  # Determine outcome type if not provided explicitly
  if (outcome_type == "auto") {
    outcome_type <- if (is.factor(outcome) || is.character(outcome)) "categorical" else "numeric"
  }

  # Extract features, excluding the outcome column
  features <- data[, -outcome_col]

  # Compute correlation matrix for features
  corr_matrix <- cor(features)

  # Identify pairs of highly correlated features
  high_corr_pairs <- which(abs(corr_matrix) > corr_threshold, arr.ind = TRUE)
  high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] != high_corr_pairs[, 2], ]

  # Keep track of features to drop
  drop_features <- c()

  # Iterate over each correlated pair
  for (i in 1:nrow(high_corr_pairs)) {
    feat1 <- colnames(corr_matrix)[high_corr_pairs[i, 1]]
    feat2 <- colnames(corr_matrix)[high_corr_pairs[i, 2]]

    # Skip if either feature is already marked for dropping
    if (feat1 %in% drop_features || feat2 %in% drop_features) next

    if (method == "correlation") {
      # Calculate correlation with the outcome (numeric or categorical)
      if (outcome_type == "numeric") {
        corr_feat1 <- abs(cor(features[, feat1], outcome))
        corr_feat2 <- abs(cor(features[, feat2], outcome))
      } else {
        # For categorical outcome, use point-biserial correlation or a similar method
        corr_feat1 <- abs(cor.test(as.numeric(outcome), features[, feat1])$estimate)
        corr_feat2 <- abs(cor.test(as.numeric(outcome), features[, feat2])$estimate)
      }

    } else if (method == "mutual_information") {
      # Calculate mutual information with the outcome
      if (outcome_type == "numeric") {
        mi_feat1 <- mutinformation(discretize(features[[feat1]]), discretize(outcome))
        mi_feat2 <- mutinformation(discretize(features[[feat2]]), discretize(outcome))
      } else {
        # Categorical outcome - direct mutual information
        mi_feat1 <- mutinformation(discretize(features[[feat1]]), as.factor(outcome))
        mi_feat2 <- mutinformation(discretize(features[[feat2]]), as.factor(outcome))
      }

      corr_feat1 <- mi_feat1
      corr_feat2 <- mi_feat2

    } else if (method == "feature_importance") {
      # Calculate feature importance using Random Forest
      if (outcome_type == "numeric") {
        model <- randomForest(x = features, y = outcome)
      } else {
        model <- randomForest(x = features, y = as.factor(outcome))
      }

      importance <- importance(model)
      corr_feat1 <- importance[row.names(importance) == feat1, ]
      corr_feat2 <- importance[row.names(importance) == feat2, ]

    } else if (method == "variance") {
      # Compare variance
      corr_feat1 <- var(features[[feat1]])
      corr_feat2 <- var(features[[feat2]])
    }

    # Drop the feature with the lower score
    if (corr_feat1 < corr_feat2) {
      drop_features <- c(drop_features, feat1)
    } else {
      drop_features <- c(drop_features, feat2)
    }
  }

  # Return the cleaned dataset without the dropped features
  return(data[, setdiff(names(data), drop_features), drop = FALSE])
}

