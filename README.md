
# `rxn.cond.class` - Classify Chemical Reaction Conditions

`rxn.cond.class` is an R package designed to classify and visualize chemical reaction conditions using both ordinal and non-ordinal models. It includes functionality for similarity-based sampling, model ranking, model evaluation, and heatmap visualization for model performance.

## Installation

To install the development version from GitHub, use:

```r
devtools::install_github("your-username/rxn.cond.class")
```

## Package Overview

`rxn.cond.class` provides several functions for classification, probability analysis, data sampling, and model evaluation with a focus on reaction condition classification.

### Key Features

- **Similarity-based sampling**: Select representative samples based on class-specific feature vectors.
- **Model Ranking**: Compare ordinal and non-ordinal models based on performance.
- **Model Evaluation**: Evaluate models with cross-validation, accuracy, and McFadden's pseudo-R².
- **Heatmap Visualization**: Visualize prediction probabilities and classification tables.

## Example Usage

### Model Search

This section demonstrates how to clean data, perform similarity-based sampling, rank models, and split data into training and testing sets.

```r
# Load data
data <- rxn.cond.class::example_training_data

# Clean and organize data
row.names(data) <- data[,2]
data$class <- as.factor(data$class)
data <- data[,-c(1:2)] # Remove name and tag

# Perform similarity-based sampling
one <- simi.sampler(data, 1) # Class 1 with itself
two <- simi.sampler(data, 2) # Class 2 with itself
three <- simi.sampler(data, 3) # Class 3 with itself
one_three <- simi.sampler(data, 1, 3) # Class 1 with Class 3
two_three <- simi.sampler(data, 2, 3) # Class 2 with Class 3

# Combine similarities
similarties <- c(union(one, one_three), union(two, two_three), three)

# Rank ordinal models
models.ordinal <- sub_model_log(data = data[similarties, ],
                                min = 2,
                                max = 2, 
                                ordinal = TRUE)

# Rank non-ordinal models
models.non.ordinal <- sub_model_log(data = data[similarties, ],
                                    min = 2,
                                    max = 2, 
                                    ordinal = FALSE)

# Define training and test sets
Train.set <- data[similarties, ]
Test.set <- data[-similarties, ]

# Load and organize external validation data
Prediction.set <- rxn.cond.class::example_validation_data
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[,-1]
Prediction.set$class <- as.factor(Prediction.set$class)
row.names(Prediction.set) <- RN
```

### Ordinal Model Example

#### Model Ranking

```r
# Present the ranked list of ordinal models
knitr::kable(models.ordinal)
```

#### Training Set

```r
# Use the first ranked ordinal model
test.form <- models.ordinal[1, 1]

# Define starting coefficients
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)

# Train the ordinal logistic regression model
test <- MASS::polr(test.form,
                   data = Train.set,
                   Hess = TRUE, 
                   control = list(maxit = 100), 
                   start = start)

# Cross-validation (smallest-group's-fold)
k.fold.log.iter(formula = test.form, 
                data = Train.set, 
                ordinal = TRUE, 
                stratify = TRUE, 
                iterations = 20, 
                verbose = TRUE)

# Leave-one-out cross-validation
k.fold.log.iter(formula = test.form, 
                data = Train.set, 
                ordinal = TRUE, 
                folds = nrow(Train.set), 
                stratify = FALSE, 
                iterations = 1, 
                verbose = TRUE)
```

#### Model Information and Visualization (Training Set)

```r
# Display model information and confusion matrix plot
mod.info(test, Train.set, TRUE, TRUE)

# Classification table plot
ct.plot(class.table, 
        plot.title = 'Training Set', 
        conformation = '1. 1st Place')

# Prediction probability heatmap
prob.heatmap(test, Train.set, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')
```

#### Test Set

```r
# Evaluate the model on the test set
mod.info(test, Test.set, FALSE, FALSE)

# Classification table plot
ct.plot(class.table, 
        plot.title = 'Test Set', 
        conformation = '1. 1st Place')

# Prediction probability heatmap
prob.heatmap(test, Test.set, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')
```

#### External Validation

```r
# Evaluate the model on the external validation set
mod.info(test, Prediction.set, FALSE)

# Classification table plot
ct.plot(class.table, 
        plot.title = 'External Validation', 
        conformation = '1. 1st Place')

# Prediction probability heatmap
prob.heatmap(test, Prediction.set, 
             plot.title = 'External Validation', 
             conformation = '1. 1st Place')
```

### Non-ordinal Model Example

#### Model Ranking

```r
# Present the ranked list of non-ordinal models
knitr::kable(models.non.ordinal)
```

#### Training Set

```r
# Use the first ranked non-ordinal model
test.form <- models.non.ordinal[1, 1]

# Train the non-ordinal multinomial regression model
test <- nnet::multinom(test.form,
                       data = Train.set,
                       maxit = 2000, 
                       trace = FALSE)

# Cross-validation (smallest-group's-fold)
k.fold.log.iter(formula = test.form, 
                data = Train.set, 
                ordinal = FALSE, 
                stratify = TRUE, 
                iterations = 20, 
                verbose = TRUE)

# Leave-one-out cross-validation
k.fold.log.iter(formula = test.form, 
                data = Train.set, 
                ordinal = FALSE, 
                folds = nrow(Train.set), 
                stratify = FALSE, 
                iterations = 1, 
                verbose = TRUE)
```

#### Model Information and Visualization (Training Set)

```r
# Display model information and confusion matrix plot
mod.info(test, Train.set, TRUE, TRUE)

# Classification table plot
ct.plot(class.table, 
        plot.title = 'Training Set', 
        conformation = '1. 1st Place')

# Prediction probability heatmap
prob.heatmap(test, Train.set, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')
```

#### Test Set

```r
# Evaluate the model on the test set
mod.info(test, Test.set, FALSE, FALSE)

# Classification table plot
ct.plot(class.table, 
        plot.title = 'Test Set', 
        conformation = '1. 1st Place')

# Prediction probability heatmap
prob.heatmap(test, Test.set, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')
```

#### External Validation

```r
# Evaluate the model on the external validation set
mod.info(test, Prediction.set, FALSE)

# Classification table plot
ct.plot(class.table, 
        plot.title = 'External Validation', 
        conformation = '1. 1st Place')

# Prediction probability heatmap
prob.heatmap(test, Prediction.set, 
             plot.title = 'External Validation', 
             conformation = '1. 1st Place')
```

## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
