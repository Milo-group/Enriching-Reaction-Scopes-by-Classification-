
# Classify Chemical Reaction Conditions with the R Package `rxn.cond.class`

`rxn.cond.class` is an R package designed to classify and visualize logistic regression classification models for chemical reaction conditions using both ordinal and non-ordinal models. It includes functionality for similarity-based sampling, model ranking, model evaluation, and heatmap visualization for model performance.

## Installation

### Installation from Github 

First, install the `remotes` package from CRAN.
The `repos = getCRANmirrors()[1,"URL"]` addition helps installation in linux interactive sessions.

```r
install.packages('remotes', repos = getCRANmirrors()[1,"URL"])
```

Once `remotes` is properly installed, use the `install_github` function to install `rxn.con.class`.
For convenience, load the package.

```r
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')

# Load
library('rxn.con.class')
```

### Package overview and Information can be found on the pakcage [GitHub page](https://github.com/barkais/rxn.cond.class/tree/main). 

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
External.set <- rxn.cond.class::example_validation_data
RN <- External.set$V1
External.set <- External.set[,-1]
External.set$class <- as.factor(External.set$class)
row.names(External.set) <- RN

# Load and organize prediction of new substrates data
Prediction.set <- rxn.cond.class::example_prediction_data
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[,-1]
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

# Train model
test <- fit_polr(formula = test.form, data = Train.set)

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
model.info <- mod.info(test, Train.set, TRUE, TRUE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Training Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.set, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')
```

#### Test Set

```r
# Evaluate the model on the test set
model.info <- mod.info(test, Test.set, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.set, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')
```

#### External Validation

```r
# Evaluate the model on the external validation set
model.info <- mod.info(test, External.set, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'External Validation', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, External.set, 
             plot.title = 'External Validation', 
             conformation = '1. 1st Place')
```

#### Prediction of New Substartes

```r
knitr::kable(cbind(predict(test, Prediction.set, 'probs') * 100,
      predicted_class = predict(test, Prediction.set, 'class')))
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
model.info <- mod.info(test, Train.set, TRUE, TRUE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Training Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.set, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')
```

#### Test Set

```r
# Evaluate the model on the test set
model.info <- mod.info(test, Test.set, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.set, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')
```

#### External Validation

```r
# Evaluate the model on the external validation set
model.info <- mod.info(test, External.set, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'External Validation', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, External.set, 
             plot.title = 'External Validation', 
             conformation = '1. 1st Place')
```

## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
