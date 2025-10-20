"
# Install
remotes::install_github('barkais/rxn.cond.class', force = TRUE)
"

# Packege loading
library('rxn.cond.class')

# Load data
ArylBorylation <- read.csv('Stevens_data_organized_for_classification.csv')
rownames(ArylBorylation) <- ArylBorylation[,162]
ArylBorylation <- ArylBorylation[,-c(27,28,159:164)] # Removes names and smiles
ArylBorylation <- ArylBorylation[,-c(116,117)] #removel of all zero column

# Convert the first 158 columns to numeric values
ArylBorylation[,c(1:154)] <- as.numeric(as.matrix(ArylBorylation[,c(1:154)]))

# data set as factor
ArylBorylation[,155] <- as.factor(as.numeric(ArylBorylation[,155]))

# Add a 'flag' column to sequentially number the rows
ArylBorylation <- plyr::mutate(ArylBorylation, flag = seq(1,nrow(ArylBorylation)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(ArylBorylation, 1)  # Sample from class 1
two <- simi.sampler(ArylBorylation, 2)  # Sample from class 2
three <- simi.sampler(ArylBorylation, 3)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(ArylBorylation, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(ArylBorylation, 2, 3)  

# Combine the similarities from the various classes into one vector
"
similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three)
                 "
similarties <- c(one,two,three)
# Define train and test data from the samples taken from groups 1 and 2
Train.data <- ArylBorylation[similarties,]
Test.data <- ArylBorylation[-similarties,]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Train.data, 
                        min = 1, 
                        max = 1, 
                        ordinal = F)
knitr::kable(models)


# Training set ------------------------------------------------------------

# Use the first ranked non-ordinal model
test.form <- models[1, 1]

# Train the non-ordinal multinomial regression model
test <- nnet::multinom(test.form,
                       data = Train.data,
                       maxit = 2000, 
                       trace = FALSE)

# Cross-validation (smallest-group's-fold)
k.fold.log.iter(formula = test.form, 
                data = Train.data, 
                ordinal = FALSE, 
                stratify = TRUE, 
                iterations = 20, 
                verbose = TRUE)

# Leave-one-out cross-validation
k.fold.log.iter(formula = test.form, 
                data = Train.data, 
                ordinal = FALSE, 
                folds = nrow(Train.data), 
                stratify = FALSE, 
                iterations = 1, 
                verbose = TRUE)


# Visualization Training --------------------------------------------------

# Display model information and confusion matrix plot
model.info <- mod.info(test, Train.data, TRUE, TRUE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Training Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')


# Test Set ----------------------------------------------------------------

# Evaluate the model on the test set
model.info <- mod.info(test, Test.data, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')

