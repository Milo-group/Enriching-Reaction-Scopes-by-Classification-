"""
# Install
remotes::install_github('barkais/rxn.cond.class', force = TRUE)
"""

# Packege loading
library('rxn.cond.class')

# Load data
Sarpong2025 <- read.csv('computed_descriptors_Sarpong2025.csv')
rownames(Sarpong2025) <- Sarpong2025[,2]
Sarpong2025 <- Sarpong2025[,-c(1,2,3)]

# Convert the first 124 columns to numeric values
Sarpong2025[,c(1:124)] <- as.numeric(as.matrix(Sarpong2025[,c(1:124)]))

# data set as factor
for (i in 125:dim(Sarpong2025)[2]) Sarpong2025[,i] <- as.factor(as.numeric(Sarpong2025[,i]))

# Add a 'flag' column to sequentially number the rows
Sarpong2025 <- plyr::mutate(Sarpong2025, flag = seq(1,nrow(Sarpong2025)))

# Perform similarity-based sampling
# Sample from group 1 (first 124 columns + flag) with sample size being 75% of group 1 size
one <- simi.sampler(Sarpong2025[, c(1:124,(dim(Sarpong2025)[2]-1) ,dim(Sarpong2025)[2])], 1, sample.size = round(sum(Sarpong2025$class == 1) * 0.75))

# Sample from group 2 (first 124 columns + flag) with sample size being 75% of group 2 size
two <- simi.sampler(Sarpong2025[, c(1:124,(dim(Sarpong2025)[2]-1) ,dim(Sarpong2025)[2])], 2, sample.size = round(sum(Sarpong2025$class == 2) * 0.75))

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- Sarpong2025[c(one, two), 125:dim(Sarpong2025)[2]]
Test.data <- Sarpong2025[-c(one, two), 125:dim(Sarpong2025)[2]]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Train.data, 
                        min = 3, 
                        max = 3, 
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
