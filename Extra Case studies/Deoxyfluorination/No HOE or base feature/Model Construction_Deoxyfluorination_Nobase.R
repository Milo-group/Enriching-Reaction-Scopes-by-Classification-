"
# Install
remotes::install_github('barkais/rxn.cond.class', force = TRUE)
"

# Packege loading
library('rxn.cond.class')

# =======================================================================
# Cl Sulfonyl -----------------------------------------------------------
# =======================================================================

# Load data
ClSulfonyl <- read.csv('Sulfonyl_3-Cl_Final.csv')
rownames(ClSulfonyl) <- ClSulfonyl[,2]
ClSulfonyl <- ClSulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
ClSulfonyl[,c(17:79)] <- as.numeric(as.matrix(ClSulfonyl[,c(17:79)]))

# data set as factor
ClSulfonyl[, c(1:16, 80)] <- lapply(ClSulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
ClSulfonyl <- plyr::mutate(ClSulfonyl, flag = seq(1,nrow(ClSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(ClSulfonyl, 1)  # Sample from class 1
two <- simi.sampler(ClSulfonyl, 2)  # Sample from class 2
three <- simi.sampler(ClSulfonyl, 3) # Sample from class 3 
four <- simi.sampler(ClSulfonyl,4) # Sample from class 4

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(ClSulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(ClSulfonyl, 2, 3) 

# Sample class 4 molecules similar to class 3
four_three <- simi.sampler(ClSulfonyl, 4, 3) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_three), 
                 union(two, two_three),
                 three,
                 union(four,four_three))
                 
# Define train and test data from the samples taken from groups 1 and 2
Train.data <- ClSulfonyl[similarties,]
Test.data <- ClSulfonyl[-similarties,]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Train.data, 
                        min = 3, 
                        max = 3, 
                        ordinal = F)
knitr::kable(models)


# Training set ------------------------------------------------------------

# Use the first ranked non-ordinal model
test.form <- models[2, 1]

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
                            conformation = '2. 2nd Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '2. 2nd Place')


# Test Set ----------------------------------------------------------------

# Evaluate the model on the test set
model.info <- mod.info(test, Test.data, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '2. 2nd Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '2. 2nd Place')



# =======================================================================
# Pyfluoro Sulfonyl -----------------------------------------------------
# =======================================================================

# Load data
PyfluoroSulfonyl <- read.csv('Sulfonyl_pyfluoro_Final.csv')
rownames(PyfluoroSulfonyl) <- PyfluoroSulfonyl[,2]
PyfluoroSulfonyl <- PyfluoroSulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
PyfluoroSulfonyl[,c(17:79)] <- as.numeric(as.matrix(PyfluoroSulfonyl[,c(17:79)]))

# data set as factor
PyfluoroSulfonyl[, c(1:16, 80)] <- lapply(PyfluoroSulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
PyfluoroSulfonyl <- plyr::mutate(PyfluoroSulfonyl, flag = seq(1,nrow(PyfluoroSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(PyfluoroSulfonyl, 1)  # Sample from class 1
two <- simi.sampler(PyfluoroSulfonyl, 2,sample.size = round(0.7*sum(PyfluoroSulfonyl$class ==2)))  # Sample from class 2
three <- simi.sampler(PyfluoroSulfonyl, 3)  # Sample from class 3
four <- simi.sampler(PyfluoroSulfonyl, 4)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(PyfluoroSulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(PyfluoroSulfonyl, 2, 3)

# Sample class 2 molecules similar to class 3
four_three <- simi.sampler(PyfluoroSulfonyl, 4, 3) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three,
                 union(four,four_three))
                 
                 
# Define train and test data from the samples taken from groups 1 and 2
Train.data <- PyfluoroSulfonyl[similarties,]
Test.data <- PyfluoroSulfonyl[-similarties,]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Train.data, 
                        min = 3, 
                        max = 3, 
                        ordinal = F)
knitr::kable(models)


# Training set ------------------------------------------------------------

# Use the first ranked non-ordinal model
test.form <- models[2, 1]

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
                            conformation = '2nd Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '2. 2nd Place')


# Test Set ----------------------------------------------------------------

# Evaluate the model on the test set
model.info <- mod.info(test, Test.data, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '2nd Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '2. 2nd Place')

# =======================================================================
# CF3 Sulfonyl ----------------------------------------------------------
# =======================================================================

# Load data
CF3Sulfonyl <- read.csv('Sulfonyl_3-CF3_Final.csv')
rownames(CF3Sulfonyl) <- CF3Sulfonyl[,2]
CF3Sulfonyl <- CF3Sulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
CF3Sulfonyl[,c(17:79)] <- as.numeric(as.matrix(CF3Sulfonyl[,c(17:79)]))

# data set as factor
CF3Sulfonyl[, c(1:16, 80)] <- lapply(CF3Sulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
CF3Sulfonyl <- plyr::mutate(CF3Sulfonyl, flag = seq(1,nrow(CF3Sulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(CF3Sulfonyl, 1)  # Sample from class 1
two <- simi.sampler(CF3Sulfonyl, 2)  # Sample from class 2
three <- simi.sampler(CF3Sulfonyl, 3)  # Sample from class 3
four <- simi.sampler(CF3Sulfonyl, 4)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(CF3Sulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(CF3Sulfonyl, 2, 3)

# Sample class 2 molecules similar to class 3
four_three <- simi.sampler(CF3Sulfonyl, 4, 3) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three,
                 union(four,four_three))


# Define train and test data from the samples taken from groups 1 and 2
Train.data <- CF3Sulfonyl[similarties,]
Test.data <- CF3Sulfonyl[-similarties,]
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

# =======================================================================
# NO2 Sulfonyl ----------------------------------------------------------
# =======================================================================

# Load data
NO2Sulfonyl <- read.csv('Sulfonyl_3-NO2_Final.csv')
rownames(NO2Sulfonyl) <- NO2Sulfonyl[,2]
NO2Sulfonyl <- NO2Sulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
NO2Sulfonyl[,c(17:79)] <- as.numeric(as.matrix(NO2Sulfonyl[,c(17:79)]))

# data set as factor
NO2Sulfonyl[, c(1:16, 80)] <- lapply(NO2Sulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
NO2Sulfonyl <- plyr::mutate(NO2Sulfonyl, flag = seq(1,nrow(NO2Sulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(NO2Sulfonyl, 1,sample.size = round(0.45*sum(NO2Sulfonyl$class == 1)))  # Sample from class 1
two <- simi.sampler(NO2Sulfonyl, 2)  # Sample from class 2
three <- simi.sampler(NO2Sulfonyl, 3)  # Sample from class 3
four <- simi.sampler(NO2Sulfonyl, 4)  # Sample from class 3

# Sample class 1 molecules similar to class 2
one_two <- simi.sampler(NO2Sulfonyl, 1, 2)  

# Sample class 2 molecules similar to class 2
three_two <- simi.sampler(NO2Sulfonyl, 3, 2)

# Sample class 2 molecules similar to class 2
four_two <- simi.sampler(NO2Sulfonyl, 4, 2) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_two), 
                 two, 
                 union(three,three_two),
                 union(four,four_two))


# Define train and test data from the samples taken from groups 1 and 2
Train.data <- NO2Sulfonyl[similarties,]
Test.data <- NO2Sulfonyl[-similarties,]

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

# =======================================================================
# PBSF Sulfonyl ---------------------------------------------------------
# =======================================================================

# Load data
PBSFSulfonyl <- read.csv('Sulfonyl_PBSF_Final.csv')
rownames(PBSFSulfonyl) <- PBSFSulfonyl[,2]
PBSFSulfonyl <- PBSFSulfonyl[,-c(1,2,81)] # Removes names and inchies

# Convert the first 158 columns to numeric values
PBSFSulfonyl[,c(16:78)] <- as.numeric(as.matrix(PBSFSulfonyl[,c(16:78)]))

# data set as factor
PBSFSulfonyl[, c(1:15, 79)] <- lapply(PBSFSulfonyl[, c(1:15, 79)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
PBSFSulfonyl <- plyr::mutate(PBSFSulfonyl, flag = seq(1,nrow(PBSFSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(PBSFSulfonyl, 1,sample.size = 0.5*sum(PBSFSulfonyl$class=="1"))  # Sample from class 1
two <- simi.sampler(PBSFSulfonyl, 2)  # Sample from class 2
three <- simi.sampler(PBSFSulfonyl, 3, sample.size = 4)  # Sample from class 3

# Sample class 1 molecules similar to class 2
one_two <- simi.sampler(PBSFSulfonyl, 1, 2)  

# Sample class 2 molecules similar to class 2
three_two <- simi.sampler(PBSFSulfonyl, 3, 2, sample.size = 4)

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_two), 
                 two, 
                 union(three,three_two))

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- PBSFSulfonyl[similarties,]
Test.data <- PBSFSulfonyl[-similarties,]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Train.data, 
                        min = 3, 
                        max = 3, 
                        ordinal = F)
knitr::kable(models)


# Training set ------------------------------------------------------------

# Use the first ranked non-ordinal model
test.form <- models[4, 1]

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
                            conformation = '4. 4th Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '4. 4th Place')


# Test Set ----------------------------------------------------------------

# Evaluate the model on the test set
model.info <- mod.info(test, Test.data, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '4. 4th Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '4. 4th Place')

