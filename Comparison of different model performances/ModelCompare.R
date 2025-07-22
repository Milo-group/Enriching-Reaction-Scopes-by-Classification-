"
# Install
remotes::install_github('barkais/rxn.cond.class', force = TRUE)
install.packages('caret')
install.packages('e1071')
install.packages('nnet')
"


# Packege loading
library('rxn.cond.class')
library(caret)
library(e1071)
library(nnet)

# Data cleaning and test-train division using similarity sampler ----------

# Load data
data <- read.csv('Training_Data.csv')

# Clean and organize data
row.names(data) <- data[,2] # Tag rows by name
data$class <- as.factor(data$class) 
data <- data[,-c(1:2)] # Remove name and tag

# Perform similarity-based sampling
one <- simi.sampler(data, 1) # Similarity of class 1 with itself
two <- simi.sampler(data, 2) # Similarity of class 2 with itself
three <- simi.sampler(data, 3) # Similarity of class 3 with itself
one_three <- simi.sampler(data, 1, 3) # Similarity of class 1 with class 3
two_three <- simi.sampler(data, 2, 3) # Similarity of class 2 with class 3

# Combine similarities
similarties <- c(union(one, one_three), union(two, two_three), three)

# Define training and test sets
Train <- data[similarties, ]
Test <- data[-similarties, ]

# Load and organize external validation data
External <- rxn.cond.class::example_validation_data
RN <- External$V1
External <- External[,-1]
External$class <- as.factor(External$class)
row.names(External) <- RN
colnames(External) <- colnames(Train[, c(2:dim(Train)[2])])

# Scaling of all data sets
scale.z <- preProcess( Train[,-c(1,dim(Train)[2])], method = c("center", "scale"))
Train.set <- predict(scale.z, Train[,-c(1,dim(Train)[2])])
Test.set <- predict(scale.z, Test[,-c(1,dim(Test)[2])])
External.set <- predict(scale.z, External[,-c(dim(External)[2])])

Train.set <-data.frame(class=Train$class,Train.set)
Test.set <-data.frame(class=Test$class,Test.set)
External.set <-data.frame(class=External$class,External.set)

# NNet (neural network) ---------------------------------------------------

set.seed(1000)

nnet_scan=tune.nnet(class ~., data = Train.set, size = seq(20, 40, 1), 
                decay = seq(0, 1, 0.01),
                linout = F, trace = F, maxit = 1000, 
                tunecontrol = tune.control(cross = 10, best.model = T))

# Tunning size and scan features of the nnet function with 10-fold cross validation.

bestSize = nnet_scan$best.parameters[1,1]
bestDecay = nnet_scan$best.parameters[1,2]

# Results under seed 1000 gave bestSize = 39; bestDecay = 0.71
# Train a final model on the entire training data after tuning

set.seed(1000)
data <- as.data.frame(Train.set)
test <- as.data.frame(Test.set)
nnet_model <- nnet(class ~., data = Train.set, size = bestSize, decay = bestDecay,
              linout = F, trace = F, maxit = 10000, MaxNWts = 9999)

# Training accuracy
train_pred_nn <- predict(nnet_model, Train.set, type = "class")
cm_nn_train <- confusionMatrix(as.factor(train_pred_nn), Train.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_nn_train$table));cat("\n");print(cm_nn_train$overall["Accuracy"])

# Test accuracy
test_pred_nn <- predict(nnet_model, Test.set, type = "class")
test_pred_nn <- factor(test_pred_nn, levels = levels(Test.set$class)) # Matches factor levels
cm_nn_test <- confusionMatrix(test_pred_nn, Test.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_nn_test$table));cat("\n");print(cm_nn_test$overall["Accuracy"])


# External validation accuracy
external_pred_nn <- predict(nnet_model, External.set, type = "class")
cm_nn_external <- confusionMatrix(as.factor(external_pred_nn), External.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_nn_external$table));cat("\n");print(cm_nn_external$overall["Accuracy"])

# SVM (support vector machine) --------------------------------------------

set.seed(1000)

svm_model <- svm(class~., data = Train.set, kernel='linear')

# Train confusion matrix:
train_pred_svm <- predict(svm_model)
cm_svm_train <- confusionMatrix(train_pred_svm, Train.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_svm_train$table));cat("\n");print(cm_svm_train$overall["Accuracy"])

# Test confusion matrix:
test_pred_svm <- predict(svm_model, newdata = Test.set)
cm_svm_test <- confusionMatrix(test_pred_svm, Test.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_svm_test$table));cat("\n");print(cm_svm_test$overall["Accuracy"])

# Predict on the External data
external_pred_svm <- predict(svm_model, newdata = External.set)
cm_svm_external <- confusionMatrix(external_pred_svm, External.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_svm_external$table));cat("\n");print(cm_svm_external$overall["Accuracy"])


# Random Forest -----------------------------------------------------------

set.seed(1000)

train.control <- trainControl(method = "cv", number = 10)
rf <- caret::train(class ~., data=Train[-1],
                     method='rf',verbose=F,
                     trControl=train.control, tuneGrid = expand.grid(mtry = c(2, 4, 6, 8, 10)))

# tunegrid manually sets a grid of values for the 'mtry' parameter to be tested.
# mtry: Number of randomly selected predictors considered at each split in a tree.

# Best amount of predictors
print(rf$bestTune)

# accuracy vs. selected predictors
plot(rf)

# Accuracy of training set
train_pred_rf <- predict(rf)
cm_rf_train <- confusionMatrix(train_pred_rf, Train$class)
cat("Confusion Matrix\n\n"); print(t(cm_rf_train$table));cat("\n");print(cm_rf_train$overall["Accuracy"])

# Predict on the test data
test_pred_rf <- predict(rf, newdata = Test[-1])
cm_rf_test <- confusionMatrix(test_pred_rf, Test$class)
cat("Confusion Matrix\n\n"); print(t(cm_rf_test$table));cat("\n");print(cm_rf_test$overall["Accuracy"])

# Predict on the External data
external_pred_rf <- predict(rf, newdata = External)
cm_rf_external <- confusionMatrix(external_pred_rf, External$class)
cat("Confusion Matrix\n\n"); print(t(cm_rf_external$table));cat("\n");print(cm_rf_external$overall["Accuracy"])

# KNN (K- nearest neighbors) ----------------------------------------------

set.seed(1000)

train.control <- caret::trainControl(method = "cv", number = 10)
knn <- caret::train(class ~., data=Train.set,
                   method='knn',
                   trControl=train.control,tuneGrid = expand.grid(k = c(1, 3, 5, 7, 9, 11)))

# tunegrid manually sets a grid of values for the 'mtry' parameter to be tested.
# mtry: Number of randomly selected predictors considered at each split in a tree.

# Best amount of predictors
print(knn$bestTune)

# accuracy vs. k
plot(knn)

# training accuracy
train_pred_knn <- predict(knn)
cm_knn_train <- confusionMatrix(train_pred_knn, Train.set$class)
cat("Confusion Matrix\n\n"); print(cm_knn_train$table);cat("\n");print(cm_knn_train$overall["Accuracy"])

# Predict on the test data
test_pred_knn <- predict(knn, newdata = Test.set)
cm_knn_test <- confusionMatrix(test_pred_knn, Test.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_knn_test$table));cat("\n");print(cm_knn_test$overall["Accuracy"])

# Predict on the External data
external_pred_knn <- predict(knn, newdata = External.set)
cm_knn_external <- confusionMatrix(external_pred_knn, External.set$class)
cat("Confusion Matrix\n\n"); print(t(cm_knn_external$table));cat("\n");print(cm_knn_external$overall["Accuracy"])
