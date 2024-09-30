```r
library(rxn.cond.class)

# Load data from a CSV file into a data frame
data <- data.frame(data.table::fread('Training_Data.csv'), check.names = F)

# Clean and organize data
row.names(data) <- data[,2]  # Set the second column as row names
data$class <- as.factor(data$class)  # Convert the 'class' column to a factor
data <- data[,-c(1:2)]  # Remove the first and second columns (name and tag)

# Similarity-based sampling for each class
one <- simi.sampler(data, 1)  # Sample from class 1
two <- simi.sampler(data, 2)  # Sample from class 2
three <- simi.sampler(data, 3)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(data, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(data, 2, 3)  

# Combine the similarities from the various classes into one vector
similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three)

# Train ordinal models using the subset of data corresponding to similarities
models.ordinal <- sub_model_log(data = data[similarties, ], 
                                min = 4, 
                                max = 4, 
                                ordinal = T)


# Define the training set from the subset of similarities
Train.data <- data[similarties, ]

# Define the test set from the remaining rows not in the similarities set
Test.data <- data[-similarties, ]
```
