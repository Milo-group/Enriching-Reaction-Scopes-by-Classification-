```r
library(rxn.cond.class)

# Import training data from a CSV file
data_bin <- data.frame(data.table::fread("Training_Data_rdkit_stereoelectronic.csv"), check.names = F) 

# Set the first column as row names
row.names(data_bin) <- data_bin[,1]  
data_bin <- data_bin[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
data_bin$class <- as.factor(data_bin$class)  

# Convert the first 50 columns to numeric values
data_bin[,c(1:50)] <- as.numeric(as.matrix(data_bin[,c(1:50)]))  

# Convert columns 51 to 64 to factors after converting to numeric
for (i in 51:64) data_bin[,i] <- as.factor(as.numeric(data_bin[,i]))

# Add a 'flag' column to sequentially number the rows
data_bin <- plyr::mutate(data_bin, flag = seq(1,nrow(data_bin)))

# Sample from group 1 using simi.sampler
one <- simi.sampler(data_bin, 1)

# Sample 21 molecules from group 2
two <- simi.sampler(data_bin, 2, sample.size = 21)

# Sample from group 3 using simi.sampler
three <- simi.sampler(data_bin, 3)

# Combine the similarities from groups 1, 2, and 3 into one vector
uni_similarties <- c(one, 
                     two, 
                     three)

# Train models on the subset of data corresponding to the combined similarities
models <- sub_model_log(data = data_bin[uni_similarties, ], 
                        min = 3, 
                        max = 3, 
                        ordinal = T)

# Present the list of models in a table format
knitr::kable(models)

# Define training data from the subset of similarities
Train.data <- data_bin[uni_similarties, ]

# Define test data from the remaining rows not in the similarities set
Test.data <- data_bin[-uni_similarties, ]
```
