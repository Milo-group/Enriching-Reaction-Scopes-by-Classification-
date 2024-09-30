```r
library(rxn.cond.class)

# Import training data from a CSV file
data_bin <- data.frame(data.table::fread('Training_Data_rdkit_stereoelectronic.csv'), check.names = F)

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

# Sample from group 1 (first 50 columns + flag) with sample size being 75% of group 1 size
one <- simi.sampler(data_bin[, c(1:50, 65)], 1, sample.size = round(sum(data_bin$class == 1) * 0.75))

# Sample from group 2 (first 50 columns + flag) with sample size being 75% of group 2 size
two <- simi.sampler(data_bin[, c(1:50, 65)], 2, sample.size = round(sum(data_bin$class == 2) * 0.75))

# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = data_bin[c(one, two), ], 
                                 min = 4, 
                                 max = 4, 
                                 ordinal = F)

# Present the list of models in a table format
knitr::kable(models)

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- data_bin[c(one, two), ]
Test.data <- data_bin[-c(one, two), ]

```
