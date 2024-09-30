```r
library(rxn.cond.class)  # Load custom package for reaction condition classification

# Import training data from CSV and convert it into a data frame
data50 <- data.frame(data.table::fread('Training_Data.csv')) 
row.names(data50) <- data50[,1]  # Set the first column as row names
data50 <- data50[,-1]  # Remove the first column after setting row names

# Convert all columns in the data frame to factors
for(i in 1:ncol(data50)) data50[, i] <- as.factor(data50[, i])

# Add a new column 'flag' that sequentially numbers the rows
data50 <- plyr::mutate(data50, flag = seq(1,nrow(data50)))

## Sampler data

# Import sampling data from CSV and convert to data frame
sampl.dat <- data.table::fread('Sampling_Data.csv') 
sampl.dat <- data.frame(sampl.dat)
row.names(sampl.dat) <- sampl.dat[,1]  # Set the first column as row names
sampl.dat <- sampl.dat[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
sampl.dat$class <- as.factor(sampl.dat$class)

# Convert the first 60 columns of sampling data to numeric values
sampl.dat[,c(1:60)] <- as.numeric(as.matrix(sampl.dat[,c(1:60)]))

# Add a 'flag' column to sequentially number the rows
sampl.dat <- plyr::mutate(sampl.dat, flag = seq(1,nrow(sampl.dat)))

# Sample from group 1 using simi.sampler, returning samples from the smallest group
one <- simi.sampler(sampl.dat, 1)
two <- simi.sampler(sampl.dat, 2)
three <- simi.sampler(sampl.dat, 3)
four <- simi.sampler(sampl.dat, 4)

# Sample from group 2 molecules that are similar to group 1
two_three <- simi.sampler(sampl.dat, 2, 3)

# Sample from group 1 molecules that are similar to group 3
one_three <- simi.sampler(sampl.dat, 1, 3)

# Sample from group 4 molecules that are similar to group 3
four_three <- simi.sampler(sampl.dat, 4, 3)

## Back to model search

# Combine similarities from the various groups sampled
similarties <- c(three,
                union(two, two_three),
                union(one, one_three),
                union(four, four_three))

# Train models using sub_model_log on the subset of data corresponding to similarities
models <- sub_model_log(data = data50[similarties,], 
                        min = 4, 
                        max = 4, 
                        ordinal = F)

# Present the list of models in a table format
knitr::kable(models)

# Define training data from the subset of similarities
Train.data <- data50[similarties,]

# Define test data from the remaining rows not in the similarities set
Test.data <- data50[-similarties,]

```
