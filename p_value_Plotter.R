# Function to plot the results of the reproducabilty of mutational ordering 
# in the presence of co-occuring mutations and non-independent mutational paths.

#

# Geno Guerra
# University of California, Berkeley
# Created on May 11, 2019


# Set the location of the output from Reproducability_Permutation_Test.R
setwd("~/Desktop/")
#setwd("~/Dropbox/Mutational_Ordering")

# This plots the histogram of p-values.
p_vals = data.frame(read.table('p_values.txt', header = FALSE))
ps = sapply(p_vals, as.numeric)
hist(ps, breaks = 75, xlab = "p-value", col = "black", main = "Distribution of p-values")
length(which(ps < 0.05))/ nrow(ps)



# This evaluates the conditional distribution of p-values based on the occurence of a specific ordering (2,3) or (3,2)
no_col = max(count.fields('paths_taken.txt'), sep = " ")
Path_Data = data.frame(read.table('paths_taken.txt',sep =" ", header = FALSE, fill = TRUE, col.names = (c(1:no_col))))

counter = 1
n_paths = 13
n_iter = nrow(Path_Data)/(n_paths+1)
p_val_vec_32 = c()
p_val_vec_23 = c()
for(i in c(1:n_iter)){
  
  # determine where the 2 occurs in the first path
  l2 = which(Path_Data[counter +1, ] == 2)
  l3 = which(Path_Data[counter +1, ] == 3)
  
  if( l2 <l3){
    p_val_vec_23 = c(p_val_vec_23, Path_Data[counter,1])
    
  }
  if(l3 < l2){
    p_val_vec_32 = c(p_val_vec_32, Path_Data[counter,1])
    
  }
  counter = counter + n_paths + 1
  
}
