# Function to evaluate the reproducabilty of mutational ordering 
# in the presence of co-occuring mutations and non-independent mutational paths.

# RMO function adpated from Toprak et al 2011

# Geno Guerra
# University of California, Berkeley
# Created on May 11, 2019


# Set location of the ordered paths file (this will also be where the output files will be written).
#setwd("~/Dropbox/Mutational_Ordering")
setwd("~/Desktop/")
# Read in the paths file
Paths = read.table("RMO_Paths.txt", header = TRUE, sep =  "")

#Set the number of independent iterations
num_samples = 500

#Set the number of random permutations per iteration
n_randomizations = 100000


# Takes as input a string path, and creates a list of values by randomly ordering co-occuring mutations.
Parse_Path <- function(string_path){
  ordered_path = c()
  # First spilt the path by commas
  split_path = strsplit(string_path, ",")
  
  # For each of the entries in split path, check for "+", if none, assign to the ordered path,
  # if at least 1 "+", randomly assign to the ordered path. 
  # Assuming paths are only coded by numbers, "," and "+"
  for(p in c(1:length(split_path[[1]]))){
    
    # If we contain a non digit (implying there must be a "+"), randomly order the digits, and assign to ordered_path
    if(grepl('\\D',split_path[[1]][p]) ){
      split_sub_path = strsplit(split_path[[1]][p],split = "+",fixed = TRUE )
      ordered_path = c(ordered_path, sample(split_sub_path[[1]])) 
      #'sample' here is a random ordering of the mutations which co-occurred.
    }else{
      # There contains no "+" character, so the ordering is clear. 
      ordered_path = c(ordered_path, split_path[[1]][p])
    }
  }
  return(as.integer(ordered_path))  
  
}
  
  
  




# Reproducibility of mutational order (RMO) score function (adapted from Toprak et al 2011)
RMO = function(List_Path,n_groups){
  # Takes in a list of paths
  
  score = 0 
  
  # This loops through paths p = (1: n_groups -1), and compares to (p+1: n_groups),
  # to ensure that each pair of paths is compared exactly once. 
  
  #Note: These nested for-loops are a death wish for run time.
  
  for(path in c(1:(n_groups-1))){
    # Get the current path 
    cur_path = List_Path[[path]]
    len_path = length(cur_path)
    for(i in c(1:(len_path-1))){
      for(j in c((i+1):len_path)){
        # This loops through all pairs of mutations in the path, in the order they appear (so p1 appears before p2)
        p1 = cur_path[i]
        p2 = cur_path[j]
        
        
        # Loop through all other paths, comparing this specific pair ordering with the pair ordering in each path.
        
        for(z in c((path+1):n_groups)){
            
          # if both mutations are in the next list (p1 and p2 in List_Path[[z]]), 
          if( (p1 %in% List_Path[[z]])&(p2 %in% List_Path[[z]]) ){
            #   check if they occur in the same order, or the opposite order. 
            if(which(List_Path[[z]] ==p1) < which(List_Path[[z]]== p2)){
              #if they occur in the same order, and are both present, score +1
            score = score + 1  
            }else{
              # if they don't occur in the same order, but are both present, score -1
              score = score - 1
            } # close else
          } # close if both in next path, z.
        }# close loop through paths after path 'path'
   
      } # close loop through p2 (j)
    } # close loop through p1 (i)
  
  } # close loop through paths. 
  
  return(score)
}

#------------ Set up the ability to hold the p_values for each random sample, and set the number of samples.--------

##pvals = mat.or.vec(num_samples,1)
# I don't use this (above) anymore as I instead write each p-value to output as soon as it is computed.


# ---Do the sampling num_samples times. 
for(s in c(1:num_samples)){

#------------ Randomly select one path from each Group. -----------

# Get the unique list of groups, and the number of groups. 
groups = unique(Paths$Group)
n_groups = length(groups)


# This holds the subset of paths where one path is randomly sampled from each group. 
Subset_Paths = data.frame(n_groups,2)
colnames(Subset_Paths) = c('Groups', 'Paths')

# Do the sub sampling. 
for(g in c(1:n_groups)){
  
  
  # Get all of the paths in this group
  indivs = which(Paths$Group ==groups[g])
  
  # if only one path, take it 
  if(length(indivs) == 1){
    indiv = indivs[1]
  }else{
    # otherwise, randomly sample 1 path from the list of paths in this group. 
      indiv = sample(indivs,1)
    }
  
  # Assign their Path row to the Subset matrix row. 
  Subset_Paths[g,] <- c(toString(Paths[indiv,1]), toString(Paths[indiv,2]))

}

#---------------Parse the path into individual characters, breaking ties by random ordering.---------------

# We use the 'list' data-type as the paths are of varying length, and this data structure allows for that. 

Path_List = list()


for(g in c(1:n_groups)){
  
  # Take the path, which is currently of 'string' type, and parse it into a vector of integers, randomly ordering any 'co-occuring' mutations (mutations that we have been unable to order)
  Path_List[[g]] = Parse_Path(Subset_Paths[g,2])
  
}
#---------------------------------Get score of true values---------------------------------

True_Score = RMO(Path_List, n_groups)

# ----------------------------Now randomize a set amount of iterations and get new scores. ----------------------------


# This list will hold all of the randomization scores
Random_Scores = mat.or.vec(n_randomizations,1)

for( r in c(1:n_randomizations)){
  Random_List =list()
  for(l in c(1:n_groups)){
      # This randomizes the mutational path for sampled individual, l.
    Random_List[[l]] = sample(Path_List[[l]])
  }
  
  Random_Scores[r] = RMO(Random_List, n_groups)
}

# ------ Store the pvalue of the true_score in a file names 'p_values.txt'------------
#  Calculate the pvalue
#pvals[s] = length(which(Random_Scores > True_Score))/n_randomizations
pval = length(which(Random_Scores > True_Score))/n_randomizations
#  Store it
write(pval, file  = 'p_values.txt', append = TRUE)
write(pval, file  = 'paths_taken.txt', append = TRUE)
lapply(Path_List, write, "paths_taken.txt", append = TRUE, ncolumns = 10)

}

# -------------------Command to plot the pvals in a histogram.------------------------
##hist(pvals, xlab = "p-value", ylab = "frequency", main = "Distribution of p-values", breaks = 30)


