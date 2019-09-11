# ModifiedRMO

This set of R scripts implements the RMO (Reproducibiilty of Mutational Order) score first proposed in Toprak et al (2012) with the new additions of accounting for non-independent mutational pathways (in the presence of shared ancestry), and partially unresolved mutational pathways. 

RMO measures the level of structure in a set of paths, giving positive weight for pairs of mutations when they occur in the same order on two different paths, and negative weight when they occur in opposite orders. A high RMO value indicates high structure, and a low RMO value indicates random ordering. RMO as defined in Toprak et al (2012) assumes independent evolution of lineages, such that all mutations are independent of other lineages, and the ordering of all mutations can be resolved, completely. 

In the presence of a species phylogeny, mutations can occur on shared branches, and multiple mutations can occur on the same branch where their order cannot be determined. We here develop a method to make use of RMO in the presence of these two complications. 

The input data to Reproducability_Permutation_Test.R, RMO_Paths.txt is a comma-separated set of mutational pathways, grouped by paths that share ancestry (mutations which arise on the branch of a shared ancestor) where the paths need not be fully resolved (a "+" symbol signifies co-occurence on the same branch). 

Reproducability_Permutation_Test.R takes as input this file, and for each of a set number of replicates, randomly chooses one pathway from each group, arbitrarily resolves unresolved mutational orderings, and computes the RMO score of this subset. Next, for specified number of permutations, the paths within this subset are randomly permuted, and a new RMO score is calculated. We count the number of randomly permuted pathways that have equal or higher RMO scores to our unpermuted subset. A p-value can be calcuated which measures the likelihood of the unpermuted subset having arisen completely by chance (aka with no structure to the pathways).  

One output to this script is p_values.txt, which is the set of p-values as described above, one for each replicate described above. 

The second output to this script is paths_taken.txt , which includes the p-values and he unpermuted subset that was sampled in that replicate. This allows the user to study conditional distributions of p-values based on which paths were sampled/how mutational orderings were resolved. 


The R script p_value_Plotter.R takes as input p_values.txt and paths_taken.txt. The first part of the script produces a histogram of p-values for visualizing the results. The second part is specific to our data set, which includes some code for studying the conditional distributions of p-values for scenarios when order (2,3) is randomly chosen as the order in group A or when order (3,2) is randomly chosen. This allows us to understand the right tail distribution of the p-value histogram. 


## References
Toprak, Erdal, et al. "Evolutionary paths to antibiotic resistance under dynamically sustained drug selection." Nature genetics 44.1 (2012): 101.
